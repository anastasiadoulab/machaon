import random
import sys

import matplotlib.pyplot as plt
import numpy as np
from src.pdbhandler import PDBHandler
import time
import requests
import os
import json
import subprocess
import urllib.parse
import urllib.request
from tqdm import tqdm
from multiprocessing import Pool
from collections import defaultdict

class Enricher:

    def __init__(self):
        plt.style.use('ggplot')

        self.file_name = ''
        self.go_cache = defaultdict()

        self._go_request_url = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:'
        self._enriched_headers = ['pdbId', 'chainId', 'chainLength', 'resolution', 'geneId', 'taxonomy',
                                 'molecularFunction', 'cellularComponent', 'biologicalProcess']
        self._request_chunk_size = 100

        self.root_disk = ''
        self.pdb_dataset_path = ''
        self._uniprot_data_path = ''
        self._pdb_info_path = ''
        self._enrichment_cache = ''
        self._rcsb_enrichment_cache = ''
        self._entrez_cache = ''

        self._uniprot_request_url = 'https://www.uniprot.org/uniprot/ACCESSION.xml'
        self.cores = 1

        self.verbose = False
        self.debugging = False

    def set_root_disk(self, root_disk):
        # Set cache paths for Machaon
        self.root_disk = root_disk
        self._uniprot_data_path = os.path.sep.join([self.root_disk, 'uniprotdata'])
        self._pdb_info_path = os.path.join(self.root_disk, 'pdbinfo')
        self._enrichment_cache = os.path.join(self.root_disk, 'enrichment')
        self._rcsb_enrichment_cache = os.path.join(self.root_disk, 'rcsbenrich')
        self._entrez_cache = os.path.join(self.root_disk, 'entrez')
        if (os.path.exists(self._uniprot_data_path) is False):
            os.makedirs(self._uniprot_data_path)
        if (os.path.exists(self._pdb_info_path) is False):
            os.makedirs(self._pdb_info_path)
        if (os.path.exists(self._enrichment_cache) is False):
            os.makedirs(self._enrichment_cache)
        if (os.path.exists(self._rcsb_enrichment_cache) is False):
            os.makedirs(self._rcsb_enrichment_cache)
        if (os.path.exists(self._entrez_cache) is False):
            os.makedirs(self._entrez_cache)

    def load_go_cache(self):
        # Load Gene Ontology cache
        self.go_cache = defaultdict()
        if (os.path.exists('go_cache.csv')):
            with open('go_cache.csv', 'r', encoding='utf-8') as cacheFile:
                for line in cacheFile:
                    parts = line.strip('\n').split('\t')
                    self.go_cache[parts[0]] = (parts[1], parts[2])

    def update_go_cache(self, term_id, term_info):
        # Update Gene Ontology cache
        self.go_cache[term_id] = term_info
        with open('go_cache.csv', 'a', encoding='utf-8') as cacheFile:
            cacheFile.write(''.join([term_id, '\t', '\t'.join(term_info), '\n']))

    def get_uniprot_data(self, uniprot_accession):
        uniprot_file_path = os.path.join(self._uniprot_data_path, ''.join([uniprot_accession, '.xml']))
        if (os.path.exists(uniprot_file_path)):
            with open(uniprot_file_path, 'r') as uniprotFile:
                web_data = uniprotFile.read()
        else:
            # Retrieve data from UniProt with a variable delay
            # for requesting on different times in parallel threads.
            time.sleep(random.randint(5, 15))
            request_url = self._uniprot_request_url.replace('ACCESSION', uniprot_accession)
            print(request_url)
            result = requests.get(request_url)
            web_data = result.text
            with open(uniprot_file_path, 'w') as uniprotFile:
                uniprotFile.write(web_data)
        return web_data

    def access_pdb_info(self, entry):
        # Retrieve meta-data from a PDB file
        structure_id, chain_id = entry
        pdb_id = '_'.join([structure_id, chain_id])
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        cache_path = ''.join([self._pdb_info_path, os.path.sep, pdb_id, '.pdb'])
        pdb_dataset_path = ''.join([self.root_disk, os.path.sep, self.pdb_dataset_path])
        if (os.path.exists(cache_path) is False):
            residues, resolution, organism, gene_names = pdbhandler.get_pdb_info(
                ''.join([pdb_dataset_path, os.path.sep, structure_id, '.pdb']), chain_id)
            with open(cache_path, 'w') as dataFile:
                output = '\t'.join([repr(residues), repr(resolution), organism, gene_names])
                dataFile.write(output)
        else:
            with open(cache_path, 'r') as dataFile:
                residues, resolution, organism, gene_names = dataFile.readline().replace('None', '-1.0').split('\t')
        return pdb_id, residues, resolution, organism, gene_names

    def retrieve_uniprot_mapping(self, entry):
        # Retrieve PDB to UniProt ID mapping
        structure_id, chain_id = entry
        pdbhandler = PDBHandler(structure_id)
        pdbhandler.root_disk = self.root_disk
        pdbhandler.verbose = self.verbose
        pdb_path = ''.join([self.root_disk, os.path.sep, self.pdb_dataset_path, os.path.sep, structure_id, '.pdb'])
        return '_'.join([structure_id, chain_id]), pdbhandler.get_uniprot_accession_number(chain_id, structure_id,
                                                                                           pdb_path)

    def fetch_enrichment(self, data, column_headers=['b-phipsi', 'w-rdist', 't-alpha']):
        data_pdb_ids = data['structureId'].values
        data_values = data[column_headers].values
        total_pdb_ids = len(data_pdb_ids)
        enrichment = defaultdict()
        enrichmentinfo = defaultdict()
        pdb_info = defaultdict()
        entrez = defaultdict()
        ids_for_request = []
        ids_for_entrez_mapping = defaultdict()
        chunk_pdb_ids = []
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdbhandler.verbose = self.verbose

        enriched_filename = self.file_name.replace('.csv', '-enriched.csv')
        if os.path.exists(enriched_filename):
            os.remove(enriched_filename)

        # Use different strategies to gather information : cached data, local resources, online web services
        print(''.join(['Total PDB IDs: ', repr(total_pdb_ids)]))
        chunks = 0
        print("Fetching enrichment data from local sources...")
        if (self.debugging is False):
            with Pool(self.cores) as pool:
                for result in tqdm(pool.imap_unordered(self.read_enrichment, data_pdb_ids), total=len(data_pdb_ids), file=sys.stdout):
                    enrichment_data, pdb_details, pdb_id = result
                    if (len(enrichment_data) > 0):
                        entrez[pdb_id] = enrichment_data[0]
                        enrichment[pdb_id] = ','.join(enrichment_data)
                    pdb_info[pdb_id] = [pdb_details[0], pdb_details[1]]
                pool.close()
                pool.join()
        else:
            for chunked_index, chunked_pdb_id in enumerate(tqdm(data_pdb_ids)):
                enrichment_data, pdb_details, chunked_pdb_id = self.read_enrichment(chunked_pdb_id)
                if (len(enrichment_data) > 0):
                    entrez[chunked_pdb_id] = enrichment_data[0]
                    enrichment[chunked_pdb_id] = ','.join(enrichment_data)
                pdb_info[chunked_pdb_id] = [pdb_details[0], pdb_details[1]]

        for index, pdb_id in enumerate(tqdm(data_pdb_ids, file=sys.stdout)):
            if (pdb_id in enrichment):
                continue
            parts = pdb_id.split('_')
            structure_id, chain_id = '_'.join(parts[:-1]), parts[-1]
            if('AF-' not in pdb_id):
                ids_for_request.append(pdb_id)
                if (len(ids_for_request) == self._request_chunk_size or index + 1 == total_pdb_ids):
                    chunks += 1
                    print(''.join(['\nFetching GO entries [chunk:', repr(chunks), ']']))
                    ids_never_requested = []
                    for request_id in ids_for_request:
                        data_path = os.path.join(self._rcsb_enrichment_cache, ''.join([request_id, '.log']))
                        if (os.path.exists(data_path)):
                            with open(data_path, 'r') as dataFile:
                                enrichmentinfo[request_id] = dataFile.readline()
                        else:
                            ids_never_requested.append(request_id)
                    if (len(ids_never_requested) > 0):
                        requested_info = self.get_enrichment_info(ids_never_requested)
                        for id_index, requested_id in enumerate(ids_never_requested):
                            enrichmentinfo[requested_id] = requested_info[id_index]
                            data_path = os.path.join(self._rcsb_enrichment_cache, ''.join([requested_id, '.log']))
                            with open(data_path, 'w') as dataFile:
                                dataFile.write(requested_info[id_index])
                        ids_for_request = []
            accession_ids = []
            predicted_pdb_ids = []
            if (pdb_id not in entrez):
                if ('AF-' not in pdb_id):
                    chunk_pdb_ids.append([structure_id, chain_id])
                else:
                    # AlphaFold PDBs
                    predicted_pdb_ids.append([structure_id, chain_id])
                    pdbhandler.structure_id = structure_id
                    pdbhandler.get_uniprot_accession_by_alphafold_pdbid()
                    accession_ids.append(pdbhandler.uniprot_accession_number)
            # Retrieve Entrez information
            if (len(predicted_pdb_ids) == self._request_chunk_size or index + 1 == total_pdb_ids):
                print("Fetching Uniprot accession numbers...")
                entrez_mappings = self.map_uniprot_to_entrez(accession_ids)
                for keyIndex, key in enumerate(predicted_pdb_ids):
                    entrez[key] = entrez_mappings[keyIndex]
            if (len(chunk_pdb_ids) == self._request_chunk_size or index + 1 == total_pdb_ids):
                print("Fetching Uniprot accession numbers...")
                if (self.debugging is False):
                    with Pool(self.cores) as pool:
                        for result in tqdm(pool.imap_unordered(self.retrieve_uniprot_mapping, chunk_pdb_ids),
                                           total=len(chunk_pdb_ids), file=sys.stdout):
                            pdb_id_owner, uniprot_id = result
                            ids_for_entrez_mapping[pdb_id_owner] = uniprot_id
                        pool.close()
                        pool.join()
                else:
                    for chunked_index, chunked_pdb_id in enumerate(tqdm(chunk_pdb_ids)):
                        pdb_id_owner, uniprot_id = self.retrieve_uniprot_mapping(chunked_pdb_id)
                        ids_for_entrez_mapping[pdb_id_owner] = uniprot_id
                chunk_pdb_ids = list(ids_for_entrez_mapping.keys())
                print(''.join(['Requesting ENTREZ entries [chunk:', repr(chunks), ']']))
                accession_ids = []
                for key in chunk_pdb_ids:
                    accession_ids.append(ids_for_entrez_mapping[key])
                entrez_mappings = self.map_uniprot_to_entrez(accession_ids)
                for keyIndex, key in enumerate(chunk_pdb_ids):
                    entrez[key] = entrez_mappings[keyIndex]
                ids_for_entrez_mapping = defaultdict()
                chunk_pdb_ids = []

        # Store the enriched version
        with open(enriched_filename, 'w', encoding='utf-8') as output_file:
            output_line = []
            output_line.append('\t'.join(self._enriched_headers))
            output_line.append('\t')
            output_line.append('\t'.join(column_headers))
            output_line.append('\n')
            output_file.write(''.join(output_line))
            for chunked_index, chunked_pdb_id in enumerate(data_pdb_ids):
                output_line = []
                if (chunked_pdb_id in enrichmentinfo):
                    enrichment[chunked_pdb_id] = ','.join(
                        [entrez[chunked_pdb_id] if chunked_pdb_id in entrez else '#', enrichmentinfo[chunked_pdb_id]])
                    data_path = os.path.join(self._enrichment_cache, ''.join([chunked_pdb_id, '.log']))
                    with open(data_path, 'w') as dataFile:
                        dataFile.write(enrichment[chunked_pdb_id])
                parts = chunked_pdb_id.split('_')
                output_line.append('\t'.join(['_'.join(parts[:-1]), parts[-1]]))
                output_line.append('\t')
                output_line.append('\t'.join([repr(x) for x in pdb_info[chunked_pdb_id]]))
                output_line.append('\t')
                output_line.append(
                    enrichment[chunked_pdb_id].replace(',', '\t') if chunked_pdb_id in enrichment else '\t'.join(
                        ['#'] * 5))
                output_line.append('\t')
                output_line.append('\t'.join(
                    [repr(x) if isinstance(x, np.floating) or isinstance(x, float) else x for x in
                     data_values[chunked_index]]))
                output_line.append('\n')
                output_file.write(''.join(output_line))

    def read_enrichment(self, pdb_id):
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        parts = pdb_id.split('_')
        structure_id, chain_id = '-'.join(parts[:-1]), parts[-1]
        pdbpath = ''.join([pdb_dataset_path, os.path.sep, structure_id, '.pdb'])
        pdbhandler.structure_id = structure_id
        pdbpath = pdbhandler.handle_alphafold_pdbid(pdbpath)
        pdb_info = pdbhandler.get_pdb_info(pdbpath, chain_id)

        uniprot_dataset_path = os.path.sep.join([self.root_disk, 'idmapping_selected.tab.gz'])
        data_path = os.path.join(self._enrichment_cache, ''.join([pdb_id, '.log']))
        enrichment = []

        # Retrieve data for enriching an entry from cached or local resources
        if (os.path.exists(data_path)):
            with open(data_path, 'r') as dataFile:
                enrichment = dataFile.readline().split(',')
        else:
            try:
                pattern = pdb_id.replace('.', ':').replace('_', ':')
                # For AlphaFold PDB ID, use uniprot accession number instead
                if('AF-' in pdb_id):
                    pattern = pdbhandler.uniprot_accession_number
                output = subprocess.run(
                    ['zgrep', '-i', pattern, '-m', '1', uniprot_dataset_path],
                    capture_output=True,
                    timeout=30)
                parts = output.stdout.decode("utf-8").split('\t')

                # 1. UniProtKB-AC
                # 2. UniProtKB-ID
                # 3. GeneID (EntrezGene)
                # 4. RefSeq
                # 5. GI
                # 6. PDB
                # 7. GO
                # 8. UniRef100
                # 9. UniRef90
                # 10. UniRef50
                # 11. UniParc
                # 12. PIR
                # 13. NCBI-taxon
                # 14. MIM
                # 15. UniGene
                # 16. PubMed
                # 17. EMBL
                # 18. EMBL-CDS
                # 19. Ensembl
                # 20. Ensembl_TRS
                # 21. Ensembl_PRO
                # 22. Additional PubMed
                # NCBI taxon id 9606, homosapiens
                enrichment.append(parts[2] if len(parts[2]) > 0 else '#')
                enrichment.append(parts[12] if len(parts[12]) > 0 else '#')
                enrichment.append(','.join(self.split_go_annotations(parts[6].split(';'))))
            except:
                if (self.verbose is True):
                    print(''.join(['\n', pdb_id, ' not found in idmapping']))
                else:
                    pass
        return enrichment, pdb_info, pdb_id

    def map_uniprot_to_entrez(self, accession_ids):
        ids_for_request = []
        entrez_ids = defaultdict()
        results = []

        # Gather the UniProt IDs that have not been mapped before
        for accession_id_index, accession_id in enumerate(accession_ids):
            if (accession_id == ''):
                continue
            data_path = os.path.join(self._entrez_cache, ''.join([accession_id, '.log']))
            if (os.path.exists(data_path)):
                with open(data_path, 'r') as dataFile:
                    entrez_ids[accession_id] = dataFile.readline()
            else:
                ids_for_request.append(accession_id_index)

        # Retrieve Entrez IDs
        if (len(ids_for_request) > 0):
            accession_ids_for_request = [accession_ids[accessionIDIndex] for accessionIDIndex in ids_for_request]
            time.sleep(3)
            request_url = 'https://rest.uniprot.org/idmapping/run'
            job_polling_url = 'https://rest.uniprot.org/idmapping/status/'
            result_url = 'https://rest.uniprot.org/idmapping/uniref/results/'
            polling_timeout = 5
            entrez_ids = defaultdict()
            parameters = {
                'from': 'UniProtKB_AC-ID',
                'to': 'GeneID',
                'ids': ','.join(accession_ids_for_request)
            }
            # Request a job
            response = requests.post(request_url, params=parameters)
            if (response.status_code == 200):
                json_data = json.loads(response.text)
                if ('jobId' in json_data):
                    # Polling for a job
                    while True:
                        job_response = requests.get(''.join([job_polling_url, json_data['jobId']]))
                        if (job_response.status_code == 200):
                            job_data = json.loads(job_response.text)
                            result_data = []
                            if ('jobStatus' in job_data):
                                # Retrieve results from finished job
                                if (job_data['jobStatus'] == 'FINISHED'):
                                    result_response = requests.get(''.join([result_url, json_data['jobId']]))
                                    if (result_response.status_code == 200):
                                        result_data = json.loads(result_response.text)
                            elif ('results' in job_data):
                                result_data = job_data
                            if ('results' in result_data):
                                for entry in result_data['results']:
                                    entrez_ids[entry['from']] = entry['to']

                        time.sleep(3)
                        polling_timeout -= 1
                        if (polling_timeout == 0):
                            break

        # Cache the retrieved information
        for accession_id in accession_ids:
            entrez_id = entrez_ids[accession_id] if (accession_id in entrez_ids) else '#'
            data_path = os.path.join(self._entrez_cache, ''.join([accession_id, '.log']))
            if (os.path.exists(data_path) is False):
                with open(data_path, 'w') as dataFile:
                    dataFile.write(entrez_id)
            results.append(entrez_id)

        return results

    def split_go_annotations(self, terms):
        annotations = defaultdict()
        term_types = ['molecular_function', 'cellular_component', 'biological_process']
        for term_type in term_types:
            annotations[term_type] = []
        # Group the GO terms in three category types: 'molecular_function', 'cellular_component', 'biological_process'
        for term_id in terms:
            if ('GO:' not in term_id):
                continue
            term_id = term_id.replace('GO:', '').strip()
            if (term_id not in self.go_cache):
                term_info = self.get_go_term(term_id)
                self.update_go_cache(term_id, term_info)
            else:
                term_info = self.go_cache[term_id]
            if (term_info[1] in annotations):
                annotations[term_info[1]].append(term_id)
        data = []
        for term_type in term_types:
            data.append('#'.join(annotations[term_type]) if len(annotations[term_type]) > 0 else '#')
        return data

    def get_enrichment_info(self, structure_ids):
        # Use RCSB GraphQL services to retrieve information and enrich the selected entries
        time.sleep(3)
        results = []
        request_url = ''.join(['https://data.rcsb.org/graphql'])
        query = '''{
                                            polymer_entity_instances(instance_ids: [PDB_CHAIN])
                      {
                            polymer_entity
                            {
                                rcsb_id
                                entity_poly
                                { 
                                    pdbx_strand_id
                                }
                                rcsb_polymer_entity_annotation 
                                {
                                    annotation_id
                                    name
                                    description
                                    type
                                }
                                rcsb_entity_source_organism 
                                {
                                    ncbi_taxonomy_id
                                }
                            }
                      }
                }'''.replace('PDB_CHAIN',
                             ','.join([''.join(['"', pdbID.replace('_', '.'), '"']) for pdbID in structure_ids]))

        print(query)
        response = requests.post(request_url, json={'query': query})
        requested_annotations = defaultdict()
        if (response.status_code == 200):
            json_data = json.loads(response.text)
            # Parse the json response
            if (json_data['data'] is not None
                    and len(json_data['data']) > 0
                    and len(json_data['data']['polymer_entity_instances']) > 0):
                for polymer_instance in json_data['data']['polymer_entity_instances']:
                    requested_annotation = ','.join(['#'] * 3)
                    taxonomy_id = '#'
                    if (polymer_instance['polymer_entity'] is not None):
                        structure_id = polymer_instance['polymer_entity']['rcsb_id'].split('_')[0]
                        chain_ids = polymer_instance['polymer_entity']['entity_poly']['pdbx_strand_id'].split(',')
                        if (polymer_instance['polymer_entity']['rcsb_polymer_entity_annotation'] is not None):
                            annotations = self.split_go_annotations([annotation['annotation_id'] for annotation in
                                                                     polymer_instance['polymer_entity'][
                                                                       'rcsb_polymer_entity_annotation']])
                            requested_annotation = ','.join(annotations)
                        if (polymer_instance['polymer_entity']['rcsb_entity_source_organism'] is not None and len(
                                polymer_instance['polymer_entity']['rcsb_entity_source_organism']) > 0):
                            taxonomy_id = polymer_instance['polymer_entity']['rcsb_entity_source_organism'][0][
                                'ncbi_taxonomy_id']
                        requested_annotation = ','.join([repr(taxonomy_id), requested_annotation])
                        for chainID in chain_ids:
                            requested_annotations['_'.join([structure_id, chainID])] = requested_annotation
                for pdb_id in structure_ids:
                    if (pdb_id in requested_annotations):
                        results.append(requested_annotations[pdb_id])
                    else:
                        results.append(','.join(['#'] * 4))
        else:
            if (self.verbose is True):
                print(response.text)
        return results

    def get_go_term(self, go_id):
        # Retrieve information on a GO term ID from EMBL-EBI QuickGO
        time.sleep(random.randint(5, 15))
        print(' '.join(['Requesting GO term :', go_id]))
        response = requests.get(url=''.join([self._go_request_url, go_id]))
        go_term = 'n/a'
        term_type = '-'
        if (response.status_code == 200):
            data = response.json()
            go_term = data['results'][0]['name']
            term_type = data['results'][0]['aspect']
            self.go_cache[go_id] = (go_term, term_type)
        else:
            if (self.verbose is True):
                print(response.text)
        return go_term, term_type
