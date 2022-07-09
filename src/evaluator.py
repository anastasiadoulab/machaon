from collections import defaultdict
from src.structurealigner import StructureAligner
import os
from src.pdbhandler import PDBHandler
import traceback
import subprocess
from src.enricher import Enricher
import bs4 as bs
from src.viralrefseq import ViralRefSeq
from rdkit import DataStructs


class Evaluator:

    def __init__(self, viral=False):
        self.reference_data = defaultdict
        self.candidate_data = defaultdict
        self.reference_pdb_path = ''
        self.candidate_pdb_path = ''
        self._aligner = StructureAligner()
        self.root_disk = ''
        self.pdb_dataset_path = ''
        self.viral = viral
        self.verbose = True
        self._max_secondary_structure_length = 1500 # alignments involving greater lengths require great amount of time
        self.override_pdb_id = ''

    def set_root_disk(self, root_disk):
        self.root_disk = root_disk
        self._aligner.set_root_disk(root_disk)

    def set_pdb_path(self, pdb_path):
        self.pdb_dataset_path = pdb_path
        self._aligner.pdb_dataset_path = pdb_path

    def load_data(self, pdb_id, chain_id, refseq_id=''):
        full_pdb_file_path = os.path.sep.join([self.pdb_dataset_path, ''.join([pdb_id, '.pdb'])])
        pdbhandler = PDBHandler(pdb_id)
        pdbhandler.root_disk = self.root_disk
        full_pdb_file_path = pdbhandler.handle_alphafold_pdbid(full_pdb_file_path)
        # Get chain's ascending order number
        chain_index = pdbhandler.get_chain_index(full_pdb_file_path, chain_id)
        pdbhandler.verbose = self.verbose
        if(self.override_pdb_id == ''):
            pdbhandler.get_uniprot_accession_number(chain_id, pdb_id, full_pdb_file_path)
        else:
            # Assign an existing PDB ID, if this is a custom-named structure
            pdbhandler.get_uniprot_accession_number(chain_id, self.override_pdb_id)
        # Get existing residues in PDB data
        available_residues = pdbhandler.get_residue_range(full_pdb_file_path, chain_id)
        # Extract 1D (sequence) and 2D structures from PDB
        secondary_pdb_structure, primary_pdb_structure = self._aligner.get_structure([pdbhandler.structure_id, chain_id],
                                                                                     available_residues)
        # Retrieve protein sequence from UniProt
        protein_sequence = self.get_protein_sequence(pdbhandler.uniprot_accession_number)
        data = {'pdbId': pdb_id, 'chainId': chain_id, 'chainIndex': chain_index, '1D': protein_sequence,
                '1DPDB': primary_pdb_structure, '2D': secondary_pdb_structure, 'fullPDBPath': full_pdb_file_path}
        # Calculate molecular fingerprint
        try:
            chains = pdbhandler.fetch_pdb_peptidic_chains(full_pdb_file_path)
            data['fingerprint'] = pdbhandler.calculate_fingerprint(chains[chain_id]) if (chains is not False) else False
        except:
            data['fingerprint'] = False
            print(traceback.format_exc())
        # Retrieve transcript from RefSeq resources
        transcript_parts = []
        if (pdbhandler.uniprot_accession_number != ''):
            if (self.viral is True):
                info = self.get_viral_information(pdbhandler.uniprot_accession_number)
                if ('refSeqId' in info):
                    viral_data_handler = ViralRefSeq(info['refSeqId'], None if 'geneId' not in info else info['geneId'])
                    viral_data_handler.set_root_disk(self.root_disk)
                    viral_data_handler.fetch_viral_genome()
                    viral_data_handler.fetch_viral_genome_info()
                    if (len(viral_data_handler.genome_info) > 1):
                        refseq_data = viral_data_handler.parse_gene_data()
                        data = {**data, **refseq_data}
                        data['geneLength'] = len(viral_data_handler.sequence)
                        data['refSeqId'] = viral_data_handler.refseq_id
                    else:
                        self.viral = False
                        refseq_id = info['refSeqId']
            if (self.viral is False):
                if (refseq_id == ''):
                    refseq_id = self._aligner.read_ref_seq_id(pdbhandler.uniprot_accession_number)
                if (refseq_id != ''):
                    transcript_parts = self._aligner.fetch_transcript_parts(refseq_id)
                if (len(transcript_parts) > 0):
                    data['5UTR'] = transcript_parts[0]
                    data['CDS'] = transcript_parts[1]
                    data['3UTR'] = transcript_parts[2]
                    data['geneLength'] = transcript_parts[3]
                    data['refSeqId'] = refseq_id
        return data

    def get_protein_sequence(self, uniprot_accession):
        # Parsing UniProt XML
        sequence = ''
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        web_data = enricher.get_uniprot_data(uniprot_accession)
        if (len(web_data) > 2):
            xml_handler = bs.BeautifulSoup(web_data, "lxml")
            xml_info = xml_handler.find("entry")
            if (xml_info is not None):
                xml_info = xml_info.findAll('sequence')
                if (xml_info is not None):
                    for sequence_tag in xml_info:
                        if (sequence_tag.get('length') is not None and sequence_tag.get('checksum') is not None):
                            sequence = sequence_tag.text.strip()
        return sequence

    def get_viral_information(self, uniprot_accession):
        # Map UniProt accesion number with Entrez and RefSeq ids
        output = defaultdict()
        # Parsing XML
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        web_data = enricher.get_uniprot_data(uniprot_accession)
        if (len(web_data) < 2):
            return output
        xml_handler = bs.BeautifulSoup(web_data, 'lxml')

        xml_info = xml_handler.find('dbreference', {'type': 'GeneID'})
        if (xml_info is not None):
            output['geneId'] = xml_info.get('id')

        xml_info = xml_handler.findAll('dbreference', {'type': 'RefSeq'})
        if (xml_info is not None):
            for refseq_reference in xml_info:
                refseq_info = refseq_reference.find('property', {'type': 'nucleotide sequence ID'})
                if (refseq_info is not None):
                    output['refSeqId'] = refseq_info.get('value').split('.')[0]
        return output

    def load_go_data(self, data_row):
        result = defaultdict()
        term_types = ['molecularFunction', 'cellularComponent', 'biologicalProcess']
        for termType in term_types:
            result[termType] = [x for x in str(data_row[termType]).split('#') if x != '']
        return result

    def align_sequences(self, sections):
        # Align 1D / 2D sequences and calculate sequence identities
        result = defaultdict()
        for section in self.reference_data:
            if (section not in sections):
                continue
            identity = -1
            no_gap_identity = -1
            gaps = -1
            score = 0
            if (self.reference_data[section] != '' and
                section in self.candidate_data and
                self.candidate_data[section] != ''):
                config = defaultdict()
                alignment_type = 'OTHER'
                if (section in ['2D', '1DPDB']):
                    config['gapChars'] = ['{', ' ']
                    config['compressGaps'] = True
                # print(section)
                if (section == '1D'):
                    alignment_type = 'PROTEIN'
                # Avoid long sequences that lead to long running alignment processes
                skip_condition = section == '2D'
                skip_condition &= ((self._max_secondary_structure_length < len(self.reference_data[section])) or
                                   (self._max_secondary_structure_length < len(self.candidate_data[section])))
                if (skip_condition is False):
                    score, alignment, alignment_result = self._aligner.align(self.reference_data[section],
                                                                             self.candidate_data[section],
                                                                             alignment_type,
                                                                             config)
                    if (score is not False):
                        identity, no_gap_identity, gaps = self._aligner.calculate_identity(alignment)
                    else:
                        score = 0
            result[section] = (score, identity, no_gap_identity, gaps)
        return result

    def compare_gene_ontology(self):
        # Find common Gene Ontology entries between reference candidates
        result = []
        for term_type in self.reference_data['GO']:
            common = [x in self.reference_data['GO'][term_type] for x in self.candidate_data['GO'][term_type]].count(True)
            if (len(self.reference_data['GO'][term_type]) > 0):
                result.append(common / len(self.reference_data['GO'][term_type]))
            else:
                result.append(0.0)
        return result

    def calculate_tm_score(self):
        # Compute Template Modeling Score with TM-Align
        result = (-1, -1)
        try:
            output = subprocess.run(
                ['./TMalign', self.reference_data['fullPDBPath'], self.candidate_data['fullPDBPath'], '-split', '2', '-ter', '1', '-outfmt',
                 '2', '-chain1_id', str(self.reference_data['chainIndex']), '-chain2_id',
                 str(self.candidate_data['chainIndex'])], capture_output=True, cwd=os.getcwd())
            output = output.stdout.decode("utf-8").replace('PDBs/', '').replace('PDBs_vir/', '').replace('.pdb:', '_')
            output = output.split('\t')
            if (len(output) > 2):
                result = (float(output[2]), float(output[3]), float(output[-2]))
        except:
            if(self.verbose is True):
                print(traceback.format_exc())
        return result

    def compute_chemical_similarity(self):
        # Compute Tanimoto's Index
        result = -1
        if (self.reference_data['fingerprint'] is not False and self.candidate_data['fingerprint'] is not False):
            result = DataStructs.FingerprintSimilarity(self.reference_data['fingerprint'],
                                                       self.candidate_data['fingerprint'])
        return result
