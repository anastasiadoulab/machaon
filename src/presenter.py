import ast
import bs4 as bs
from src.enricher import Enricher
import pandas as pd
import numpy as np

from src.executionhandler import ExecutionHandler
from src.pdbhandler import PDBHandler
import os
import plotly.express as px
import seaborn as sns
from wordcloud import WordCloud
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import re
import pickle

plt.rcParams.update({'font.size': 8})


class Presenter:

    def __init__(self):
        self.root_disk = ''
        self.data_path = ''
        self.plots_output_path = ''
        self.go_output_path = ''
        self.verbose = False
        self._dataKeys = ['chainLength', 'pdbId', 'geneId', 'chainId', 'domain', 'alignedRange', 'referenceSegmentIndex',
                          'referenceRange', 'resolution', 'b-phipsi', 'w-rdist', 't-alpha', '3D-score', '1D-identity', '1D-identity-ng',
                          '1D-identity-gaps', '1D-content', '1DPDB-identity', '1DPDB-identity-ng', '1DPDB-identity-gaps',
                          '1DPDB-content', '2D-identity', '2D-identity-ng', '2D-identity-gaps', '2D-content', '5UTR-score',
                          '5UTR-identity', '5UTR-identity-ng', '5UTR-identity-gaps', '5UTR-content', 'CDS-score',
                          'CDS-identity', 'CDS-identity-ng', 'CDS-identity-gaps', 'CDS-content', '3UTR-score', '3UTR-identity',
                          '3UTR-identity-ng', '3UTR-identity-gaps',  '3UTR-content', 'chemSim', 'molecularFunctionSim',
                          'cellularComponentSim',  'biologicalProcessSim', 'geneLength', 'refSeqId']
        self._display_titles = {'3D-score': '3D Similarity', '1D-identity': '1D identity', '2D-identity': '2D identity',
                               '5UTR-identity': '5\'-UTR identity', '3UTR-identity': '3\'-UTR identity',
                               'CDS-identity': 'CDS identity', 'chemSim': 'Chemical similarity',
                                'molecularFunctionSim': 'Common functions',
                               'cellularComponentSim': 'Common locations', 'biologicalProcessSim': 'Common pathways'}
        self.taxonomy_depth = 100
        self.report_styling = \
            '<!DOCTYPE html><head><style> .column, .description, .item-separator {	float: left; width: 100%;} .edges {	width: 35%;}.middle {   width: 30%;} \
             .container{	padding: 10px;	overflow: visible;} .odd{border-color: bisque;} li{ word-break: break-word; } \
             hr { border-radius: 5px; border: 2px dashed silver; width: 85%;} .row_index { width: max-content; position: relative; \
             margin-top: -28px; font-size: 22pt; font-weight: bolder; color: silver;} body { font-family: \'Roboto\',sans-serif;} \
             h1{color: gray;} .goinfo {font-size: smaller;} .goheader{font-size: initial;}  @media print{ .description, .goinfo {display: none !important;} \
             .middle, .edges { width: 100%;}  ul, li {display: inline; padding: 2px;} .container { margin-bottom: 10px; } .column { padding: 5px; } } \
            </style><meta charset="UTF-8"></head><body><h1>'

        self.comparison_mode = ''
        self.fold_names = {'H': 'α-helix', 'B': 'β-bridge residue', 'b': 'β-bridge residue', 'E': 'extended strand', 'G': '3-helix', 'I': 'π-helix',
                           'P': 'PPII-helix', 'T': 'turn', 'S': 'bend', '.': 'loop', ' ': 'loop', 'C': 'loop', '?' : '?'}
        self.cores = 1
        self.debugging = False

        self.vector_format_output = False
        self.tiff_format_output = False

    def get_information(self, uniprot_accession):
        output = defaultdict()
        # Parsing UniProt metadata
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        web_data = enricher.get_uniprot_data(uniprot_accession)
        if (len(web_data) < 2):
            return output
        xml_handler = bs.BeautifulSoup(web_data, "lxml")
        output['uniprotId'] = uniprot_accession
        xml_info = xml_handler.find('protein')
        if (xml_info is not None):
            output['protein'] = xml_info.text.strip().split('\n')[0]
        xml_info = xml_handler.find('gene')
        if (xml_info is not None):
            output['gene'] = xml_info.text.strip().split('\n')[0]
        xml_info = xml_handler.find('organism')
        if (xml_info is not None):
            xml_info = xml_info.find('name', {'type': 'scientific'})
            if (xml_info is not None):
                output['organism'] = xml_info.text.strip().replace('\n', '|')
        xml_info = xml_handler.findAll('sequence')
        if (xml_info is not None):
            for sequence_tag in xml_info:
                if (sequence_tag.get('length') is not None and sequence_tag.get('checksum') is not None):
                    output['sequenceLength'] = sequence_tag.get('length')
        xml_info = xml_handler.find('comment', {'type': 'function'})
        if (xml_info is not None):
            output['function'] = xml_info.text.strip()
        xml_info = xml_handler.find('comment', {'type': 'subunit'})
        if (xml_info is not None):
            output['subunit'] = xml_info.text.strip()
        return output

    def create_node_label(self, name, population):
        # Create a label for a node in the graph visualization
        # displaying the species name and its total occurrences in the results
        viral_suffix_index = name.find('vir')
        if (viral_suffix_index != -1 and name[viral_suffix_index - 1] == ' '):
            viral_suffix_index = -1
        if (viral_suffix_index != -1):
            name = ''.join([name[:viral_suffix_index], '-\n', name[viral_suffix_index:]])
        name = name.replace(' ', '\n')
        return ''.join([name, '\n(', repr(population), ')']) if population > 1 else name

    def display_taxonomy(self, entry, file_name=''):
        pdb_id, chain_id, naming_extension = entry
        if (file_name == ''):
            if (len(naming_extension) > 0 and '_' not in naming_extension):
                naming_extension = ''.join(['_', naming_extension])
            file_name = os.path.join(self.data_path,
                                    ''.join([pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name,
                           '. You can also provide a filename explicitly as a second argument.']))
        evaluated_data = pd.read_csv(file_name, sep='\t')
        aggregated_taxonomy = defaultdict()
        top_roots = []
        max_depth = 0
        # Distinguish hierarchies of different taxonomies (e.g. Eukaryotes)
        for index, row in evaluated_data.iterrows():
            if (index != 0):
                pdbhandler = PDBHandler()
                pdbhandler.root_disk = self.root_disk
                pdbhandler.verbose = self.verbose
                pdbhandler.get_uniprot_accession_number(row['chainId'], row['pdbId'])
                if (pdbhandler.uniprot_accession_number != ''):
                    entry = pdbhandler.uniprot_accession_number, aggregated_taxonomy, top_roots, max_depth
                    aggregated_taxonomy, top_roots, max_depth = self.get_lineage_aggregated_information(entry)
        taxa = list(aggregated_taxonomy.keys())
        # Construct a graph for each hierarchy
        # Use a different color for each depth level
        color_palette = sns.color_palette('colorblind', self.taxonomy_depth)
        for top_node in top_roots:
            network_graph = nx.Graph()
            total_nodes = 0
            for taxon in taxa:
                node_label = self.create_node_label(taxon, aggregated_taxonomy[taxon]['population'])
                top_root = aggregated_taxonomy[taxon]['topRoot']
                if (top_node != top_root):
                    continue
                if (aggregated_taxonomy[taxon]['root'] != '-'):
                    root_taxon = aggregated_taxonomy[taxon]['root']
                    parent_node_label = self.create_node_label(root_taxon, aggregated_taxonomy[root_taxon]['population'])
                    network_graph.add_edge(node_label, parent_node_label, color=color_palette[aggregated_taxonomy[taxon]['level']],
                               weight=2 if aggregated_taxonomy[taxon]['level'] < 2 else 1)
                    total_nodes += 1
            fig, ax = plt.subplots()
            # fig = plt.figure(figsize=(15, 10) if total_nodes > 50 else (8, 8))
            fig.set_figwidth(15  if total_nodes > 50 else 8)
            fig.set_figheight(10 if total_nodes > 50 else 8)
            # Arrange nodes by total population
            if (total_nodes > 50 or 'Viruses' in top_roots):
                graph_layout = nx.kamada_kawai_layout(network_graph)
            else:
                graph_layout = nx.spring_layout(network_graph, k=0.12, iterations=80)
            # Each depth has a thinner branch than one of the previous one
            edges = network_graph.edges()
            colors = [network_graph[u][v]['color'] for u, v in edges]
            weights = [network_graph[u][v]['weight'] for u, v in edges]

            nx.draw_networkx(network_graph, pos=graph_layout, with_labels=True, node_color='black', font_size=8, node_size=3, edge_color=colors,
                    width=weights)
            # Store plot in different image formats
            output_filename = os.path.join(self.plots_output_path,
                                          ''.join([pdb_id, '_', chain_id, naming_extension, '-', top_node, '-tree']))
            plt.grid(False)
            # fig.set_facecolor('white')
            # ax = plt.axes()

            # # ax.axis('off')
            ax.set_facecolor('white')
            plt.savefig(''.join([output_filename, '.png']), format='png', dpi=300, bbox_inches='tight')
            if (self.vector_format_output is True):
                plt.savefig(''.join([output_filename, '.eps']), format='eps', dpi=300, bbox_inches='tight')
            if (self.tiff_format_output is True):
                plt.savefig(''.join([output_filename, '.tiff']), format='tiff', dpi=300, bbox_inches='tight')
            plt.clf()
            plt.figure(figsize=plt.rcParams.get('figure.figsize'))

    def get_lineage_aggregated_information(self, entry):
        # Parse lineage information from UniPro metadata in XML
        uniprot_accession, aggregated_taxonomy, top_roots, max_depth = entry
        lineage = defaultdict() if aggregated_taxonomy is None else aggregated_taxonomy
        roots = [] if top_roots is None else top_roots
        max_depth = 0 if max_depth is None else max_depth
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        web_data = enricher.get_uniprot_data(uniprot_accession)
        if (len(web_data) > 2):
            xml_handler = bs.BeautifulSoup(web_data, "lxml")
            xml_info = xml_handler.find('lineage')
            if (xml_info is not None):
                xml_info = xml_handler.findAll('taxon')
                if (xml_info is not None):
                    root = None
                    top_root = None
                    if (max_depth < len(xml_info)):
                        max_depth = len(xml_info)
                    # Parse taxonomy up to a depth level
                    for tagIndex, taxonTag in enumerate(xml_info):
                        if (self.taxonomy_depth < tagIndex + 1):
                            break
                        taxon_name = taxonTag.text.strip()
                        if (taxon_name not in lineage):
                            lineage[taxon_name] = defaultdict()
                            lineage[taxon_name]['population'] = 0
                        if ('level' not in lineage[taxon_name]):
                            lineage[taxon_name]['level'] = tagIndex
                        # Find the top root over every taxonomy list in the XML
                        # (multiple taxonomy hierarchies might be present)
                        if (root is not None):
                            if ('root' not in lineage[taxon_name]):
                                lineage[taxon_name]['root'] = root
                                lineage[taxon_name]['topRoot'] = top_root
                            if (root != lineage[taxon_name]['root'] and self.verbose is True):
                                print(''.join(['Different tree roots :', root, ' ', lineage[taxon_name]['root']]))
                        else:
                            lineage[taxon_name]['root'] = '-'
                            lineage[taxon_name]['topRoot'] = taxon_name
                            top_root = taxon_name
                            if (taxon_name not in roots):
                                roots.append(taxon_name)
                        root = taxon_name
                        lineage[taxon_name]['population'] += 1
        return lineage, roots, max_depth

    def create_word_cloud(self, entry, file_name=''):
        pdb_id, chain_id, naming_extension = entry
        execution_handler = ExecutionHandler(self.cores, 1 * 60)
        if (file_name == ''):
            if (len(naming_extension) > 0 and '_' not in naming_extension):
                naming_extension = ''.join(['_', naming_extension])
            file_name = os.path.join(self.data_path,
                                    ''.join([pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name,
                           '. You can also provide a filename explicitly as a second argument.']))
        evaluated_data = pd.read_csv(file_name, sep='\t')
        text = []
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdbhandler.verbose = self.verbose

        # Accumulate a word label for each protein in the results
        input_entries = [(row, pdbhandler) for index, row in evaluated_data.iterrows() if index != 0]
        results = []
        if self.debugging is False:
            results = execution_handler.parallelize(self.collect_words, input_entries)
        else:
            for entry in input_entries:
                results.append(self.collect_words(entry))
        if len(results) > 0:
            text = [item for sublist in results for item in sublist]
            # Calculate occurrences for each accumulated label
            frequencies = Counter(text)
            # Construct a word cloud
            wc = WordCloud(collocations=False).generate_from_frequencies(frequencies)
            # Create a visualization and store it in different formats
            plt.axis("off")
            plt.imshow(wc, interpolation="bilinear")
            output_filename = os.path.join(self.plots_output_path,
                                          ''.join([pdb_id, '_', chain_id, naming_extension, '-wordcloud']))
            plt.savefig(''.join([output_filename, '.png']), format='png', dpi=300, bbox_inches='tight')
            if (self.vector_format_output is True):
                plt.savefig(''.join([output_filename, '.eps']), format='eps', dpi=300, bbox_inches='tight')
            if (self.tiff_format_output is True):
                plt.savefig(''.join([output_filename, '.tiff']), format='tiff', dpi=300, bbox_inches='tight')
            plt.clf()
            plt.figure(figsize=plt.rcParams.get('figure.figsize'))

    def collect_words(self, entry):
        row, pdbhandler = entry
        words = []
        uniprot_accession_number = pdbhandler.get_uniprot_accession_number(row['chainId'], row['pdbId'])
        # Keep gene name for human proteins or organism name for proteins from different origin
        if (pdbhandler.uniprot_accession_number != ''):
            requested = self.get_information(uniprot_accession_number)
            if ('organism' in requested):
                organisms = requested['organism'].lower().split('|')
                if len(organisms) > 0:
                    if ('homo sapiens' in organisms):
                        if ('gene' in requested):
                            parts = requested['gene'].split('|')
                            if len(parts) > 0:
                                words = [parts[0]]
                    else:
                        words = organisms
        return words

    def display_structure_id(self, structure_id):
        # Check if the current entry is a predicted model
        matched = re.findall("MGYP[0-9]{12}", structure_id)
        alphafold = structure_id.startswith('AF-')
        structure_link = 'https://www.rcsb.org/structure/'
        identifier_label = 'PDB ID'
        if alphafold is True:
            structure_link = 'https://alphafold.ebi.ac.uk/entry/'
            identifier_label = 'AF ID'
        if len(matched) > 0:
            structure_link = 'https://esmatlas.com/resources/detail/'
            identifier_label = 'MGnify ID'
        return structure_link, identifier_label

    def create_report(self, entry, file_name='', reduced_data=False):
        execution_handler = ExecutionHandler(self.cores, 1 * 60)
        custom_file = file_name != ''
        pdb_id, chain_id, naming_extension = entry
        if (len(naming_extension) > 0 and '_' not in naming_extension):
            naming_extension = ''.join(['_', naming_extension])
        if (file_name == ''):
            file_name = os.path.join(self.data_path, ''.join([pdb_id, '_', chain_id, naming_extension,
                                     '-merged-enriched_eval.csv' if reduced_data is False else '-merged-notenriched.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name, '. You can also provide a filename explicitly as a second argument.']))
        candidates_data = pd.read_csv(file_name, sep='\t')
        # Construct an HTML presentation of the results
        output_html = [self.report_styling]
        title = ['Structural Comparison Report for ', pdb_id, '_', chain_id, naming_extension.replace('-metrics', '')]
        if 'domain' in candidates_data.columns:
            title.append(' - domains')
        elif 'referenceSegmentIndex' in candidates_data.columns:
            title.append('  - segments')
        else:
            title.append(' - whole structures')
        total_records = len(candidates_data.index)
        title.append(''.join([' (total: ', repr(total_records - 1 if custom_file is False else total_records), ')']))
        output_html.append(''.join(title+['</h1><div class="container"><div class="item-separator"><hr></div></div>']))
        # Full report
        if reduced_data is False:
            uniprot_info = defaultdict(list)
            if(custom_file is False):
                # Reference protein placeholders
                uniprot_info['ids'].append('#')
                uniprot_info['protein_names'].append('#')
                uniprot_info['gene_names'].append('#')
                uniprot_info['organism_names'].append('#')
            input_entries = [(row, index if custom_file is False else index + 1)
                             for index, row in candidates_data.iterrows() if index != 0 or custom_file is True]
            # Create a row with all computed and retrieved information for each protein in the results
            results = []
            if self.debugging is False:
                results = execution_handler.parallelize(self.create_report_row, input_entries)
            else:
                for entry in input_entries:
                    results.append(self.create_report_row(entry))
            for html_row, row_data in results:
                uniprot_info['ids'].append(row_data['uniprotId'] if 'uniprotId' in row_data else '#')
                uniprot_info['protein_names'].append(row_data['protein'] if 'protein' in row_data else '#')
                uniprot_info['gene_names'].append(row_data['gene'] if 'gene' in row_data else '#')
                uniprot_info['organism_names'].append(row_data['organism'] if 'organism' in row_data else '#')
                output_html.append(html_row)
            candidates_data['proteinName'] = uniprot_info['protein_names']
            candidates_data['geneName'] = uniprot_info['gene_names']
            candidates_data['organismName'] = uniprot_info['organism_names']
            candidates_data.to_csv(file_name.replace('.csv', '_full.csv'), sep='\t', index=False)
        else:
            # Create a row with all computations for each protein in the final cluster
            input_entries = [(index, row) for index, row in candidates_data.iterrows() if index != 0
                             or custom_file is True]
            results = []
            if self.debugging is False:
                results = execution_handler.parallelize(self.create_short_report_row, input_entries)
            else:
                for entry in input_entries:
                    results.append(self.create_short_report_row(entry))
            output_html.append('\n'.join(results))
        output_html.append('</body></html>')
        with open(file_name.replace('.csv', '_report.html'), 'w') as reportFile:
            reportFile.write(''.join(output_html))

    def create_short_report_row(self, entry):
        index, row = entry
        # Split structure id for display
        parts = []
        if 'structureId' in row:
            parts = row['structureId'].split('_')
        else:
            parts = list(row[['pdbId', 'chainId']])
        pdb_id = ''
        chain_id = ''
        if len(parts) > 1:
            pdb_id = parts[0]
            chain_id = parts[1]
        # Handle predicted structures
        alphafold = pdb_id.startswith('AF-')
        structure_link, identifier_label = self.display_structure_id(pdb_id)
        return ''.join(['<div class="container', ' odd' if index % 2 == 0 else '', '">',
                        '<div class="column edges"><div class="row_index">', repr(index),
                        '</div><ul>'
                        # PDB Specific information
                        '<li><b>', identifier_label, ':</b> <a href="', structure_link,
                         pdb_id if alphafold is False else pdb_id.split('-')[1], '">', pdb_id,
                        '</a> | <b>Chain:</b> ',
                        chain_id, '</li>',
                        self.create_list_float_item(row, 'b-phipsi', 'b-phipsi'),
                        self.create_list_float_item(row, 'w-rdist', 'w-rdist'),
                        self.create_list_float_item(row, 't-alpha', 't-alpha'),
                        '</ul></div></div><div class="container"><div class="item-separator">',
                        '<hr></div></div>',
                        '<div class="container"><div class="item-separator"><hr></div></div>'])

    def expose_info(self, data_row):
        # Create a dictionary containing all available information for a protein in the results
        info = defaultdict()
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        enricher.load_go_cache()
        term_types = ['molecularFunction', 'cellularComponent', 'biologicalProcess']
        for term_type in term_types:
            term_ids = [x for x in str(data_row[term_type]).split('#') if x != '']
            terms = []
            for term_id in term_ids:
                if (term_id not in enricher.go_cache):
                    term_info = enricher.get_go_term(term_id)
                    if (term_info[1] != '-'):
                        enricher.update_go_cache(term_id, term_info)
                else:
                    term_info = enricher.go_cache[term_id]
                terms.append(term_info[0])
            if (len(terms) > 0):
                info[term_type] = terms

        info = {**info, **self.retrieve_info(data_row, self._dataKeys)}
        return info

    def retrieve_info(self, data_row, data_keys):
        # Retrieve all accumulated data for a protein in the results
        result = defaultdict()
        for key in data_keys:
            if (key not in data_row):
                continue
            if (np.isreal(data_row[key]) and np.isnan(data_row[key]) == False):
                if (key not in ['chainLength', 'geneLength', 'resolution', 'b-phipsi', 'w-rdist', 't-alpha']
                        and 'identity-gaps' not in key and 'content' not in key):
                    if ((key in ['3D-score', 'resolution'] or 'identity' in key) and data_row[key] < 0):
                        continue
                    result[key] = repr(round(100 * data_row[key], 2))
                else:
                    if ('chainLength' in key or 'gaps' in key or 'geneLength' in key):
                        if(int(data_row[key]) > 0):
                            result[key] = repr(int(data_row[key]))
                    elif ('resolution' in key):
                        if(int(data_row[key]) > 0):
                            result[key] = repr(data_row[key])
                    elif (key in ['b-phipsi', 'w-rdist', 't-alpha']):
                        result[key] = repr(round(data_row[key], 6))
                    else:
                        result[key] = repr(data_row[key])
            elif (type(data_row[key]) == str):
                if(data_row[key] != '-'):
                    result[key] = data_row[key]
        return result

    # Format how the residue position ranges are displayed in the HTML report
    @staticmethod
    def sanitize_ranges(data):
        return data.replace('"', '').replace("'", '').replace(", ", '-').replace("]-[", ', ').replace('[', '').replace(
            ']', '')

    # Create an HTML output for an item in the metadata of a protein in the results
    @staticmethod
    def create_report_list_item(info, key, label, unit=''):
        return ''.join(['<li><b>', label, ':</b> ', ' '.join([info[key], unit]) if key in info else ' N/A', '</li>'])

    # Create an HTML output for a float item in the metadata of a protein in the results
    @staticmethod
    def create_list_float_item(info, key, label):
        return ''.join(['<li><b>', label, ':</b> ', repr(info[key]) if key in info else ' N/A', '</li>'])

    # Create an HTML link for an item in the metadata of a protein in the results
    @staticmethod
    def create_list_item_link(ref, info, key, label, unit=''):
        if key in info and info[key] != '#':
            return ''.join(['<li><b>', label, ':</b> <a href="', ref, info[key].split(',')[0].split(';')[0].strip(), '">', ' '.join([info[key], unit]), '</a></li>'])
        else:
            return ''.join(['<li><b>', label, ':</b> N/A</li>'])

    @staticmethod
    def create_content_list_item(info, key, label):
        output =  ' N/A'
        if key in info:
            if info[key] not in ['\'-\'', '-1.0']:
                output = info[key]
        return ''.join(['<li><b>', label, ':</b> ', output, '</li>'])

    # Create an HTML output for a protein in the results (a row in the HTML report)
    def create_report_row(self, entry):
        data_row, row_index = entry
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdbhandler.verbose = self.verbose
        pdbhandler.get_uniprot_accession_number(data_row['chainId'], data_row['pdbId'])
        requested = defaultdict()
        if (pdbhandler.uniprot_accession_number != ''):
            requested = self.get_information(pdbhandler.uniprot_accession_number)
        info = self.expose_info(data_row)
        info = {**requested, **info}

        # Genomic information label
        sequence_length_label = 'Sequence length'
        if ('refSeqId' in info):
            if (len(info['refSeqId']) > 1):
                sequence_length_label = ''.join(['Genomic' if info['refSeqId'][:2] == 'NC' else 'Transcript', ' sequence length'])

        # Check if the current entry is a predicted model
        alphafold = info['pdbId'].startswith('AF-')
        structure_link, identifier_label = self.display_structure_id(info['pdbId'])

        output_html = ' '.join([''.join(['<div class="container', ' odd' if row_index % 2 == 0 else '', '">',
        '<div class="column edges"><div class="row_index">', repr(row_index), '</div><ul>']),
        # Protein-level information
        self.create_report_list_item(info, 'protein', 'Protein name'),
        self.create_report_list_item(info, 'organism', 'Organism'),
        self.create_list_item_link('https://www.uniprot.org/uniprotkb/', info, 'uniprotId', 'Uniprot Accession Number'),
        self.create_report_list_item(info, 'sequenceLength', 'Protein sequence length', 'aa'),
        self.create_report_list_item(info, '1D-identity', '1D identity (%)'),
        self.create_report_list_item(info, '1D-identity-ng', '1D identity (%) [Gaps excluded]'),
        self.create_report_list_item(info, '1D-identity-gaps', '1D identity - Alignment Gaps'),
        self.create_content_list_item(info, '1D-content', '1D aligned content (&lt;aminoacid&gt;:%)'),
        self.create_report_list_item(info, 'molecularFunctionSim', 'Common reported functions (%)'),
        self.create_report_list_item(info, 'cellularComponentSim', 'Common reported locations (%)'),
        self.create_report_list_item(info, 'biologicalProcessSim', 'Common reported processes (%)'),
        '</ul></div><div class="column middle"><ul>',
        # PDB Specific information
        self.create_list_item_link(structure_link, info if alphafold is False else {'pdbId' : info['pdbId'].split('-')[1]}, 'pdbId', identifier_label),
        self.create_report_list_item(info, 'chainId', 'Chain'),
        self.create_report_list_item(info, 'chainLength', 'Crystallized protein length' if identifier_label == 'PDB ID' else 'Protein length', 'aa'),
        self.create_report_list_item(info, 'resolution', 'Resolution', 'Å'),
        ''.join(['<li><b>Associated domain:</b> ', info['domain'], '</li>']) if 'domain' in info else '',
        ''.join(['<li><b>Alinged residues range:</b> ', self.sanitize_ranges(repr(info['alignedRange'])), '</li>']) if 'alignedRange' in info else '',
        ''.join(['<li><b>Aligned to segment part (indices):</b> ', self.sanitize_ranges(repr(info['referenceSegmentIndex'])).replace('-', ', '), '</li>'])
        if 'referenceSegmentIndex' in info else '',
        ''.join(['<li><b>Alinged residues range of reference:</b> ', self.sanitize_ranges(repr(info['referenceRange'])), '</li>']) if 'referenceRange' in info else '',
        self.create_report_list_item(info, 'b-phipsi', 'b-phipsi' if self.comparison_mode != '' else ''.join(['b-phipsi (', self.comparison_mode, ')'])),
        self.create_report_list_item(info, 'w-rdist', 'w-rdist' if self.comparison_mode != '' else ''.join(['w-rdist (', self.comparison_mode, ')'])),
        self.create_report_list_item(info, 't-alpha', 't-alpha' if self.comparison_mode != '' else ''.join(['t-alpha (', self.comparison_mode, ')'])),
        self.create_report_list_item(info if 'chemSim' in info and float(info['chemSim']) > 0 else {}, 'chemSim', 'Chemical similarity (Tanimoto Index) (%)'),
        self.create_report_list_item(info, '1DPDB-identity', '1D identity (%) [PDB]'),
        self.create_report_list_item(info, '1DPDB-identity-ng', '1D identity (%) [Gaps excluded][PDB]'),
        self.create_report_list_item(info, '1DPDB-identity-gaps', '1D identity - Alignment Gaps [PDB]'),
        self.create_content_list_item(info, '1DPDB-content', '1D aligned content [PDB] (&lt;aminoacid&gt;:%)'),
        self.create_report_list_item(info, '2D-identity', '2D identity (%) [PDB]'),
        self.create_report_list_item(info, '2D-identity-ng', '2D identity (%) [Gaps excluded][PDB]'),
        self.create_report_list_item(info, '2D-identity-gaps', '2D identity - Alignment Gaps [PDB]'),
        self.create_content_list_item(info, '2D-content', '2D aligned content [PDB] (&lt;2D-fold&gt;:%)'),
        self.create_report_list_item(info, '3D-score', '3D similarity (TM-Score) (%) [PDB]'),
        '</ul></div><div class="column edges"><ul>',
        # Gene/Transcript specific information
        self.create_report_list_item(info, 'gene', 'Gene name'),
        self.create_list_item_link('https://www.ncbi.nlm.nih.gov/gene/', info, 'geneId', 'Entrez ID') if 'geneId' in info else '',
        self.create_list_item_link('https://www.ncbi.nlm.nih.gov/nuccore/', info, 'refSeqId', 'RefSeq ID'),
        self.create_report_list_item(info, 'geneLength', sequence_length_label),
        '<li><b>5-UTR|CDS|3-UTR identity (%):</b> ',
        ' | '.join([info[key] if key in info else 'N/A' for key in ['5UTR-identity', 'CDS-identity', '3UTR-identity']]),
        '</li>',
        '<li><b>5-UTR|CDS|3-UTR identity (%) [Gaps excluded]:</b> ',
        ' | '.join([info[key] if key in info else 'N/A' for key in ['5UTR-identity-ng', 'CDS-identity-ng', '3UTR-identity-ng']]),
        '</li>',
        '<li><b>5-UTR|CDS|3-UTR identity [Alignment Gaps]:</b> ',
        ' | '.join([info[key] if key in info and int(info[key]) > -1 else 'N/A' for key in ['5UTR-identity-gaps', 'CDS-identity-gaps', '3UTR-identity-gaps']]),
        '</li>',
        self.create_content_list_item(info, '5UTR-content', '5-UTR aligned content (&lt;base&gt;:%)'),
        self.create_content_list_item(info, 'CDS-content', 'CDS aligned content (&lt;base&gt;:%)'),
        self.create_content_list_item(info, '3UTR-content', '3-UTR aligned content (&lt;base&gt;:%)'),
        '</ul></div><div class="description"><br/><b>Uniprot Description:</b><br/><br/>',
        info['function'] if 'function' in info else 'N/A',
        ''.join(['<br/><br/>', info['subunit']]) if 'subunit' in info else 'N/A',
        '<br/><br/>',
        '<b class="goheader">Gene Ontology Information:</b><br/><br/>',
        '</div><div class="column edges goinfo">',
        '<u>Molecular Function</u>',
        ''.join(['<ul>'] + [''.join(['<li>', term, '</li>'])for term in info['molecularFunction']] + ['</ul>']) if 'molecularFunction' in info else '<br/><br/>N/A',
        '<br/><br/><br/><br/>',
        '</div>',
        '<div class="column middle goinfo">',
        '<u>Location</u>',
        ''.join(['<ul>'] + [''.join(['<li>', term, '</li>']) for term in info['cellularComponent']] + ['</ul>']) if 'cellularComponent' in info else '<br/><br/>N/A',
        '</div><div class="column edges goinfo">',
        '<u>Biological process</u>',
        ''.join(['<ul>'] + [''.join(['<li>', term, '</li>']) for term in info['biologicalProcess']] + ['</ul>']) if 'biologicalProcess' in info else '<br/><br/>N/A',
        '<br/></div></div><div class="container"><div class="item-separator"><hr></div></div>'])


        return ''.join(output_html), info

    # Create a polar plot that aggregates the values of a metric computed for each protein in the results
    def plot_results_polar(self, entry, file_name=''):
        pdb_id, chain_id, naming_extension, section = entry
        if (len(naming_extension) > 0 and '_' not in naming_extension):
            naming_extension = ''.join(['_', naming_extension])
        if (file_name == ''):
            file_name = os.path.join(self.data_path,
                                    ''.join([pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name,
                           '. You can also provide a filename explicitly as a second argument.']))
        evaluated_data = pd.read_csv(file_name, sep='\t')
        evaluated_data['mergedId'] = evaluated_data['pdbId'] + '_' + evaluated_data['chainId']
        evaluated_data = evaluated_data.loc[:, ['mergedId', section]]
        total_set_size = len(evaluated_data.index) - 1
        evaluated_data = evaluated_data.dropna()
        evaluated_data = evaluated_data[evaluated_data[section] >= 0]
        evaluated_data[section] *= 100.0
        heatmap_title = ''.join('Percentage (%)')
        evaluated_data[heatmap_title] = evaluated_data[section]
        subset_size = repr(len(evaluated_data[section].index))
        evaluated_data.sort_values(by=heatmap_title, inplace=True, ascending=False)
        plot_title = ''.join([self._display_titles[section], ' for ', subset_size, ' proteins in the final set (total: ',
                             repr(total_set_size), ')'])
        figure = px.bar_polar(evaluated_data, r=heatmap_title, theta='mergedId',
                              color=heatmap_title, template="plotly", barnorm='percent',
                              color_discrete_sequence=sns.color_palette('colorblind'), range_color=[0.0, 100.0]
                              ,range_r=[0.0, 100.0])
        figure.update_layout(
            font_size=12,
            polar=dict(
                radialaxis=dict(showticklabels=False, ticks=''),
                angularaxis=dict(showticklabels=False, ticks='')
            ),
            autosize=False,
            width=500,
            height=500,
            title={
                'text': plot_title,
                'xanchor': 'center',
                'yanchor': 'bottom',
                'y': 0.10,
                'x': 0.5},
            margin=dict(l=30, r=30, t=60, b=80)
        )

        output_filename = os.path.join(self.plots_output_path,
                                      ''.join([pdb_id, '_', chain_id, naming_extension, '_', section]))
        figure.write_image(''.join([output_filename, '.png']), format='png', scale=5)
        if (self.vector_format_output is True):
            try:
                figure.write_image(''.join([output_filename, '.eps']), format='eps', scale=5)
            except:
                pass
            figure.write_image(''.join([output_filename, '.pdf']), format='pdf', scale=5)
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))

    @staticmethod
    def parse_content(row_value):
        if ':' in repr(row_value):
            return ast.literal_eval(row_value)
        else:
            return False

    def plot_2D_coverage(self, entry, file_name=''):
        # Retrieve the results
        pdb_id, chain_id, naming_extension = entry
        if (file_name == ''):
            if (len(naming_extension) > 0 and '_' not in naming_extension):
                naming_extension = ''.join(['_', naming_extension])
            file_name = os.path.join(self.data_path,
                                     ''.join( [pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name,
                           '. You can also provide a filename explicitly as a second argument.']))
        evaluated_data = pd.read_csv(file_name, sep='\t')
        content_dict = defaultdict(list)

        # Parse 2D content percetnages from alignments
        data_dicts = [self.parse_content(row) for row in evaluated_data['2D-content'].tolist()]
        for data_dict in data_dicts:
            if data_dict is not False:
                for key in data_dict:
                    content_dict[key].append(data_dict[key])

        # Find the maximum occurrence for plotting (axis x limit)
        max_occurrence = -1
        for key in content_dict:
            content_dict[key].sort()
            if max_occurrence < len(content_dict[key]):
                max_occurrence = len(content_dict[key])

        # Print the percentages for each 2D fold type
        fold_types = content_dict.keys()
        if len(fold_types) > 0:
            color_palette = sns.color_palette('colorblind', len(fold_types))
            figure, axes = plt.subplots(len(fold_types), figsize=(15, 15))
            sub_handles = []
            sub_labels = []
            for type_index, fold_type in enumerate(fold_types):
                axes[type_index].plot([x for x in range(len(content_dict[fold_type]))], content_dict[fold_type],
                                      color=color_palette[type_index], label=self.fold_names[fold_type] if fold_type in self.fold_names else '?')
                axes[type_index].set_xlabel('Sample size')
                axes[type_index].set_ylim(0, 100.0)
                axes[type_index].set_xlim(0, max_occurrence)
                plt.subplots_adjust(hspace=0.5)
                handles, labels = axes[type_index].get_legend_handles_labels()
                sub_handles.extend(handles)
                sub_labels.extend(labels)

            figure.add_subplot(111, frameon=False)
            plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
            plt.ylabel('Aligned content(%)')

            figure.legend(sub_handles, sub_labels, loc='lower left', fancybox=True, framealpha=1.0,
                          bbox_to_anchor=(0.05, 0.11))

            # Export to multiple formats according to the configuration
            output_filename = os.path.join(self.plots_output_path,
                                           ''.join([pdb_id, '_', chain_id, naming_extension, '-2Dfold_coverage']))
            plt.tight_layout()
            plt.savefig(''.join([output_filename, '.png']), format='png', dpi=300, bbox_inches='tight')
            if (self.vector_format_output is True):
                plt.savefig(''.join([output_filename, '.eps']), format='eps', dpi=300, bbox_inches='tight')
            if (self.tiff_format_output is True):
                plt.savefig(''.join([output_filename, '.tiff']), format='tiff', dpi=300, bbox_inches='tight')
            plt.clf()
            plt.figure(figsize=plt.rcParams.get('figure.figsize'))

    def plot_2D_matches(self, entry):
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdbhandler.verbose = self.verbose
        pdb_id, chain_id, naming_extension, aggregated_mask, unavailable_residues, plot_output_path = entry
        if (len(naming_extension) > 0 and '_' not in naming_extension):
            naming_extension = ''.join(['_', naming_extension])

        if np.sum(aggregated_mask) > 0:
            # Plot unavailable residues
            plt.scatter(unavailable_residues, [0] * len(unavailable_residues), color='black', marker='|')

            # Plot 2D matches
            legend_elements = []
            plt.ylabel('2D Matches')
            plt.xlabel('Residue positions')
            if (len(unavailable_residues) > 0):
                legend_elements.append(plt.Rectangle((0, 0), 1, 1, label='Missing', fc='black'))
                plt.legend(handles=legend_elements, facecolor='white', framealpha=0.5, loc='best', handlelength=1,
                           handleheight=1.125)
            color_palette = sns.color_palette('colorblind', 1)
            plt.plot(range(1, len(aggregated_mask) + 1), aggregated_mask, color=color_palette[0])

            # Export to multiple formats according to the configuration
            output_filename = os.path.join(plot_output_path,
                                           ''.join([pdb_id, '_', chain_id, naming_extension, '-2Dmatches']))
            plt.tight_layout()
            plt.savefig(''.join([output_filename, '.png']), format='png', dpi=300, bbox_inches='tight')
            if (self.vector_format_output is True):
                plt.savefig(''.join([output_filename, '.eps']), format='eps', dpi=300, bbox_inches='tight')
            if (self.tiff_format_output is True):
                plt.savefig(''.join([output_filename, '.tiff']), format='tiff', dpi=300, bbox_inches='tight')
            plt.clf()
            plt.figure(figsize=plt.rcParams.get('figure.figsize'))
