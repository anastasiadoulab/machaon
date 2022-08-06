import bs4 as bs
from src.enricher import Enricher
import pandas as pd
import numpy as np
from src.pdbhandler import PDBHandler
import os
import plotly.express as px
import seaborn as sns
from wordcloud import WordCloud
from collections import Counter
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict

plt.rcParams.update({'font.size': 8})


class Presenter:

    def __init__(self):
        self.root_disk = ''
        self.data_path = ''
        self.plots_output_path = ''
        self.verbose = False
        self._dataKeys = ['chainLength', 'pdbId', 'chainId', 'domain', 'alignedRange', 'referenceSegmentIndex',
                          'referenceRange', 'resolution', 'b-phipsi', 'w-rdist', 't-alpha', '3D-score', '1D-identity', '1D-identity-ng',
                          '1D-identity-gaps', '1DPDB-identity', '1DPDB-identity-ng', '1DPDB-identity-gaps',
                          '2D-identity', '2D-identity-ng', '2D-identity-gaps', '5UTR-score', '5UTR-identity', '5UTR-identity-ng',
                          '5UTR-identity-gaps', 'CDS-score', 'CDS-identity', 'CDS-identity-ng', 'CDS-identity-gaps',
                          '3UTR-score', '3UTR-identity', '3UTR-identity-ng', '3UTR-identity-gaps',  'chemSim', 'molecularFunctionSim',
                          'cellularComponentSim',  'biologicalProcessSim', 'geneLength', 'refSeqId']
        self._display_titles = {'3D-score': '3D Similarity', '1D-identity': '1D identity', '2D-identity': '2D identity',
                               '5UTR-identity': '5\'-UTR identity', '3UTR-identity': '3\'-UTR identity',
                               'CDS-identity': 'CDS identity', 'chemSim': 'Chemical similarity',
                                'molecularFunctionSim': 'Common functions',
                               'cellularComponentSim': 'Common locations', 'biologicalProcessSim': 'Common pathways'}
        self.taxonomy_depth = 100
        self._uniprot_data_path = None

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
        self._uniprot_data_path = os.path.join(self.root_disk, 'uniprotdata')
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
            plt.figure(figsize=(15, 10) if total_nodes > 50 else (8, 8))
            # Arrange nodes by total population
            if (total_nodes > 50 or 'Viruses' in top_roots):
                graph_layout = nx.kamada_kawai_layout(network_graph)
            else:
                graph_layout = nx.spring_layout(network_graph, k=0.12, iterations=80)
            # Each depth has a thinner branch than one of the previous one
            edges = network_graph.edges()
            colors = [network_graph[u][v]['color'] for u, v in edges]
            weights = [network_graph[u][v]['weight'] for u, v in edges]
            nx.draw(network_graph, graph_layout, with_labels=True, node_color='black', font_size=8, node_size=3, edge_color=colors,
                    width=weights)
            # Store plot in different image formats
            output_filename = os.path.join(self.plots_output_path,
                                          ''.join([pdb_id, '_', chain_id, naming_extension, '-', top_node, '-tree']))
            plt.savefig(''.join([output_filename, '.png']), format='png', dpi=300, bbox_inches='tight')
            plt.savefig(''.join([output_filename, '.eps']), format='eps', dpi=300, bbox_inches='tight')
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
        if (file_name == ''):
            if (len(naming_extension) > 0 and '_' not in naming_extension):
                naming_extension = ''.join(['_', naming_extension])
            file_name = os.path.join(self.data_path,
                                    ''.join([pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name,
                           '. You can also provide a filename explicitly as a second argument.']))
        self._uniprot_data_path = os.path.join(self.root_disk, 'uniprotdata')
        evaluated_data = pd.read_csv(file_name, sep='\t')
        text = []
        # Accumulate a word label for each protein in the results
        for index, row in evaluated_data.iterrows():
            if (index != 0):
                pdbhandler = PDBHandler()
                pdbhandler.root_disk = self.root_disk
                pdbhandler.verbose = self.verbose
                pdbhandler.get_uniprot_accession_number(row['chainId'], row['pdbId'])
                # Keep gene name for human proteins or organism name for proteins from different origin
                if (pdbhandler.uniprot_accession_number != ''):
                    requested = self.get_information(pdbhandler.uniprot_accession_number)
                    if('organism' in requested):
                        organisms = requested['organism'].lower().split('|')
                        if ('homo sapiens' in organisms):
                            if('gene' in requested):
                                text.append(requested['gene'].split('|')[0])
                        else:
                            text.extend(organisms)
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
        plt.savefig(''.join([output_filename, '.eps']), format='eps', dpi=300, bbox_inches='tight')
        plt.clf()
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))

    def create_report(self, entry, file_name=''):
        custom_file = file_name != ''
        pdb_id, chain_id, naming_extension = entry
        if (len(naming_extension) > 0 and '_' not in naming_extension):
            naming_extension = ''.join(['_', naming_extension])
        if (file_name == ''):
            file_name = os.path.join(self.data_path, ''.join([pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval.csv']))
        if (os.path.exists(file_name) is False):
            print(''.join(['The file does not exist: ', file_name, '. You can also provide a filename explicitly as a second argument.']))
        self._uniprot_data_path = os.path.join(self.root_disk, 'uniprotdata')
        evaluated_data = pd.read_csv(file_name, sep='\t')
        # Construct an HTML presentation of the results
        output_html = [
            '<!DOCTYPE html><head><style> .column, .description, .item-separator {	float: left; width: 100%;} .edges {	width: 35%;}.middle {   width: 30%;}',
            ' .container{	padding: 10px;	overflow: visible;} .odd{border-color: bisque;} li{ word-break: break-word; }',
            'hr { border-radius: 5px; border: 2px dashed silver; width: 85%;} .row_index { width: max-content; position: relative;'
            'margin-top: -28px; font-size: 22pt; font-weight: bolder; color: silver;} body { font-family: \'Roboto\',sans-serif;}',
            ' h1{color: gray;} .goinfo {font-size: smaller;} .goheader{font-size: initial;}  @media print{ .description, .goinfo {display: none !important;}',
            ' .middle, .edges { width: 100%;}  ul, li {display: inline; padding: 2px;} .container { margin-bottom: 10px; } .column { padding: 5px; } }',
            '</style><meta charset="UTF-8"></head><body><h1>']
        title = ['Structural Comparison Report for ', pdb_id, '_', chain_id, naming_extension.replace('-metrics', '')]
        if 'domain' in evaluated_data.columns:
            title.append(' - domains')
        elif 'referenceSegmentIndex' in evaluated_data.columns:
            title.append('  - segments')
        else:
            title.append(' - whole structures')
        total_records = len(evaluated_data.index)
        title.append(''.join([' (total: ', repr(total_records - 1 if custom_file is False else total_records), ')']))
        output_html.append(''.join(title))
        output_html.append('</h1><div class="container"><div class="item-separator"><hr></div></div>')
        uniprot_info = defaultdict(list)
        if(custom_file is False):
            # Reference protein placeholders
            uniprot_info['ids'].append('#')
            uniprot_info['protein_names'].append('#')
            uniprot_info['gene_names'].append('#')
            uniprot_info['organism_names'].append('#')
        # Create a row with all computed and retrieved information for each protein in the results
        for index, row in evaluated_data.iterrows():
            # Skip reference protein
            if (index == 0 and custom_file is False):
                continue
            html_row, row_data = self.create_report_row(row, index if custom_file is False else index + 1)
            uniprot_info['ids'].append(row_data['uniprotId'] if 'uniprotId' in row_data else '#')
            uniprot_info['protein_names'].append(row_data['protein'] if 'protein' in row_data else '#')
            uniprot_info['gene_names'].append(row_data['gene'] if 'gene' in row_data else '#')
            uniprot_info['organism_names'].append(row_data['organism'] if 'organism' in row_data else '#')
            output_html.append(html_row)
            output_html.append('<div class="container"><div class="item-separator"><hr></div></div>')
        output_html.append('</body></html>')
        evaluated_data['uniprotId'] = uniprot_info['ids']
        evaluated_data['proteinName'] = uniprot_info['protein_names']
        evaluated_data['geneName'] = uniprot_info['gene_names']
        evaluated_data['organismName'] = uniprot_info['organism_names']
        evaluated_data.to_csv(file_name.replace('.csv', '_full.csv'), sep='\t', index=False)
        with open(file_name.replace('.csv', '_report.html'), 'w') as reportFile:
            reportFile.write(''.join(output_html))

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
                if (key not in ['chainLength', 'resolution', 'b-phipsi', 'w-rdist', 't-alpha'] and 'identity-gaps' not in key):
                    if ((key in ['3D-score', 'resolution'] or 'identity' in key) and data_row[key] < 0):
                        continue
                    result[key] = repr(round(100 * data_row[key], 2))
                else:
                    if ('chainLength' in key or 'gaps' in key):
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
    def sanitize_ranges(self, data):
        return data.replace('"', '').replace("'", '').replace(", ", '-').replace("]-[", ', ').replace('[', '').replace(
            ']', '')

    # Create an HTML output for an item in the metadata of a protein in the results
    def create_report_list_item(self, info, key, label, unit=''):
        return ''.join(['<li><b>', label, ':</b> ', ' '.join([info[key], unit]) if key in info else 'N/A', '</li>'])

    # Create an HTML output for a protein in the results (a row in the HTML report)
    def create_report_row(self, data_row, row_index):
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdbhandler.verbose = self.verbose
        pdbhandler.get_uniprot_accession_number(data_row['chainId'], data_row['pdbId'])
        requested = defaultdict()
        if (pdbhandler.uniprot_accession_number != ''):
            requested = self.get_information(pdbhandler.uniprot_accession_number)
        info = self.expose_info(data_row)
        info = {**requested, **info}

        output_html = [''.join(['<div class="container', ' odd' if row_index % 2 == 0 else '', '">'])]
        output_html.append(''.join(['<div class="column edges"><div class="row_index">', repr(row_index), '</div><ul>']))
        # Protein-level information
        output_html.append(self.create_report_list_item(info, 'protein', 'Protein name'))
        output_html.append(self.create_report_list_item(info, 'organism', 'Organism'))
        output_html.append(self.create_report_list_item(info, 'uniprotId', 'Uniprot Accession Number'))
        output_html.append(self.create_report_list_item(info, 'sequenceLength', 'Protein sequence length', 'aa'))
        output_html.append(self.create_report_list_item(info, '1D-identity', '1D identity (%)'))
        output_html.append(self.create_report_list_item(info, '1D-identity-ng', '1D identity (%) [Gaps excluded]'))
        output_html.append(self.create_report_list_item(info, '1D-identity-gaps', '1D identity - Alignment Gaps'))
        output_html.append(self.create_report_list_item(info, 'molecularFunctionSim', 'Common reported functions (%)'))
        output_html.append(self.create_report_list_item(info, 'cellularComponentSim', 'Common reported locations (%)'))
        output_html.append(self.create_report_list_item(info, 'biologicalProcessSim', 'Common reported processes (%)'))
        output_html.append('</ul></div><div class="column middle"><ul>')
        # PDB Specific information
        output_html.append(self.create_report_list_item(info, 'pdbId', 'PDB ID'))
        output_html.append(self.create_report_list_item(info, 'chainId', 'Chain'))
        output_html.append(self.create_report_list_item(info, 'chainLength', 'Crystallized protein length', 'aa'))
        output_html.append(self.create_report_list_item(info, 'resolution', 'Resolution', 'Å'))
        if ('domain' in info):
            output_html.append(''.join(['<li><b>Associated domain:</b> ', info['domain'], '</li>']))
        if ('alignedRange' in info):
            output_html.append(''.join(['<li><b>Alinged residues range:</b> ', self.sanitize_ranges(repr(info['alignedRange'])), '</li>']))
        if ('referenceSegmentIndex' in info):
            output_html.append(''.join(['<li><b>Aligned to segment part (indices):</b> ', self.sanitize_ranges(repr(info['referenceSegmentIndex'])).replace('-', ', '), '</li>']))
        if ('alignedRange' in info):
            output_html.append(''.join(['<li><b>Alinged residues range of reference:</b> ', self.sanitize_ranges(repr(info['referenceRange'])), '</li>']))
        output_html.append(self.create_report_list_item(info, 'b-phipsi', 'b-phipsi'))
        output_html.append(self.create_report_list_item(info, 'w-rdist', 'w-rdist'))
        output_html.append(self.create_report_list_item(info, 't-alpha', 't-alpha'))
        output_html.append(self.create_report_list_item(info if 'chemSim' in info and float(info['chemSim']) > 0 else {},
                                                        'chemSim', 'Chemical similarity (Tanimoto Index) (%)'))
        output_html.append(self.create_report_list_item(info, '1DPDB-identity', '1D identity (%) [PDB]'))
        output_html.append(self.create_report_list_item(info, '1DPDB-identity-ng', '1D identity (%) [Gaps excluded][PDB]'))
        output_html.append(self.create_report_list_item(info, '1DPDB-identity-gaps', '1D identity - Alignment Gaps [PDB]'))
        output_html.append(self.create_report_list_item(info, '2D-identity', '2D identity (%) [PDB]'))
        output_html.append(self.create_report_list_item(info, '2D-identity-ng', '2D identity (%) [Gaps excluded][PDB]'))
        output_html.append(self.create_report_list_item(info, '2D-identity-gaps', '2D identity - Alignment Gaps [PDB]'))
        output_html.append(self.create_report_list_item(info, '3D-score', '3D similarity (TM-Score) (%) [PDB]'))
        output_html.append('</ul></div><div class="column edges"><ul>')
        # Gene/Transcript specific information
        output_html.append(self.create_report_list_item(info, 'gene', 'Gene name'))
        output_html.append(self.create_report_list_item(info, 'refSeqId', 'RefSeq ID'))
        sequence_length_label = 'Sequence length'
        if('refSeqId' in info):
            if(len(info['refSeqId'])> 1):
                sequence_length_label = ''.join([ 'Genomic' if info['refSeqId'][:2] == 'NC' else 'Transcript', ' sequence length'])
        output_html.append(self.create_report_list_item(info, 'geneLength', sequence_length_label))
        output_html.append('<li><b>5-UTR|CDS|3-UTR identity (%):</b> ')
        for keyIndex, key in enumerate(['5UTR-identity', 'CDS-identity', '3UTR-identity']):
            if (key in info):
                output_html.append(info[key])
            else:
                output_html.append('N/A')
            if (keyIndex < 2):
                output_html.append(' | ')
        output_html.append('</li>')
        output_html.append('<li><b>5-UTR|CDS|3-UTR identity (%) [Gaps excluded]:</b> ')
        for keyIndex, key in enumerate(['5UTR-identity-ng', 'CDS-identity-ng', '3UTR-identity-ng']):
            if (key in info):
                output_html.append(info[key])
            else:
                output_html.append('N/A')
            if (keyIndex < 2):
                output_html.append(' | ')
        output_html.append('</li>')
        output_html.append('<li><b>5-UTR|CDS|3-UTR identity [Alignment Gaps]:</b> ')
        for keyIndex, key in enumerate(['5UTR-identity-gaps', 'CDS-identity-gaps', '3UTR-identity-gaps']):
            gaps_output = 'N/A'
            if (key in info):
                if(int(info[key]) > -1):
                    gaps_output = info[key]
            output_html.append(gaps_output)
            if (keyIndex < 2):
                output_html.append(' | ')
        output_html.append('</li>')
        output_html.append('</ul></div><div class="description"><br/><b>Uniprot Description:</b><br/><br/>')
        if ('function' in info):
            output_html.append(info['function'])
        else:
            output_html.append('N/A')
        if ('subunit' in info):
            output_html.append('<br/><br/>')
            output_html.append(info['subunit'])
        output_html.append('<br/><br/>')

        output_html.append('<b class="goheader">Gene Ontology Information:</b><br/><br/>')
        output_html.append('</div><div class="column edges goinfo">')
        output_html.append('<u>Molecular Function</u>')
        if ('molecularFunction' in info):
            output_html.append('<ul>')
            for term in info['molecularFunction']:
                output_html.append('<li>')
                output_html.append(term)
                output_html.append('</li>')
            output_html.append('</ul>')
        else:
            output_html.append('<br/><br/>N/A')
        output_html.append('<br/><br/><br/><br/>')
        output_html.append('</div>')
        output_html.append('<div class="column middle goinfo">')
        output_html.append('<u>Location</u>')
        if ('cellularComponent' in info):
            output_html.append('<ul>')
            for term in info['cellularComponent']:
                output_html.append('<li>')
                output_html.append(term)
                output_html.append('</li>')
            output_html.append('</ul>')
        else:
            output_html.append('<br/><br/>N/A')
        output_html.append('</div><div class="column edges goinfo">')
        output_html.append('<u>Biological process</u>')
        if ('biologicalProcess' in info):
            output_html.append('<ul>')
            for term in info['biologicalProcess']:
                output_html.append('<li>')
                output_html.append(term)
                output_html.append('</li>')
            output_html.append('</ul>')
        else:
            output_html.append('<br/><br/>N/A')
        output_html.append('<br/></div></div>')

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
        figure.write_image(''.join([output_filename, '.eps']), format='eps', scale=5)
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))
