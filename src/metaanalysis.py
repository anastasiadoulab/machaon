import os
import traceback
import pandas as pd
from src.evaluator import Evaluator
import io
from src.structurealigner import StructureAligner
from src.presenter import Presenter
import pickle
import matplotlib.pyplot as plt
import numpy as np
from src.enricher import Enricher
import seaborn as sns
from src.pdbhandler import PDBHandler
import random
from collections import defaultdict
from src.executionhandler import ExecutionHandler
import uuid
from matplotlib.lines import Line2D
from collections import Counter
import subprocess

plt.rcParams.update({'font.size': 8})
random.seed(42)



class MetaAnalysis:

    def __init__(self):
        self.root_disk = ''
        self.debugging = False
        self.reference_pdbid = []
        self.reference_chain_id = []
        self.reference_segments = []
        self.is_reference_viral = True
        self.output_path = ''
        self.pdb_dataset_path = ''
        self.go_output_path = ''
        self.cores = 1
        self.override_pdb_id = ''
        self.alignment_levels = ['primary', 'secondary', 'mixed', 'hydrophobicity']
        self.go_alignment_level = 'secondary'
        self.segments = []
        self.known_areas = []
        self.viral_content = True
        self.vector_format_output = False
        self.tiff_format_output = False
        self.alignment_backend = ''
        self.verbose = False

    def evaluate(self, entry):
        # Extended comparisons of a finalist protein with the reference.
        # - 1D: protein sequence
        # - 1DPDB: protein sequence from PDB
        # - 2D: secondary structure from PDB
        # - 5UTR/CDS/3UTR: parts of transcript sequence
        # - 3D: TM-Score
        # - Tanimoto Index
        # - Common GO terms of 3 categories: 'molecularFunction', 'cellularComponent', 'biologicalProcess'
        reference_data, data_row, is_viral, tmalign_available = entry
        id_parts = data_row['pdbId'], data_row['chainId']
        evaluator = Evaluator(self.alignment_backend, is_viral)
        evaluator.set_root_disk(self.root_disk)
        evaluator.set_pdb_path(os.path.sep.join([self.root_disk, self.pdb_dataset_path]))
        evaluator.reference_data = reference_data
        evaluator.candidate_data = evaluator.load_data(*id_parts)
        evaluator.candidate_data['GO'] = evaluator.load_go_data(data_row)
        result = False
        try:
            sections = ['1D', '1DPDB', '2D', '5UTR', 'CDS', '3UTR']
            compared = evaluator.align_sequences(sections)
            sections.append('3D')
            if(tmalign_available):
                compared['3D'] = evaluator.calculate_tm_score()
            output = [id_parts[0], id_parts[1]]
            for section in sections:
                if (section in compared):
                    output.extend([repr(x) for x in compared[section]])
                else:
                    # every section has 5 values: score, identity,
                    # no gaps identity, identity gaps & content
                    output.extend(['-1'] * 5)
            output.append(repr(evaluator.compute_chemical_similarity()))
            output.extend([repr(x) for x in evaluator.compare_gene_ontology()])
            if('geneLength' in evaluator.candidate_data):
                output.append(repr(evaluator.candidate_data['geneLength']))
                output.append(evaluator.candidate_data['refSeqId'])
            else:
                output.extend(['-', '-'])
            result = '\t'.join(output)
        except:
            print(''.join(['Failed: ', repr(id_parts)]))
            print(traceback.format_exc())
        return result

    def evaluate_results(self):
        # Extended comparisons for each entry in the results run in parallel
        reference_evaluator = Evaluator(self.alignment_backend)
        reference_evaluator.set_root_disk(self.root_disk)
        reference_evaluator.set_pdb_path(os.path.sep.join([self.root_disk, self.pdb_dataset_path]))
        reference_evaluator.override_pdb_id = self.override_pdb_id
        for reference_segment in self.reference_segments:
            reference_evaluator.viral = self.is_reference_viral
            reference_protein = (self.reference_pdbid, self.reference_chain_id)
            input_entries = []
            naming_extension = reference_segment
            if ('_site' in naming_extension):
                naming_extension = ''.join([reference_segment, '-metrics'])
            elif (len(naming_extension) > 0):
                naming_extension = ''.join(['_', naming_extension])
            output_filepath = ''.join([self.output_path, os.path.sep, reference_protein[0], '_', reference_protein[1], naming_extension, '-merged-enriched_eval.csv'])
            if(os.path.exists(output_filepath)):
                print('Results are already evaluated for this entry: ', ''.join([reference_protein[0], '_', reference_protein[1], naming_extension]),
                      '. If you did not expect this outcome please delete this file and execute again:\n', output_filepath)
                continue
            reference_data = reference_evaluator.load_data(*reference_protein)

            # Log the reference data
            data_log = { key: [reference_data[key]] for key in ['pdbId', 'chainId', '1D', '2D', '1DPDB', 'CDS', '5UTR', '3UTR', 'geneLength', 'refSeqId'] if key in reference_data}
            if '2D' in data_log:
                if len(data_log['2D'][0]) > 0:
                    struct_content = Counter(data_log['2D'][0])
                    total = sum(struct_content.values(), 0.0)
                    struct_content = {k: (v/total) * 100.0 for k, v in struct_content.items() if k != '-'}
                    data_log['2DContent'] = [repr(struct_content)]
            data_file_name = ''.join(
                [self.output_path, os.path.sep, reference_protein[0], '_', str(reference_protein[1]), naming_extension,
                 '-data.csv'])
            pd.DataFrame(data_log).to_csv(data_file_name, index=False, sep='\t')

            # Prepare the input for evaluation
            file_name = ''.join([self.output_path, os.path.sep, reference_protein[0], '_', str(reference_protein[1]), naming_extension, '-merged-enriched.csv'])
            data = pd.read_csv(file_name, sep='\t')
            results = []
            tmalign_available = os.path.exists("./TMalign")
            for rowIndex, row in data.iterrows():
                if (rowIndex == 0):
                    reference_data['GO'] = reference_evaluator.load_go_data(row)
                else:
                    input_entries.append([reference_data, row, self.viral_content, tmalign_available])

            # Evaluation
            if (self.debugging is False):
                execution_handler = ExecutionHandler(self.cores, 12 * 60)
                results = execution_handler.parallelize(self.evaluate, input_entries)
            else:
                fails = 0
                for input_entry in input_entries:
                    result = self.evaluate(input_entry)
                    if (result is not False):
                        results.append(result)
                    else:
                        fails += 1
                print(''.join(['fails : ', repr(fails)]))

            # Store computed data to csv
            evaluation_data = pd.read_csv((io.StringIO('\n'.join(results))),
                                          names=['pdbId', 'chainId', '1D-score', '1D-identity', '1D-identity-ng',
                                                 '1D-identity-gaps', '1D-content', '1DPDB-score', '1DPDB-identity',
                                                 '1DPDB-identity-ng', '1DPDB-identity-gaps', '1DPDB-content',
                                                 '2D-score', '2D-identity', '2D-identity-ng', '2D-identity-gaps',
                                                 '2D-content', '5UTR-score', '5UTR-identity', '5UTR-identity-ng',
                                                 '5UTR-identity-gaps', '5UTR-content', 'CDS-score', 'CDS-identity',
                                                 'CDS-identity-ng', 'CDS-identity-gaps', 'CDS-content', '3UTR-score',
                                                 '3UTR-identity', '3UTR-identity-ng', '3UTR-identity-gaps', '3UTR-content',
                                                 '3D-score', '3D-score-cand', 'length', 'chemSim', 'molecularFunctionSim',
                                                 'cellularComponentSim', 'biologicalProcessSim', 'geneLength', 'refSeqId'],
                                          sep='\t')
            evaluated_data = pd.merge(data, evaluation_data, how='left', on=['pdbId', 'chainId'])
            evaluated_data.to_csv(output_filepath, index=False, sep='\t')

    @staticmethod
    def store_property_alignment(entry):
        # Append the alignment result to a dictionary
        aggregation_dict, pdb_id, chain_id, alignment_info, sequence, reference_residue_range, total_residue_range, aligner = entry
        aggregation_dict['mask'].append(alignment_info['mask'])
        alignment_details = ''
        identity = -1
        if(len(alignment_info['mask']) > 0):
            identity, no_gap_identity, _ = aligner.calculate_identity(alignment_info['result'])
            if(identity > 0):
                alignment_details = ''.join(['PDB ID: ', pdb_id, '.', chain_id, ' | Alignment (sequence identity: ',
                                            repr(round(identity * 100, 4)), '%):\n',
                                            '\n'.join(alignment_info['result'])])
        aggregation_dict['identity'].append(identity)
        aggregation_dict['alignment'].append(alignment_details)
        aggregation_dict['structure_id'].append(''.join([pdb_id, '_', chain_id]))

    def align_candidate_by_property(self, entry):
        # Perform alignments on user-specified levels (e.g. secondary structure) between proteins having
        # a GO property of interest and the reference
        row, reference_sequences, reference_residues = entry
        alignments = defaultdict()
        reference_residue_range = [x for x in range(reference_residues['fullRange'][0], reference_residues['fullRange'][-1] + 1)]
        aligner = StructureAligner(self.alignment_backend)
        aligner.set_root_disk(self.root_disk)
        pdb_data_root = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        aligner.pdb_dataset_path = pdb_data_root
        pdbhandler = PDBHandler(row['pdbId'])
        full_pdb_file_path = os.path.join(pdb_data_root, '.'.join([row['pdbId'], 'pdb']))
        residue_range = pdbhandler.get_residue_range(full_pdb_file_path, row['chainId'])
        available_residues = residue_range['fullRange']
        total_residue_range = [x for x in range(available_residues[0], available_residues[-1] + 1)]
        sequences = aligner.get_all_sequences([pdbhandler.structure_id, row['chainId']], residue_range)
        config = {'gapChars': ['{', ' ']}
        property_alignments = defaultdict()
        if(len(sequences) > 0):
            # Perform alignments for each specified level
            for alignment_level in self.alignment_levels:
                if(alignment_level not in sequences):
                    continue
                alignments[alignment_level] = defaultdict()
                # Choose scoring strategies for alignments
                scoring_strategy = 'OTHER'
                if(alignment_level == 'primary'):
                    scoring_strategy = 'PROTEIN'
                result = aligner.align(reference_sequences[alignment_level], sequences[alignment_level], scoring_strategy, config)
                alignments[alignment_level]['score'], alignments[alignment_level]['alignment'], alignments[alignment_level]['result'] = result
                alignments[alignment_level]['mask'] = aligner.get_alignment_mask(alignments[alignment_level]['alignment']) \
                                                      if alignments[alignment_level]['alignment'] is not False else []
            # Aggregate alignments by GO type
            for property_type in ['molecularFunction', 'cellularComponent', 'biologicalProcess']:
                go_ids = row[property_type].split('#')
                for goId in go_ids:
                    if (len(goId.strip('')) == 0):
                        continue
                    if (property_type not in property_alignments):
                        property_alignments[property_type] = defaultdict()
                    if (goId not in property_alignments[property_type]):
                        property_alignments[property_type][goId] = defaultdict()
                        for alignment_level in self.alignment_levels:
                            if (alignment_level not in alignments):
                                continue
                            if (alignment_level not in property_alignments[property_type][goId]):
                                property_alignments[property_type][goId][alignment_level] = defaultdict(list)
                            entry = property_alignments[property_type][goId][alignment_level], row['pdbId'], row['chainId'], alignments[alignment_level], \
                                    sequences[alignment_level], reference_residue_range, total_residue_range, aligner
                            self.store_property_alignment(entry)
        return property_alignments, ['_'.join([row['pdbId'], row['chainId']]), alignments]

    def aggregate_go_information(self):
        # Perform alignments in parallel for each protein in the results
        # and construct a dictionary with this information, indexed by GO properties
        if (os.path.exists(self.go_output_path) is False):
            os.makedirs(self.go_output_path)
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        aligner = StructureAligner(self.alignment_backend)
        aligner.set_root_disk(self.root_disk)
        aligner.pdb_dataset_path = pdb_dataset_path
        for referenceSegment in self.reference_segments:
            reference_protein = (self.reference_pdbid, self.reference_chain_id)
            naming_extension = referenceSegment
            if ('_site' in naming_extension):
                naming_extension = ''.join([referenceSegment, '-metrics'])
            elif (len(naming_extension) > 0):
                naming_extension = ''.join(['_', naming_extension])
            output_filename = ''.join([reference_protein[0], '_', str(reference_protein[1]), naming_extension])
            output_filepath = ''.join([self.go_output_path, os.path.sep, output_filename, '-go-aggregate.dict'])
            if (os.path.exists(output_filepath) is True):
                print('GO information is already aggregated for the entry: ', output_filename,
                      '. If you did not expect this outcome please delete this file and execute again:\n',
                      output_filepath)
                continue
            file_name = ''.join([self.output_path, os.path.sep, reference_protein[0], '_', str(reference_protein[1]),
                                 naming_extension, '-merged-enriched.csv'])
            data = pd.read_csv(file_name, sep='\t')
            results = []
            pdbhandler = PDBHandler()
            full_pdb_file_path = os.path.join(pdb_dataset_path, '.'.join([self.reference_pdbid, 'pdb']))
            available_residues = pdbhandler.get_residue_range(full_pdb_file_path, self.reference_chain_id)
            reference_sequences = aligner.get_all_sequences(reference_protein, available_residues)
            self.alignment_levels = list(reference_sequences.keys())
            if(len(reference_sequences) < 1):
                print('No data on reference is available for analysis.\n')
                return False
            input_entries = [[row, reference_sequences, available_residues] for index, row in data.iterrows() if index > 0]
            if (self.debugging is True):
                for entry in input_entries:
                    results.append(self.align_candidate_by_property(entry))
            else:
                execution_handler = ExecutionHandler(self.cores, 30)
                results = execution_handler.parallelize(self.align_candidate_by_property, input_entries)
            go_alignment_info = defaultdict()
            alignment_info = defaultdict()
            for result_index, entry_alignment in enumerate(results):
                structured_result, raw_result = entry_alignment
                # Connecting the proteins' alignments with their Gene Ontology data for future queries
                for property_type in structured_result:
                    if (property_type not in go_alignment_info):
                        go_alignment_info[property_type] = defaultdict()
                    for go_id in structured_result[property_type]:
                        if (go_id not in go_alignment_info[property_type]):
                            go_alignment_info[property_type][go_id] = defaultdict()
                        for alignment_level in structured_result[property_type][go_id]:
                            if (alignment_level not in go_alignment_info[property_type][go_id]):
                                go_alignment_info[property_type][go_id][alignment_level] = defaultdict(list)
                            for key in structured_result[property_type][go_id][alignment_level].keys():
                                go_alignment_info[property_type][go_id][alignment_level][key].extend(structured_result[property_type][go_id][alignment_level][key])
                # Keeping the raw results for further inspection
                alignment_info[raw_result[0]] = raw_result[1]
            # Store GO-connected alignments
            with open(output_filepath, 'wb') as handle:
                pickle.dump(go_alignment_info, handle, protocol=pickle.HIGHEST_PROTOCOL)
            # Store alignments
            output_filepath = ''.join([self.go_output_path, os.path.sep, output_filename, '-alignments.dict'])
            with open(output_filepath, 'wb') as handle:
                pickle.dump(alignment_info, handle, protocol=pickle.HIGHEST_PROTOCOL)

    def get_available_go_information(self, entry):
        # Get coverage by the results for each GO type
        pdb_id, chain_id, naming_extension = entry
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        enricher.load_go_cache()
        if (len(naming_extension) > 0):
            if ('_' not in naming_extension):
                naming_extension = ''.join(['_', naming_extension])
        with open(os.path.join(self.go_output_path, ''.join([pdb_id, '_', chain_id, naming_extension, '-go-aggregate.dict'])), 'rb') as handle:
            alignment_info = pickle.load(handle)
        for propertyType in alignment_info:
            go_ids = alignment_info[propertyType]
            info = [(enricher.go_cache[go_id][0], len(go_ids[go_id][self.go_alignment_level]['identity'])) if self.go_alignment_level in go_ids[go_id] else None for go_id in go_ids]
            info = [x for x in info if x is not None]
            if(len(info) > 0):
                info = pd.DataFrame(info, columns=['propertyName', 'coverage']).sort_values(by='coverage', ascending=False)
                if (len(go_ids) > 0):
                    info.to_csv(os.path.join(self.go_output_path, ''.join(
                        [pdb_id, '_', chain_id, naming_extension, '_', propertyType, '-go-aggregate.csv'])), index=False,
                                sep='\t')
            else:
                print('No relevant data were found. Please try with different alignment level. Currently selected level is: ', self.go_alignment_level)

    def plot_correlations(self, entry):
        aggregated, average_span, characterized_positions, known_mask, figure_filename, available_residues = entry
        color_palette = sns.color_palette('colorblind', 4)

        # Plot matches with blue color
        plt.plot(range(1, len(aggregated) + 1), aggregated, color=color_palette[0], linewidth=0.5)
        data_frame = pd.DataFrame(aggregated, index=range(1, len(aggregated) + 1), columns=['matches'])

        # Compute exponentially weighted moving average (EWM)
        data_frame['ewm'] = data_frame.ewm(span=average_span).mean()
        known_areas_report = []
        if (known_mask is not None):
            data_frame['known'] = known_mask
        exp_moving_average = data_frame['ewm'].to_list()

        # Counter EWM lagging in the start and end
        for i in range(1, 3):
            exp_moving_average[i - 1] = data_frame.iloc[::-1]['matches'].ewm(span=average_span).mean().to_list()[
                len(data_frame) - i]
        exp_moving_average[len(exp_moving_average) - 1] = pd.DataFrame(
            data_frame.iloc[len(data_frame) - average_span + int(0.1 * average_span):len(data_frame) - 1][
                'matches'].to_list() + data_frame.iloc[0:int(0.1 * average_span)]['matches'].to_list(),
            columns=['matches'])['matches'].ewm(span=average_span).mean().to_list()[-1]
        data_frame['ewm'] = exp_moving_average
        # Plot EWM (yellow color)
        plt.plot(range(1, len(exp_moving_average) + 1), exp_moving_average, color=color_palette[1])

        # Select peaks according to their distance from EWM and their number of matches
        above_average = data_frame[data_frame['matches'] >= data_frame['ewm']].copy()
        above_average['avgdist'] = above_average['matches'] - above_average['ewm']
        above_average = above_average[
            above_average['avgdist'] > (np.mean(above_average['avgdist']) + np.std(above_average['avgdist']))]
        above_average['position'] = above_average.index
        close_to_max = data_frame[data_frame['matches'] >= (data_frame['matches'].max() * 0.9)].copy()
        close_to_max['position'] = close_to_max.index
        selected_positions = pd.concat([above_average, close_to_max])
        selected_positions.drop_duplicates(subset=['position'], keep='first', inplace=True)

        # Mark the selected peaks with green color
        plt.scatter(selected_positions.index.values.tolist(), selected_positions['matches'].to_list(),
                    color=color_palette[2])

        # Mark peaks belonging to areas with a known characteristic as specified by the user
        # (Orange dots)
        characterized_residues = []
        if (known_mask is not None):
            characterized_residues = selected_positions[selected_positions['known'] > 0].copy()
            plt.scatter(characterized_residues.index.values.tolist(), characterized_residues['matches'].to_list(),
                        color=color_palette[3])
        plt.xlabel('Residues')
        plt.ylabel('2D Matches')

        # Plot unavailable residues
        unavailable_residues = [x for x in range(1, len(aggregated) + 1) if x not in available_residues['fullRange']]
        plt.scatter(unavailable_residues, [0] * len(unavailable_residues), color='black', marker='|')

        # Create legends for the plot
        legend_elements = []
        legend_elements.append(
            Line2D([], [], marker='o', color='w', label='Selected', markerfacecolor=color_palette[2],
                   markeredgewidth=0.0, linewidth=0.0, markersize=10))
        if (known_mask is not None):
            legend_elements.append(
                Line2D([], [], marker='o', color='w', label='Marked', markerfacecolor=color_palette[3],
                       markeredgewidth=0.0, linewidth=0.0, markersize=10))
        legend_elements.append(plt.Rectangle((0, 0), 1, 1, label='EWMA', fc=color_palette[1]))
        if(len(unavailable_residues) > 0):
            legend_elements.append(plt.Rectangle((0, 0), 1, 1, label='Missing', fc='black'))
        plt.legend(handles=legend_elements, facecolor='white', framealpha=0.5,  loc='best', handlelength=1,
                   handleheight=1.125)

        # Store plots and related outputs
        plt.savefig(''.join([self.go_output_path, os.path.sep, figure_filename, '.png']), format='png', dpi=300)
        ax = plt.gca()
        ax.set_rasterized(True)
        if(self.vector_format_output is True):
            plt.savefig(''.join([self.go_output_path, os.path.sep, figure_filename, '.eps']), format='eps', dpi=300)
        if(self.tiff_format_output is True):
            plt.savefig(''.join([self.go_output_path, os.path.sep, figure_filename, '.tiff']), format='tiff', dpi=300)
        plt.clf()
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))
        if (known_mask is not None):
            known_areas_report.append(
                ''.join(['Total residues in provided known sites of reference protein: ', repr(characterized_positions)]))
            known_areas_report.append(
                ''.join([' Correlated residues in areas specified as known sites: ', repr(len(characterized_residues)), ' | ',
                         repr(round(len(characterized_residues) / characterized_positions * 100, 2)),
                         '%. \nResidue numbers:']))
            known_areas_report.append(','.join([repr(x) for x in characterized_residues.index.values.tolist()]))
            known_areas_report.append(
                ''.join(['Total correlated residues with this property: ', repr(len(selected_positions)), '\n',
                         repr(round(len(characterized_residues) / len(selected_positions) * 100, 2)),
                         '% of the correlated residues are in areas specified as known sites \nCorrelated residues for the selected property:']))
            known_areas_report.append(','.join([repr(x) for x in selected_positions.index.values.tolist()]))
        plt.clf()
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))
        return ''.join(known_areas_report)

    def find_batch_structure_matches(self, entry):
        pdb_id, chain_id, naming_extension, plot_output_path = entry
        pdbhandler = PDBHandler()
        pdb_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path, ''.join([pdb_id, '.pdb'])])
        available_residues = pdbhandler.get_residue_range(pdb_path, chain_id)
        pdbhandler.get_uniprot_accession_number(chain_id, pdb_id)
        evaluator = Evaluator(self.alignment_backend)
        evaluator.set_root_disk(self.root_disk)

        # Retrieve protein sequence length
        protein_sequence = evaluator.get_protein_sequence(pdbhandler.uniprot_accession_number)
        sequence_length = len(protein_sequence)
        aggregated_mask = np.array([0] * sequence_length)  # initialize mask

        reference_structure_id = ''.join([pdb_id, '_', chain_id, naming_extension])
        data_filepath = ''.join([self.go_output_path, os.path.sep, reference_structure_id, '-alignments.dict'])
        if os.path.exists(data_filepath) is True:

            with open(data_filepath, 'rb') as handle:
                alignment_info = pickle.load(handle)

            alignment_dict = {}

            # Collect alignments
            for structure_id in alignment_info:
                if structure_id not in alignment_dict:
                    if 'secondary' in alignment_info[structure_id]:
                        if alignment_info[structure_id]['secondary']['mask'] is not None:
                            alignment_frequencies = alignment_info[structure_id]['secondary']['mask']
                            padded_mask = [alignment_frequencies[available_residues['fullRange'].index(index + 1)]
                                           if (index + 1) in available_residues['fullRange']
                                              and available_residues['fullRange'].index(index + 1) < len(
                                           alignment_frequencies)  else 0 for index in range(sequence_length)]
                            aggregated_mask = aggregated_mask + padded_mask

            # Fetch unavailable residues
            pdb_path = os.path.join(self.root_disk, self.pdb_dataset_path, ''.join([pdb_id, '.pdb']))
            available_residues = pdbhandler.get_residue_range(pdb_path, chain_id)
            unavailable_residues = [x for x in range(1, len(aggregated_mask) + 1) if x not in available_residues['fullRange']]

            # Plot result
            presenter = Presenter()
            presenter.root_disk = self.root_disk
            presenter.verbose = self.debugging
            presenter.tiff_format_output = self.tiff_format_output
            presenter.vector_format_output = self.vector_format_output
            entry =  pdb_id, chain_id, naming_extension, aggregated_mask, unavailable_residues, plot_output_path
            presenter.plot_2D_matches(entry)

    def perform_go_meta_analysis(self, entry):
        pdb_id, chain_id, naming_extension, property_type, properties, search_terms = entry
        go_search = len(search_terms) != 0
        if(go_search):
            go_types = ['biologicalProcess', 'molecularFunction', 'cellularComponent']
        else:
            go_types = [property_type]
        enricher = Enricher()
        enricher.set_root_disk(self.root_disk)
        enricher.load_go_cache()
        pdbhandler = PDBHandler()
        evaluator = Evaluator(self.alignment_backend)
        evaluator.set_root_disk(self.root_disk)
        encountered = []
        alignment_frequencies = defaultdict()
        pdb_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path, ''.join([pdb_id, '.pdb'])])
        available_residues = pdbhandler.get_residue_range(pdb_path, chain_id)
        pdbhandler.get_uniprot_accession_number(chain_id, pdb_id if self.override_pdb_id == '' else self.override_pdb_id)

        # Retrieve protein sequence length
        protein_sequence = evaluator.get_protein_sequence(pdbhandler.uniprot_accession_number)
        sequence_length = len(protein_sequence)
        average_span = int(0.3 * sequence_length)
        aggregated = np.array([0] * sequence_length)  # initialize mask

        # Set mask for residue positions with known properties
        known_mask = None
        known_sites = []
        if (len(self.known_areas)>0):
            known_sites = sorted(list(set([x for y in self.known_areas for x in y])))
            known_mask = [1 if x in known_sites else 0 for x in range(1, sequence_length + 1)]

        if (len(naming_extension) > 0):
            if ('_' not in naming_extension):
                naming_extension = ''.join(['_', naming_extension])
        with open(os.path.join(self.go_output_path, ''.join([pdb_id, '_', chain_id, naming_extension, '-go-aggregate.dict'])), 'rb') as handle:
            alignment_info = pickle.load(handle)
        analysis_details_output = []
        checked_structure_ids = []
        go_analysis_filename = ''
        # Search stored GO term related alignments of proteins in the results
        # for a user-specified property of interest
        for property_type in go_types:
            total = len(alignment_info[property_type].keys())
            analysis_details_output.append(''.join(['###########', property_type, '\n']))
            for go_id in alignment_info[property_type]:
                total -= 1
                # if the cache does not include this GO id
                # fetch its associated term
                if (go_id not in enricher.go_cache):
                    term_info = enricher.get_go_term(go_id)
                    enricher.update_go_cache(go_id, term_info)
                search_condition = False
                # Compare user defined partial/full GO terms
                # with current GO term
                if go_search:
                    for search_term in search_terms:
                        if search_term in enricher.go_cache[go_id][0]:
                            search_condition = True
                            break
                if (enricher.go_cache[go_id][0] in properties or search_condition):
                    structure_ids = alignment_info[property_type][go_id][self.go_alignment_level]['structure_id']
                    excluded_indices = [id_index for id_index, structure_id in enumerate(structure_ids) if structure_id in checked_structure_ids]
                    checked_structure_ids.extend([structure_id for id_index, structure_id in enumerate(structure_ids) if id_index not in excluded_indices])
                    encountered.append(enricher.go_cache[go_id][0])
                    if (len(excluded_indices) != 0):
                        alignment_frequencies[go_id] = np.sum(np.array(list(filter(None, [alignment_mask for mask_index, alignment_mask
                                                              in enumerate(alignment_info[property_type][go_id][self.go_alignment_level]['mask'])
                                                              if mask_index not in excluded_indices]))), axis=0)
                    else:
                        alignment_frequencies[go_id] = np.sum(np.array(list(filter(None, alignment_info[property_type][go_id][self.go_alignment_level]['mask']))), axis=0)
                    # Log related information for further inspection
                    stats_output = 'Aggregated alignment statistics: \n'
                    for alignment_level in alignment_info[property_type][go_id]:
                        stats_output = ''.join([stats_output, ' Median ',
                                                alignment_level.replace('secondary', 'secondary structure'),
                                                ' sequence identity: ',
                               repr(round(np.median(np.array(alignment_info[property_type][go_id][alignment_level]['identity'])), 4)), '\n'])
                    if isinstance(alignment_frequencies[go_id], list) is False:
                        if(len(alignment_frequencies[go_id].shape) != 0):
                            padded_mask = [alignment_frequencies[go_id][available_residues['fullRange'].index(index + 1)]
                                           if (index + 1) in available_residues['fullRange']
                                           and available_residues['fullRange'].index(index + 1) < len(alignment_frequencies[go_id])
                                           else 0
                                           for index in range(sequence_length)]
                            assertion_index = random.randint(0, len(available_residues['fullRange']) - 1)
                            assert padded_mask[available_residues['fullRange'][assertion_index] - 1] == alignment_frequencies[go_id][
                                assertion_index]
                            aggregated = aggregated + padded_mask
                            analysis_details_output.append(''.join(['********* ', enricher.go_cache[go_id][0], '\n\n',
                                                                    stats_output, '\n\n\nAlignment information: \n']))
                            for alignment_level in self.alignment_levels:
                                analysis_details_output.append(''.join(['\n****** ', alignment_level, '\n\n\n',
                                '\n\n\n'.join([x for x in alignment_info[property_type][go_id][alignment_level]['alignment'] if len(x) > 0]), '\n\n']))
                # Terminate search if every property of interest is found or
                # there are not more available properties to search for
                if ((go_search is False and len(encountered) == len(properties)) or total == 0):
                    if(go_search is True):
                        if(go_analysis_filename == ''):
                            go_analysis_filename = ''.join([pdb_id, '_', chain_id, naming_extension, '_',
                                                            '_'.join([search_term.replace(os.path.sep, '')
                                                                      for search_term in search_terms]), '_',
                                                            uuid.uuid4().hex])
                    else:
                        go_analysis_filename = ''.join([pdb_id, '_', chain_id, naming_extension, '_', properties[0].replace(' ', '_'),  '_', uuid.uuid4().hex])
                    entry = aggregated, average_span, len(known_sites), known_mask, go_analysis_filename, available_residues
                    if(len(encountered) > 0):
                        print(''.join(['Chosen properties found among the candidates: ', repr(encountered)]))
                    else:
                        print('No alignments were found for this property.')
                    break
        # Log aggregated search information
        if (np.sum(aggregated) > 0):
            known_sites_report = self.plot_correlations(entry)
            if (len(known_sites_report) > 0):
                analysis_details_output.append(''.join(['\n\nKnown sites report: \n', known_sites_report]))
            analysis_details_output.append(''.join(['\n\nChosen properties found among the candidates: ',  ','.join(encountered)]))
            analysis_details_output.append(''.join(['\n\nThe candidates: ', ','.join(checked_structure_ids)]))

            # Create mini report
            output_filepath = ''.join([self.output_path, os.path.sep, pdb_id, '_', chain_id, naming_extension, '-merged-enriched_eval_full.csv'])
            results = pd.read_csv(output_filepath, sep='\t')
            filter_mask = None
            for struct_id in checked_structure_ids:
                pdbid, chainid = struct_id.split('_')
                if(filter_mask is None):
                    filter_mask = (results['pdbId'] == pdbid) & (results['chainId'] == chainid)
                else:
                    filter_mask = filter_mask | ((results['pdbId'] == pdbid) & (results['chainId'] == chainid))
            results = results[filter_mask]
            report_path = os.path.join(self.go_output_path, ''.join([go_analysis_filename, '-pres.csv']))
            results.to_csv(report_path, index=False, sep='\t')
            presenter = Presenter()
            presenter.root_disk = self.root_disk
            presenter.verbose = self.debugging
            presenter.tiff_format_output = self.tiff_format_output
            presenter.vector_format_output = self.vector_format_output
            presenter.create_report([pdb_id, chain_id, naming_extension], file_name=report_path)

            # Output protein-protein interactions
            if os.path.exists(os.path.join(self.go_output_path, ''.join([go_analysis_filename, '_common_interactions.tsv']))) is False:
                self.retrieve_ppi(pdbhandler.uniprot_accession_number, results['uniprotId'].to_list(), go_analysis_filename)

        with open(os.path.join(self.go_output_path, ''.join([go_analysis_filename, '.log'])), 'w') as logFile:
            logFile.write('\n'.join(analysis_details_output))

    def retrieve_ppi(self, reference_uniprot_id, candidate_uniprot_ids, go_analysis_filename):
        input_entries = [[reference_uniprot_id, candidate_uniprot_id] for candidate_uniprot_id in candidate_uniprot_ids]
        results = []

        # Collect the interactions
        if (self.debugging is False):
            execution_handler = ExecutionHandler(self.cores, 12 * 60)
            results = execution_handler.parallelize(self.fetch_common_interactors, input_entries)
        else:
            for entry in input_entries:
                results.append(self.fetch_common_interactors(entry))

        # Convert the output into a pandas DataFrame
        data = [line.split("\t") for output in results for line in output.split('\n') if len(line)> 0]
        # if there are any interactors
        if len(data) > 0:
            df = pd.DataFrame(data, columns=["ID", "ProteinA", "ProteinB"])

            # Make sure the main proteins are in the second column
            main_proteins = [reference_uniprot_id] + candidate_uniprot_ids
            df["Protein"] = df.apply(
                lambda row: row["ProteinB"] if row["ProteinB"] in main_proteins else row["ProteinA"], axis=1)
            df["CommonInteractor"] = df.apply(
                lambda row: row["ProteinA"] if row["Protein"] == row["ProteinB"] else row["ProteinB"], axis=1)

            # Drop the original columns and rearrange
            df = df[["ID", "Protein", "CommonInteractor"]]

            # Sort by the CommonInteractor
            df = df.sort_values(by=["CommonInteractor", "Protein"])

            df.to_csv(os.path.join(self.go_output_path, ''.join([go_analysis_filename, '_common_interactions.tsv'])), sep="\t", index=False)

            df = df.drop_duplicates(subset=['Protein', 'CommonInteractor'], keep='first')
            df = df[df['Protein'] != reference_uniprot_id]


            # Create a dictionary that counts appearances of each common interactor
            interactor_count = df['CommonInteractor'].value_counts().to_dict()
            sorted_interactor_count = dict(sorted(interactor_count.items(), key=lambda item: item[1], reverse=True))

            # Save the dictionary to a TSV file
            with open(os.path.join(self.go_output_path, ''.join([go_analysis_filename, '_interactor_count.tsv'])), 'w') as file:
                file.write('CommonInteractor\tCount\n')
                for k, v in sorted_interactor_count.items():
                    file.write(f"{k}\t{v}\n")

    def fetch_common_interactors(self, entry):
        reference_uniprot_id, candidate_uniprot_id = entry
        biogrid_data_path = os.path.join(self.root_disk, 'uniprot-interactors.tsv')
        interactions = ''

        # if there is a BioGrid dataset available
        if os.path.exists(biogrid_data_path):
            # Retrieve common interactors
            try:
                # Construct the command
                awk_script = f"""awk '
                  NR==FNR {{ 
                    if ($2=="{reference_uniprot_id}" && $3!="{candidate_uniprot_id}") proteins[$3]++; 
                    else if ($3=="{reference_uniprot_id}" && $2!="{candidate_uniprot_id}") proteins[$2]++;
                    else if ($2=="{candidate_uniprot_id}" && $3!="{reference_uniprot_id}") proteins2[$3]++;
                    else if ($3=="{candidate_uniprot_id}" && $2!="{reference_uniprot_id}") proteins2[$2]++;
                    next; 
                  }} 
                  ($2 in proteins && $2 in proteins2 && ($3=="{reference_uniprot_id}" || $3=="{candidate_uniprot_id}")) || 
                  ($3 in proteins && $3 in proteins2 && ($2=="{reference_uniprot_id}" || $2=="{candidate_uniprot_id}"))
                ' {biogrid_data_path} {biogrid_data_path} """
                command = f"timeout 80s {awk_script}"
                if (self.verbose is True):
                    print(' '.join(command))
                output = subprocess.run(command, shell=True, capture_output=True)
                interactions = output.stdout.decode("utf-8")
            except:
                if (self.verbose is True):
                    print(traceback.format_exc())
        return interactions
