from src.pdbhandler import PDBHandler
import numpy as np
import pandas as pd
import os
import time
import traceback
import io
from src.scanner import Scanner
from src.structurealigner import StructureAligner
from src.bhattacharyyadistance import BhattacharyyaDistance
from scipy.stats import wasserstein_distance
from src.executionhandler import ExecutionHandler

# Constrained mode of the comparison method (segment comparisons)

class SegmentScanner(Scanner):

    def __init__(self):
        super().__init__()
        self.minimum_residue_range = 3
        self._column_headers = ['structureId', 'b-phipsi', 'w-rdist', 't-alpha', 'alignedRange',
                                'referenceSegmentIndex', 'referenceRange', 'alignedSequences']
        self.selected_alignment_level = 'mixed'
        self.alignment_levels = ['primary', 'secondary', 'mixed', 'hydrophobicity']
        self.alignment_backend = ''
        self.segments = []

    # Not applicable for this use case: features are extracted on the fly
    def load_features(self, pdb_id, chain_id, metric_index, **kwargs):
        print('Site scanning does not cache data.')
        raise Exception('This method is not available.')

    # Not applicable for this use case: features are extracted on the fly
    def collect_features(self, cores=10):
        print('Site scanning does not cache data.')
        raise Exception('This method is not available.')

    # Not applicable for this use case: features are extracted on the fly
    def create_feature_filelist(self):
        print('Site scanning does not cache data.')
        raise Exception('This method is not available.')

    def compute_metrics(self, entry):
        reference_angle_data, reference_distances, reference_triangle_metric, reference_structure, segments, pdb_path = entry
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        full_pdb_path = os.path.sep.join([pdb_dataset_path, pdb_path])
        pdb_id = pdb_path.split('.')[0]
        if (self.pdb_validation is False):
            chains = pdbhandler.get_protein_chain_ids(full_pdb_path)
        else:
            chains = pdbhandler.fetch_pdb_peptidic_chains(full_pdb_path)
        results = []
        if chains is not False and len(chains) > 0:
            for chainID in chains:
                bphipsi = None
                wrdist = None
                talpha = None
                try:
                    # Extract features from the specified segments
                    result = self.extract_features([segments, reference_structure, pdb_path, chainID])
                    aligned_sequences, residue_range_selection, reference_segments, reference_segment_alignments, features = result
                    angle_data, residue_distances, triangle_metric = features
                    # Compute metrics for the specified segments
                    if (len(angle_data) > 1):
                        bphipsi = BhattacharyyaDistance.multivariate_compare(reference_angle_data, angle_data)

                    if (len(residue_distances) > 1):
                        wrdist = wasserstein_distance(reference_distances, residue_distances)

                    if (triangle_metric > 0 and reference_triangle_metric > 0):
                        talpha = np.exp(abs(reference_triangle_metric - triangle_metric)) - 1
                except:
                    raise
                    residue_range_selection = []
                    reference_segments = []
                    reference_segment_alignments = []
                    aligned_sequences = []
                    print(''.join([pdb_id, '|', chainID, '\n', traceback.format_exc(), '\n']))
                results.append(['_'.join([pdb_id, chainID]), bphipsi, wrdist, talpha, repr(residue_range_selection),
                                repr(reference_segments),
                                repr(reference_segment_alignments), repr(aligned_sequences)])
        return results

    def extract_features_from_segment(self, entry):
        pdb_path, segment, site_index, chain_id, available_residues = entry
        angle_data = []
        residue_distances = []
        triangle_metric = -1
        path_parts = pdb_path.split(os.path.sep)
        [pdb_id, _] = path_parts[-1].split('.')
        pdbhandler = PDBHandler(pdb_id)
        pdbhandler.root_disk = self.root_disk
        # Find contiguous areas of the segment and discard spatial outliers
        print("Clustering the residue positions of preset segment.")
        parts = pdbhandler.process_segment(segment, available_residues)
        participating_residues = []
        output_result_path = ''.join([self.metrics_output_path, os.path.sep, self.reference_pdb_id, '_',
                                      self.reference_chain_id, '_site', repr(site_index), '-parts.csv'])
        pd.DataFrame(parts.items(), columns=['partIndex', 'residues']).sort_values(by='partIndex').to_csv(output_result_path, index=False)
        # For each contiguous part, extract features
        for residue_range in parts:
            residue_selection = list(range(min(parts[residue_range]), max(parts[residue_range])))
            if (len(residue_selection) < self.minimum_residue_range):
                continue
            participating_residues.extend(residue_selection)
            angles, status = self.extract_feature(chain_id, '', pdb_path, 'angles', residue_selection=residue_selection, raw_data=True)
            if (status is not False):
                angle_data.extend(angles)
            else:
                raise Exception('Bad reference data (angles).')
            distances, status = self.extract_feature(chain_id, '', pdb_path, 'distances', residue_selection=residue_selection)
            if (status is not False):
                residue_distances.extend(list(distances))
            else:
                raise Exception('Bad reference data (distances).')
        if (len(participating_residues) > 0):
            pdbhandler.residue_selection = participating_residues
            triangle_metric, status = self.extract_feature(chain_id, '', pdb_path, 'triangles', residue_selection=participating_residues)
            if (status is False):
                raise Exception('Bad reference data (triangles).')
        return parts, np.array(angle_data), np.array(residue_distances), triangle_metric

    def remove_overlaps(self, ranges):
        # Combine overlapping residue position ranges
        ranges.sort()
        while True:
            combined = []
            overlapping = []
            for index in range(1, len(ranges)):
                if (ranges[index - 1][1] > ranges[index][0]):
                    overlapping.append(index - 1)
                    overlapping.append(index)
                    combined.append(
                        [min(ranges[index - 1][0], ranges[index][0]), max(ranges[index - 1][1], ranges[index][1])])
            overlapping = list(set(overlapping))
            for index in range(len(ranges)):
                if (index not in overlapping):
                    if (ranges[index][1] - ranges[index][0] > 1):
                        if (ranges[index] not in combined):
                            combined.append(ranges[index])
            ranges = combined
            if (len(overlapping) == 0):
                break
        return ranges

    def extract_features(self, entry):
        reference_segment_parts, reference_structure, pdb_path, chain_id = entry
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        full_pdb_file_path = os.path.sep.join([pdb_dataset_path, pdb_path])
        path_parts = pdb_path.split(os.path.sep)
        [pdb_id, _] = path_parts[-1].split('.')
        pdbhandler = PDBHandler(pdb_id)
        pdbhandler.root_disk = self.root_disk
        aligner = StructureAligner()
        aligner.set_root_disk(self.root_disk)
        aligner.pdb_dataset_path = pdb_dataset_path
        angle_data = []
        residue_range_selection = []
        aligned_range_selection = []
        residue_distances = []
        triangle_metric = -1
        participating_residues = []
        aligned_reference_segment_indices = []
        reference_segments_alignments = []
        available_residues = pdbhandler.get_residue_range(full_pdb_file_path, chain_id)
        candidate_sequences = aligner.get_all_sequences((pdbhandler.structure_id, chain_id), available_residues)
        aligned_sequences = []
        output = []
        # Find parts of a candidate protein that match a part of the reference segment
        # in a level specified by the configuration (e.g. protein sequence)
        if (len(candidate_sequences) > 0
            and self.selected_alignment_level in candidate_sequences
            and candidate_sequences[self.selected_alignment_level] != ''):
            total_residue_range = [x for x in range(available_residues['fullRange'][0], available_residues['fullRange'][-1] + 1)]
            for label in reference_segment_parts:
                reference_start = reference_segment_parts[label][0] - reference_structure['residues'][0]
                reference_end = reference_segment_parts[label][-1] - reference_structure['residues'][0] + 1
                reference_part = reference_structure[self.selected_alignment_level][reference_start:reference_end]
                gap_substitutes = ['{', ' ']
                score, alignment, alignment_result = aligner.local_align(reference_part, candidate_sequences[
                    self.selected_alignment_level], gap_substitutes)
                if (score is False):
                    continue
                gaps = alignment_result[2].count('-')
                aligned_range, actual_range, range_length = aligner.get_alignment_positions(alignment, candidate_sequences[self.selected_alignment_level],
                                                                                                  total_residue_range)
                # Check whether the validity of the aligned positions that are extracted from the alignment
                if(aligner.backend == 'parasail'):
                    assert alignment_result[2].replace('-', '') == \
                           candidate_sequences[self.selected_alignment_level][aligned_range[0]:aligned_range[1]].replace('-', gap_substitutes[1])
                if(len(aligned_range) < 2):
                    continue
                aligned_parts = []
                for alignment_level in self.alignment_levels:
                    if alignment_level in candidate_sequences:
                        aligned_parts.append(candidate_sequences[alignment_level][aligned_range[0]:aligned_range[1]])
                if ((range_length - gaps) >= self.minimum_residue_range):
                    residue_range_selection.append(actual_range)
                    aligned_reference_segment_indices.append(label)
                    # Retrieve aligned reference parts for logging
                    reference_relative_range, reference_gaps = aligner.get_reference_sequence_alignment(alignment_result)
                    if (aligner.backend == 'parasail'):
                        assert alignment_result[0].replace('-', '') == \
                           reference_structure[self.selected_alignment_level][reference_start+reference_relative_range[0]:
                                                                              reference_start+reference_relative_range[1]].replace('-', gap_substitutes[0])
                    for alignment_level in self.alignment_levels:
                        if alignment_level in candidate_sequences:
                            aligned_parts.append(reference_structure[alignment_level][reference_start+reference_relative_range[0]:reference_start+reference_relative_range[1]])
                    reference_residue_range = [x for x in range(reference_structure['residues'][0],
                                                                reference_structure['residues'][-1] + 1)][
                                              reference_start:reference_end]
                    reference_residue_range = reference_residue_range[
                                              reference_relative_range[0]:reference_relative_range[1]]
                    # Sanity check of alignment's details
                    candidate_part_length = range_length + gaps
                    reference_part_length = len(reference_residue_range) + reference_gaps
                    assert gaps >= 0 and reference_gaps >= 0 and (candidate_part_length == reference_part_length), \
                        ''.join(['different total alignment lengths: ', repr(candidate_part_length), ' | ', repr(reference_part_length)])
                    reference_segments_alignments.append([reference_residue_range[0], reference_residue_range[-1]])
                    aligned_sequences.append(aligned_parts)
            # Extract features and compute metrics between the aligned parts of the compared proteins
            if (len(residue_range_selection) > 0):
                aligned_range_selection = residue_range_selection.copy()
                residue_range_selection = self.remove_overlaps(residue_range_selection)
                status = True
                for residue_range in residue_range_selection:
                    residue_selection = list(range(residue_range[0], residue_range[1]))
                    if (len(residue_selection) < self.minimum_residue_range):
                        continue
                    participating_residues.extend(residue_selection)
                    angles, status = self.extract_feature(chain_id, '', full_pdb_file_path, 'angles',
                                                  residue_selection=residue_selection, raw_data=True)
                    if (status is not False):
                        angle_data.extend(angles)
                    else:
                        break
                    distances, status = self.extract_feature(chain_id, '', full_pdb_file_path, 'distances',
                                                     residue_selection=residue_selection)
                    if (status is not False):
                        residue_distances.extend(list(distances))
                    else:
                        break
                if (status is False):
                    angle_data = []
                    residue_distances = []
                elif (len(participating_residues) > 0):
                    pdbhandler.residue_selection = participating_residues
                    triangle_metric, status = self.extract_feature(chain_id, '', full_pdb_file_path, 'triangles',
                                                           residue_selection=participating_residues)
        # Construct output with full information on the alignments
        output.append(aligned_sequences)
        output.append(aligned_range_selection)
        output.append(aligned_reference_segment_indices)
        output.append(reference_segments_alignments)
        output.append([np.array(angle_data), np.array(residue_distances), triangle_metric])
        return output

    def scan_candidates(self):
        if (os.path.exists(self.metrics_output_path) is False):
            os.makedirs(self.metrics_output_path)
        aligner = StructureAligner()
        pdbhandler = PDBHandler(self.reference_pdb_id)
        pdbhandler.root_disk = self.root_disk
        aligner.set_root_disk(self.root_disk)
        aligner.pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        reference_pdb_path = os.path.sep.join(
            [self.root_disk, self.pdb_dataset_path, ''.join([self.reference_pdb_id, '.pdb'])])
        reference_pdb_path = pdbhandler.handle_alphafold_pdbid(reference_pdb_path)
        available_residues = pdbhandler.get_residue_range(reference_pdb_path, self.reference_chain_id)
        reference_structure = aligner.get_all_sequences((pdbhandler.structure_id, self.reference_chain_id),
                                                        available_residues)
        if(self.selected_alignment_level not in reference_structure):
            print('Your selected alignment level (', self.selected_alignment_level, ') is not available for your reference protein. Please select a different one.')
            exit(1)
        self.alignment_levels = list(reference_structure.keys())
        reference_structure['residues'] = available_residues['fullRange']
        # For each reference segment in the configuration
        for site_index, segment in enumerate(self.segments):
            print(''.join(['Scanning for segment ', repr(site_index)]))
            output_result_path = ''.join([self.metrics_output_path, os.path.sep, self.reference_pdb_id, '_', self.reference_chain_id, '_site',
                                          repr(site_index), '-metrics.csv'])
            if (os.path.exists(output_result_path) is True):
                print(''.join(['Metrics from are already computed. If you wish to recompute them delete the following file and execute the method again:\n',
                               output_result_path]))
                continue
            entry = reference_pdb_path, segment, site_index, self.reference_chain_id, reference_structure['residues']
            reference_segment_parts, reference_angles, reference_distances, reference_triangle_metric = self.extract_features_from_segment(entry)
            entries = [(reference_angles, reference_distances, reference_triangle_metric, reference_structure,
                        reference_segment_parts, pdbFilename) for pdbFilename in self.candidates]
            final_results = []
            # Compute metrics in parallel
            print('Scanning candidates')
            if (self.debugging is False):
                execution_handler = ExecutionHandler(self.cores, 10 * 60)
                results = execution_handler.parallelize(self.compute_metrics, entries)
            else:
                results = []
                for candidate in entries:
                    result = self.compute_metrics(candidate)
                    results.append(result)
            # Filter out empty results
            for chainsResults in results:
                if isinstance(chainsResults, list):
                    for chainResult in chainsResults:
                        final_results.append(chainResult)
            # Construct output
            output_result_path = ''.join([self.metrics_output_path, os.path.sep, self.reference_pdb_id, '_',
                                          self.reference_chain_id, '_site', repr(site_index), '-metrics.csv'])
            output_string = '\n'.join(['\t'.join([repr(x) if isinstance(x, np.floating) or isinstance(x, float) else x for x in result]) \
                            for result in np.array(final_results) if None not in result])
            data = pd.read_csv(io.StringIO(output_string), names=self._column_headers, sep='\t')
            data.sort_values(by=['b-phipsi', 'w-rdist', 't-alpha'],
                             ascending=[True, True, True]).to_csv(output_result_path, index=False)

    def set_candidates(self):
        # Include all PDB files in the dataset path into the candidate proteins set
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        self.candidates = [filename for filename in os.listdir(pdb_dataset_path) if
                           os.path.isfile(os.path.join(pdb_dataset_path, filename))]
