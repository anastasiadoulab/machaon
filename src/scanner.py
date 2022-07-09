from collections import defaultdict
from src.pdbhandler import PDBHandler
from src.bhattacharyyadistance import BhattacharyyaDistance
import numpy as np
import pandas as pd
import os
import time
import traceback
import pickle
from scipy.stats import wasserstein_distance
import io
import subprocess
import math
from src.executionhandler import ExecutionHandler

# Standard mode of the comparison method (whole structure comparisons)

class Scanner:

    def __init__(self):
        self.root_disk = ''
        self.debugging = False
        self.verbose = False
        self.reference_pdb_id = ''
        self.reference_chain_id = ''
        self.chunk_size = 100000
        self.candidates = defaultdict
        self.features_path = ''
        self.features_file_list = ''
        self.pdb_dataset_path = ''
        self.metrics_output_path = ''
        self._feature_file_suffixes = ['angles', 'distances', 'triangles']
        self._feature_file_suffixes_dict = {k: v for v, k in enumerate(self._feature_file_suffixes)}
        self._column_headers = ['b-phipsi', 'w-rdist', 't-alpha']
        self._feature_filelist_header = '\t'.join(['structureId', 'chainId', 'metric'])
        self.cores = 1
        self.update_features = []
        self.recompute_features = False
        self.extend_feature_data = False
        self.pdb_validation = False

    # Load pre-computed features for a protein required for the calculation of a specified metric
    def load_features(self, pdb_id, chain_id, metric_index, naming_extension=''):
        store_data_path = os.path.join(self.root_disk, self.features_path)
        naming_extension = ''.join([naming_extension, '_']) if naming_extension != '' else ''
        features_data_path = os.path.join(store_data_path, ''.join(
            [pdb_id.replace('_', '-'), '_', chain_id, '_', naming_extension, self._feature_file_suffixes[metric_index], '.pkl']))
        metric_data = False
        if (os.path.exists(features_data_path)):
            with open(features_data_path, 'rb') as handle:
                metric_data = pickle.load(handle)
            if (metric_index == 0):
                metric_data = np.array(metric_data['phiPsiPairs'])
        return metric_data

    def compute_metrics(self, entry):
        reference_metric_data, comparison_info = entry
        naming_extension = ''
        if (len(comparison_info) < 4):
            pdb_id, chain_id, metric_index = comparison_info
        else:
            pdb_id, chain_id, metric_index, naming_extension = comparison_info
        # Load the required features
        metric_data = self.load_features(*comparison_info)
        result = False
        # Compute the metric
        if (metric_data is not False and metric_data is not None):
            try:
                if (metric_index == 0):
                    if (len(metric_data) > 2):
                        result = BhattacharyyaDistance.multivariate_compare(reference_metric_data, metric_data)
                elif (metric_index == 1):
                    result = math.log10(wasserstein_distance(reference_metric_data, metric_data) + 1)
                elif (metric_index == 2):
                    result = np.exp(abs(reference_metric_data - metric_data)) - 1
            except:
                print(''.join([pdb_id, '|', chain_id, '\n', traceback.format_exc()]))
        # Construct output
        if (result is not False):
            output = ['_'.join([pdb_id, chain_id])]
            if (len(comparison_info) > 3):
                output.append(naming_extension)
            output.append(result)
            result = output

        return result

    def collect_features(self):
        # Check if filelist exists
        store_path = os.path.sep.join([self.root_disk, ''.join([self.features_path, '.csv'])])
        store_path_exists = os.path.exists(store_path)
        # Re-compute the features by user's choice
        if (self.recompute_features is True):
            self.update_features = self._feature_file_suffixes.copy()
        # Extract features from the PDB dataset that was set in the configuration
        if (store_path_exists is False or len(self.update_features) > 0 or self.extend_feature_data is True):
            os.makedirs(store_path.replace('.csv', ''), exist_ok=True)
            pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
            pdb_filenames = [filename for filename in os.listdir(pdb_dataset_path) if
                            os.path.isfile(os.path.join(pdb_dataset_path, filename))]
            input_set = pdb_filenames
            if (self.debugging is False):
                # Extract features in parallel
                execution_handler = ExecutionHandler(self.cores, 12 * 60)
                execution_handler.parallelize(self.extract_features, input_set)
            else:
                for input_data in input_set:
                    self.extract_features(input_data)
            self.create_feature_filelist()
        else:
            print('Features have already been computed. Please check if there is a directory misconfiguration if you were not expecting this outcome and execute again.')

    def create_feature_filelist(self):
        # Create a file listing all filenames of the stored features for each PDB chain
        try:
            output_dir = os.path.sep.join([self.root_disk, self.features_path])
            output = subprocess.run(['ls', output_dir], capture_output=True,
                                    timeout=30)
            lines = output.stdout.decode("utf-8").split('\n')
            data_info = [self._feature_filelist_header]
            for line in lines:
                if (len(line) > 2):
                    parts = line.split('_')
                    suffix = parts[-1].replace('.pkl', '')
                    if (suffix not in self._feature_file_suffixes_dict):
                        continue
                    output_line = [parts[0], parts[1], str(self._feature_file_suffixes_dict[suffix])]
                    # if there is an extended naming (e.g. domain name)
                    if (len(parts) > 3):
                        output_line.append(parts[2])
                    data_info.append('\t'.join(output_line))
            with open(''.join([output_dir, '.csv']), 'w') as datalist:
                datalist.write('\n'.join(data_info))
        except:
            print(traceback.format_exc())

    def extract_feature(self, chain_id, output_path, full_pdbfile_path, feature, **kwargs):
        # Extract features for a PDB chain
        feature_data = False
        condition = False
        try:
            pdbhandler = PDBHandler()
            pdbhandler.root_disk = self.root_disk
            pdbhandler.verbose = self.verbose
            if 'residueSelection' in kwargs:
                pdbhandler.residue_selection = kwargs['residueSelection']
            if (os.path.exists(output_path) is False or feature in self.update_features):
                if (feature == 'angles'):
                    pdbhandler.fetch_angles(full_pdbfile_path, chain_id)
                    feature_data = pdbhandler.chain_angles[chain_id]
                    condition = len(feature_data['psi']) > 2 and len(feature_data['phi']) > 2
                elif (feature == 'distances'):
                    pdbhandler.calculate_residue_distances(full_pdbfile_path, chain_id)
                    feature_data = pdbhandler.residue_distances
                    condition = len(pdbhandler.residue_distances) != 0
                elif (feature == 'triangles'):
                    pdbhandler.load_points(full_pdbfile_path, chain_id)
                    if (len(pdbhandler.points) > 2):
                        triangles = pdbhandler.get_mesh_triangles(pdbhandler.points)
                        if (triangles is not None):
                            feature_data = np.log(np.asarray(triangles).shape[0])
                    condition = feature_data is not False
                if (condition):
                    if ('noStore' not in kwargs):
                        with open(output_path, 'wb') as handle:
                            pickle.dump(feature_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
                else:
                    feature_data = False
        except:
            feature_data = False
            print(traceback.format_exc())
        if (feature_data is False):
            with open(os.path.join('logs', ''.join([time.strftime('%Y-%m-%d', time.localtime()), '.log'])),
                      'a', encoding='utf-8') as fail_log_file:
                fail_log_file.write(
                    ' '.join(['feature_extraction', chain_id, output_path, full_pdbfile_path, feature, '\n']))
            if(self.verbose is True):
                print(' '.join(
                    [full_pdbfile_path.split('.')[-2].split(os.path.sep)[-1], chain_id, ': ', feature, ' were not extracted [',
                     output_path, ']']))
        return feature_data

    def extract_features(self, entry):
        # Extract features for each chain that represents a protein in the PDB file
        result = True
        store_data_path = os.path.sep.join([self.root_disk, self.features_path])
        pdb_path = entry
        pdbhandler = PDBHandler()
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        full_pdbfile_path = os.path.sep.join([pdb_dataset_path, pdb_path])

        if(self.pdb_validation is False):
            chains = pdbhandler.get_protein_chain_ids(full_pdbfile_path)
        else:
            chains = pdbhandler.fetch_pdb_peptidic_chains(full_pdbfile_path)
        if chains is False or len(chains) == 0:
            result = False
            if(self.verbose is True):
                print(pdb_path, ' | This file will be ignored.')
        else:
            # Extract each required feature (angles, distances, alpha shape triangles) from each PDB chain
            for chain_id in chains:
                if(len(chain_id.strip()) == 0):
                    if (self.verbose is True):
                        print(pdb_path, 'contains an empty chain identifier. This chain will be ignored.')
                    continue
                output_path = os.path.sep.join([store_data_path, pdb_path.replace('_', '-').replace('.pdb', ''.join(['_', chain_id]))])
                for feature_index, feature in enumerate(self._feature_file_suffixes):
                    feature_output_path = ''.join([output_path, '_', self._feature_file_suffixes[feature_index], '.pkl'])
                    # Skip if this feature is not set to be recomputed when recomputing is enabled
                    if (len(self.update_features) > 0 and feature not in self.update_features):
                        continue
                    # New extractions only
                    if(self.extend_feature_data is True and os.path.exists(feature_output_path) is True):
                        continue
                    self.extract_feature(chain_id, feature_output_path, full_pdbfile_path, feature)
        return result

    def scan_candidates(self):
        # Compute metrics in parallel for each candidate protein in the dataset
        if (os.path.exists(self.metrics_output_path) is False):
            os.makedirs(self.metrics_output_path)
        for metric_index, metric in enumerate(self._feature_file_suffixes):
            print('*', self._column_headers[metric_index])
            output_result_path = ''.join([self.metrics_output_path, os.path.sep, self.reference_pdb_id,
                                          '_', self.reference_chain_id, '_', self._column_headers[metric_index], '.csv'])
            if (os.path.exists(output_result_path) is True):
                print(''.join([self._column_headers[metric_index], ' is already computed. If you wish to recompute it, delete',
                               ' the following file and execute the method again:\n', output_result_path]))
                continue
            final_results = []
            reference_metric_data = self.load_features(self.reference_pdb_id, self.reference_chain_id, metric_index)
            # Check if there are reference data for the comparison
            if (reference_metric_data is False or
                reference_metric_data is None or
                (isinstance(reference_metric_data, np.ndarray) and len(reference_metric_data) < 2)):
                print(''.join(['Features for ', metric, ' are missing. Please inspect the corresponding folder: a triplet of pickle files for at least one PDBID should be present.']))
            else:
                # Compute metrics performing paired measurements for the data of the reference protein
                # and each candidate protein in the dataset
                entries = [(reference_metric_data, comparisonInfo) for comparisonInfo in self.candidates[metric_index]]
                results = []
                if (self.debugging is False):
                    execution_handler = ExecutionHandler(self.cores, 12 * 60)
                    results = execution_handler.parallelize(self.compute_metrics, entries)
                else:
                    for candidate in entries:
                        result = self.compute_metrics(candidate)
                        if (result is not False):
                            results.append(result)
                final_results.extend(results)
                # Store the results
                data = pd.read_csv((io.StringIO('\n'.join(['\t'.join(result) for result in np.array(final_results)]))),
                                   names=['structureId', self._column_headers[metric_index]], sep='\t')
                data.sort_values(by=[self._column_headers[metric_index]], ascending=[True]).to_csv(output_result_path,
                                                                                                  index=False)

    def set_candidates(self):
        # Load the list of the candidate proteins that will be compared with the reference protein
        self.candidates = defaultdict()
        filelist_info_path = os.path.join(self.root_disk, self.features_file_list)
        data = pd.read_csv(filelist_info_path, sep='\t')
        if('int' in data['structureId'].dtype.name):
            data['structureId'] = data['structureId'].apply(str)
        for i in range(len(self._feature_file_suffixes)):
            self.candidates[i] = data[(data['metric'] == i)].to_records(index=False).tolist()
