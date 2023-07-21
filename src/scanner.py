from collections import defaultdict
from src.pdbhandler import PDBHandler
from src.bhattacharyyadistance import BhattacharyyaDistance
import numpy as np
import pandas as pd
import os
import time
import traceback
import pickle
import io
import subprocess
import math
import glob
from src.executionhandler import ExecutionHandler
from scipy.linalg import det
from src.protos import features_pb2

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
        self._feature_filelist_headers = ['structureId', 'chainId']
        self.cores = 1
        self.update_features = []
        self.recompute_features = False
        self.extend_feature_data = False
        self.pdb_validation = False
        self.store_format = '.proto' # '.proto', '.pkl'
        self.metric_indices = []
        self.execution_cutoff = 12 * 60

    # Load pre-computed features for a protein required for the calculation of a specified metric
    def load_features(self, pdb_id, chain_id, naming_extension=''):
        store_data_path = os.path.join(self.root_disk, self.features_path)
        naming_extension = ''.join(['_', naming_extension]) if naming_extension != '' else ''
        features_data_path = os.path.join(store_data_path, ''.join(
            [pdb_id.replace('_', '-'), '_', chain_id, naming_extension, self.store_format]))
        metric_data = False
        if (os.path.exists(features_data_path)):
            if (self.store_format == '.pkl'):
                if(os.path.exists(features_data_path)):
                    with open(features_data_path, 'rb') as handle:
                        metric_data = pickle.load(handle)
            else:
                stored_data = features_pb2.Features()
                with open(features_data_path, 'rb') as proto_file:
                    stored_data.ParseFromString(proto_file.read())
                metric_data = defaultdict()
                metric_data['mean'] = np.array(stored_data.phipsi.mean)
                if(np.array(stored_data.phipsi.cov).shape[0] != 4):
                    print(pdb_id, chain_id, naming_extension)
                    return False
                metric_data['cov'] = np.array(stored_data.phipsi.cov).reshape(2, 2)
                metric_data['det_cov'] = stored_data.phipsi.det_cov
                metric_data['distances'] = np.array(stored_data.distances)
                metric_data['triangles'] = stored_data.triangles
        else:
            # backwards compatibility
            for metric_index in range(len(self._column_headers)):
                if (metric_index not in self.metric_indices):
                    continue
                features_data_path = os.path.join(store_data_path, ''.join(
                    [pdb_id.replace('_', '-'), '_', chain_id, '_', naming_extension,
                     self._feature_file_suffixes[metric_index], '.pkl']))
                if (os.path.exists(features_data_path)):
                    with open(features_data_path, 'rb') as handle:
                        metric_data = pickle.load(handle)
                    feature_data = defaultdict()
                    if (metric_index == 0):
                        metric_data = np.array(metric_data['phiPsiPairs'])
                        angle_data = np.array(metric_data)
                        feature_data['mean'] = np.mean(angle_data, axis=0)
                        feature_data['cov'] = np.cov(angle_data.T)
                        feature_data['det_cov'] = det(feature_data['cov'])
                        if (feature_data['det_cov'] <= 0):
                            feature_data = False
                        metric_data = feature_data
                    if (metric_index == 1):
                        feature_data['distances'] = metric_data
                    if (metric_index == 2):
                        feature_data['triangles'] = metric_data
                    metric_data = feature_data
        return metric_data

    def wasserstein_distance(self, u_values, v_values):
        u_sorter = np.argsort(u_values)
        v_sorter = np.argsort(v_values)

        all_values = np.concatenate((u_values, v_values))
        all_values.sort(kind='mergesort')

        # Compute the differences between pairs of successive values of u and v.
        deltas = np.diff(all_values)

        # Get the respective positions of the values of u and v among the values of
        # both distributions.
        u_cdf_indices = u_values[u_sorter].searchsorted(all_values[:-1], 'right')
        v_cdf_indices = v_values[v_sorter].searchsorted(all_values[:-1], 'right')

        # Calculate the CDFs of u and v.
        u_cdf = u_cdf_indices / u_values.size
        v_cdf = v_cdf_indices / v_values.size

        return np.sum(np.multiply(np.abs(u_cdf - v_cdf), deltas))

    def compute_metrics(self, entry):
        reference_metric_data, comparison_info, metric_indices = entry
        # Load the required features
        metric_data = self.load_features(*comparison_info)
        output = False
        if (metric_data is not False):
            output = ['_'.join([comparison_info[0], comparison_info[1]])]
            if len(comparison_info) > 2:
                output.append(comparison_info[2])
            for metric_index in metric_indices:
                result = False
                # Compute the metric
                try:
                    if (metric_index == 0):
                        result = BhattacharyyaDistance.multivariate_compare_group(reference_metric_data, metric_data)
                    elif (metric_index == 1):
                        result = math.log10(self.wasserstein_distance(reference_metric_data['distances'], metric_data['distances']) + 1)
                    elif (metric_index == 2):
                        result = np.exp(abs(reference_metric_data['triangles'] - metric_data['triangles'])) - 1
                except:
                    print(''.join([output[0], '\n', traceback.format_exc()]))
                    result = False
                # Construct output
                if (result is not False):
                    output.append(result)
                else:
                    output = False
                    break
        return output

    def collect_features(self, pdb_filepaths=[]):
        # Check if filelist exists
        filelist_path = os.path.sep.join([self.root_disk, ''.join([self.features_path, '.csv'])])
        filelist_exists = os.path.exists(filelist_path)
        # Re-compute the features by user's choice
        if (self.recompute_features is True):
            self.update_features = self._feature_file_suffixes.copy()
        # Extract features from the PDB dataset that was set in the configuration
        if (filelist_exists is False or len(self.update_features) > 0 or self.extend_feature_data is True):
            os.makedirs(filelist_path.replace('.csv', ''), exist_ok=True)
            pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
            if len(pdb_filepaths) == 0:
                input_set = [filename for filename in os.listdir(pdb_dataset_path) if
                                os.path.isfile(os.path.join(pdb_dataset_path, filename))]
            else:
                input_set = pdb_filepaths
            if (self.debugging is False):
                # Extract features in parallel
                execution_handler = ExecutionHandler(self.cores, self.execution_cutoff)
                execution_handler.parallelize(self.extract_features, input_set)
            else:
                for input_data in input_set:
                    self.extract_features(input_data)
            if len(pdb_filepaths) == 0 or filelist_exists is False:
                self.create_feature_filelist()
            else:
                self.update_feature_filelist(pdb_filepaths)
        else:
            print('Features have already been computed. Please check if there is a directory misconfiguration if you were not expecting this outcome and execute again.')

    def create_feature_filelist(self, output_dir=''):
        # Create a file listing all filenames of the stored features for each PDB chain
        try:
            if output_dir == '':
                output_dir = os.path.sep.join([self.root_disk, self.features_path])
            output = subprocess.run(['timeout', '30s', 'ls', output_dir], capture_output=True)
            lines = output.stdout.decode("utf-8").split('\n')
            data_info = ['\t'.join(self._feature_filelist_headers)]
            for line in lines:
                if (len(line) > 2):
                    output_line = self.parse_feature_filename(line)
                    if len(output_line) > 0:
                        data_info.append('\t'.join(output_line))
            with open(''.join([output_dir, '.csv']), 'w') as datalist:
                datalist.write('\n'.join(data_info))
                datalist.write('\n')
        except:
            print(traceback.format_exc())

    def parse_feature_filename(self, filename):
        parts = filename.split('_')
        parsed = []
        if len(parts) > 1:
            parsed = ['_'.join(parts[:-1]), parts[-1].replace(self.store_format, '')]
        return parsed

    def update_feature_filelist(self, pdb_filepaths, folder_list_path=''):
        # Update existing file list with newly added structure ids
        if folder_list_path == '':
            folder_list_path = os.path.sep.join([self.root_disk, self.features_path])
        list_path = ''.join([folder_list_path, '.csv'])
        new_entries = []
        # Check for new proto files
        for file_path in pdb_filepaths:
            file_name = file_path.split(os.sep)
            if len(file_name) > 0:
                file_name = file_name[-1].replace('.pdb', '*').replace('_', '-')
                for data_path in glob.glob(''.join([folder_list_path, os.sep, file_name])):
                    if(data_path.endswith(self.store_format)):
                        data_filename = data_path.split(os.sep)
                        if len(data_filename) > 0:
                            data_filename = data_filename[-1]
                            output_line = self.parse_feature_filename(data_filename)
                            if len(output_line) > 1:
                                new_entries.append(output_line)
        # Append to existing file list
        data = pd.read_csv(list_path, sep='\t')
        data = pd.concat([data, pd.DataFrame(new_entries, columns=self._feature_filelist_headers)], ignore_index=True)
        data = data.drop_duplicates().reset_index(drop=True)
        data.sort_values(by=self._feature_filelist_headers, ascending=[True] * len(self._feature_filelist_headers))\
            .to_csv(list_path, index=False, sep='\t')

    def extract_feature(self, chain_id, output_path, full_pdbfile_path, feature, **kwargs):
        # Extract features for a PDB chain
        feature_data = False
        status = False
        try:
            pdbhandler = PDBHandler()
            pdbhandler.root_disk = self.root_disk
            pdbhandler.verbose = self.verbose
            if 'residue_selection' in kwargs:
                pdbhandler.residue_selection = kwargs['residue_selection']
            if (feature in self.update_features or len(self.update_features) == 0):
                if (feature == 'angles'):
                    pdbhandler.fetch_angles(full_pdbfile_path, chain_id)
                    feature_data = pdbhandler.chain_angles[chain_id]
                    status = len(feature_data['psi']) > 2 and len(feature_data['phi']) > 2
                    if(status):
                        angle_data = feature_data['phiPsiPairs']
                        if 'raw_data' not in kwargs:
                            feature_data = defaultdict()
                            angle_data = np.array(angle_data)
                            feature_data['mean'] = np.mean(angle_data, axis=0)
                            feature_data['cov'] = np.cov(angle_data.T)
                            assert feature_data['cov'].shape[0] == 2
                            feature_data['det_cov'] = det(feature_data['cov'])
                            if(feature_data['det_cov'] <= 0):
                                feature_data = False
                        else:
                            feature_data = angle_data
                elif (feature == 'distances'):
                    pdbhandler.calculate_residue_distances(full_pdbfile_path, chain_id)
                    feature_data = pdbhandler.residue_distances
                    status = len(pdbhandler.residue_distances) != 0
                elif (feature == 'triangles'):
                    pdbhandler.load_points(full_pdbfile_path, chain_id)
                    if (len(pdbhandler.points) > 2):
                        triangles = pdbhandler.get_mesh_triangles(pdbhandler.points)
                        if (triangles is not None):
                            feature_data = np.log(np.asarray(triangles).shape[0])
                    status = feature_data is not False
                else:
                    feature_data = False
        except:
            feature_data = False
            if (self.verbose is True):
                print(traceback.format_exc())
        if (feature_data is False):
            with open(os.path.join('logs', ''.join([time.strftime('%Y-%m-%d', time.localtime()), '_feat_extraction.log'])),
                      'a', encoding='utf-8') as fail_log_file:
                fail_log_file.write(
                    ' '.join(['feature_extraction', chain_id, output_path, full_pdbfile_path, feature, '\n']))
            if(self.verbose is True):
                print(' '.join(
                    [full_pdbfile_path.split('.')[-2].split(os.path.sep)[-1], chain_id, ': ', feature, ' were not extracted [',
                     output_path, ']']))
        return feature_data, status

    def store_features(self, feature_data, output_path):
        if(self.store_format == '.pkl'):
            if (os.path.exists(output_path)):
                with open(output_path, 'rb') as handle:
                    data_dict = pickle.load(handle)
            else:
                data_dict = {}
            for feature in feature_data:
                data_dict[feature] = feature_data[feature]
            with open(output_path, 'wb') as handle:
                pickle.dump(data_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
        else:
            if(math.isinf(feature_data['triangles']) is False and feature_data['angles'] is not False):
                features = features_pb2.Features()
                if (os.path.exists(output_path)):
                    with open(output_path, 'rb') as proto_file:
                        features.ParseFromString(proto_file.read())
                if('angles' in feature_data):
                    features.phipsi.mean.extend(list(feature_data['angles']['mean'].reshape(-1)))
                    features.phipsi.cov.extend(list(feature_data['angles']['cov'].reshape(-1)))
                    features.phipsi.det_cov = feature_data['angles']['det_cov']
                if('distances' in feature_data):
                    features.distances.extend(list(feature_data['distances']))
                if ('triangles' in feature_data):
                    features.triangles = feature_data['triangles']
                with open(output_path, 'wb') as proto_file:
                    proto_file.write(features.SerializeToString())

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
                feature_data = defaultdict()
                status = False
                for feature_index, feature in enumerate(self._feature_file_suffixes):
                    feature_output_path = ''.join([output_path, self.store_format])
                    # Skip if this feature is not set to be recomputed when recomputing is enabled
                    if (len(self.update_features) > 0 and feature not in self.update_features and self.extend_feature_data is False):
                        continue
                    # New extractions only
                    if(self.extend_feature_data is True and os.path.exists(feature_output_path) is True):
                        continue
                    feature_data[feature], status = self.extract_feature(chain_id, feature_output_path, full_pdbfile_path, feature)
                    if (status is False):
                        break
                if (status is not False):
                    self.store_features(feature_data, feature_output_path)

        return result

    def scan_candidates(self):
        # Compute metrics in parallel for each candidate protein in the dataset
        if (os.path.exists(self.metrics_output_path) is False):
            os.makedirs(self.metrics_output_path)
        self.metric_indices = []
        output_result_paths = defaultdict()
        for metric_index, metric in enumerate(self._feature_file_suffixes):
            output_result_paths[metric_index] = ''.join([self.metrics_output_path, os.path.sep, self.reference_pdb_id,
                                          '_', self.reference_chain_id, '_', self._column_headers[metric_index], '.csv'])
            if (os.path.exists(output_result_paths[metric_index]) is True):
                print(''.join([self._column_headers[metric_index], ' is already computed. If you wish to recompute it, delete',
                               ' the following file and execute the method again:\n', output_result_paths[metric_index]]))
                continue
            self.metric_indices.append(metric_index)
        if(len(self.metric_indices) > 0):
            # Load reference features
            reference_metric_data = self.load_features(self.reference_pdb_id, self.reference_chain_id)
            if reference_metric_data is False:
                raise Exception('Reference features are not available.')
            # Compute metrics performing paired measurements for the data of the reference protein
            # and each candidate protein in the dataset
            entries = [(reference_metric_data, comparisonInfo, self.metric_indices) for comparisonInfo in self.candidates]
            results = []
            if (self.debugging is False):
                execution_handler = ExecutionHandler(self.cores , self.execution_cutoff)
                results = execution_handler.parallelize(self.compute_metrics, entries)
            else:
                for candidate in entries:
                    result = self.compute_metrics(candidate)
                    if (result is not False):
                        results.append(result)
            results = np.array(results)
            # Output each metric in a csv file
            for metric_index in self.metric_indices:
                # Store the results
                data = pd.read_csv((io.StringIO('\n'.join(['\t'.join(result) for result in results[:, [0]+[self.metric_indices.index(metric_index)+1]]]))),
                                   names=['structureId', self._column_headers[metric_index]], sep='\t')
                data.sort_values(by=[self._column_headers[metric_index]],
                                 ascending=[True]).to_csv(output_result_paths[metric_index], index=False)

    def set_candidates(self):
        # Load the list of the candidate proteins that will be compared with the reference protein
        self.candidates = defaultdict()
        filelist_info_path = os.path.join(self.root_disk, self.features_file_list)
        data = pd.read_csv(filelist_info_path, sep='\t')
        if(len(data.columns) > len(self._feature_filelist_headers)): #backwards compatibility for pickle
            data = data.iloc[:,:len(self._feature_filelist_headers)]
        data = data.drop_duplicates().reset_index(drop=True)
        if('int' in data['structureId'].dtype.name):
            data['structureId'] = data['structureId'].apply(str)
        self.candidates = data.to_records(index=False).tolist()
