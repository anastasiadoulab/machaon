from src.pdbhandler import PDBHandler
import pandas as pd
import os
import io
from src.scanner import Scanner
from src.executionhandler import ExecutionHandler
from collections import defaultdict
import numpy as np

# Constrained mode of the comparison method (domain comparisons)
class DomainScanner(Scanner):

    def __init__(self):
        super().__init__()
        self._feature_filelist_headers = ['structureId', 'chainId', 'domain']

    def extract_features(self, entry):
        result = True
        store_data_path = os.path.sep.join([self.root_disk, self.features_path])
        pdb_path = entry
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        full_pdbfile_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path, pdb_path])
        if (self.pdb_validation is False):
            chains = pdbhandler.get_protein_chain_ids(full_pdbfile_path)
        else:
            chains = pdbhandler.fetch_pdb_peptidic_chains(full_pdbfile_path)
        if chains is False or len(chains) == 0:
            result = False
        else:
            # For each PDB chain of the PDB file
            for chainID in chains:
                if (len(chainID.strip()) == 0):
                    print(pdb_path, 'contains an empty chain identifier. This chain will be ignored.')
                    continue
                available_residues = pdbhandler.get_residue_range(full_pdbfile_path, chainID)
                residue_range = available_residues['fullRange']
                pdbhandler.verbose = self.debugging
                pdbhandler.get_uniprot_accession_number(chainID, pdb_path.split('.')[0], full_pdbfile_path)
                if (pdbhandler.uniprot_accession_number != ''):
                    pdbhandler.get_domain_information()
                    # If there is domain information available, determine the residue positions
                    # and extract the required features
                    for domain in pdbhandler.domains:
                        name, start, end = domain
                        if (residue_range[0] > int(end) or residue_range[-1] < int(start)):
                            continue
                        name = pdbhandler.sanitize_domain_name(name)
                        residue_selection = list(range(int(start), int(end) + 1))
                        output_path = os.path.sep.join([store_data_path, pdb_path.replace('_', '-').replace('.pdb', ''.join(['_', chainID]))])
                        feature_data = defaultdict()
                        status = False
                        for feature_index, feature in enumerate(self._feature_file_suffixes):
                            feature_output_path = ''.join([output_path, '_', name, self.store_format])
                            # Skip if this feature is not set to be recomputed when recomputing is enabled
                            if (len(self.update_features) > 0 and feature not in self.update_features):
                                continue
                            # New extractions only
                            if (os.path.exists(feature_output_path) is True):
                                os.remove(feature_output_path)
                            if (self.extend_feature_data is True and os.path.exists(feature_output_path) is True):
                                continue
                            feature_data[feature], status = self.extract_feature(chainID, feature_output_path, full_pdbfile_path, feature, residue_selection=residue_selection)
                            if(status is False):
                                break
                        if(status is not False):
                            self.store_features(feature_data, feature_output_path)
        return result

    def scan_candidates(self):
        if (os.path.exists(self.metrics_output_path) is False):
            os.makedirs(self.metrics_output_path)
        pdbhandler = PDBHandler()
        pdbhandler.root_disk = self.root_disk
        pdbhandler.structure_id = self.reference_pdb_id
        pdbhandler.verbose = self.debugging
        pdb_dataset_path = os.path.sep.join([self.root_disk, self.pdb_dataset_path])
        full_pdbfile_path = ''.join([pdb_dataset_path, os.path.sep, self.reference_pdb_id, '.pdb'])
        pdbhandler.get_uniprot_accession_number(self.reference_chain_id, self.reference_pdb_id, full_pdbfile_path)
        # Get domain information for each candidate and load its corresponding features
        if (pdbhandler.uniprot_accession_number != ''):
            # For each reference domain
            pdbhandler.get_domain_information()
            for domainInfo in pdbhandler.domains:
                name, start, end = domainInfo
                domain = pdbhandler.sanitize_domain_name(name)
                print("Processing for: ", domain)
                reference_metric_data = self.load_features(self.reference_pdb_id, self.reference_chain_id,
                                                           naming_extension=domain)
                if (reference_metric_data is False):
                    print(''.join(
                        ['Failed to extract data from ', self.reference_pdb_id, ' - ', self.reference_chain_id, ' - ', domain]))
                    continue
                self.metric_indices = []
                output_result_paths = defaultdict()
                for metric_index, metric in enumerate(self._feature_file_suffixes):
                    output_result_paths[metric_index] = ''.join([self.metrics_output_path, os.path.sep, self.reference_pdb_id,
                                                  '_', self.reference_chain_id, '_', domain, '_', self._column_headers[metric_index], '.csv'])
                    if (os.path.exists(output_result_paths[metric_index]) is True):
                        print(''.join([self._column_headers[metric_index], ' is already computed. If you wish to recompute it, delete',
                                       ' the following file and execute the method again:\n', output_result_paths[metric_index]]))
                        continue
                    self.metric_indices.append(metric_index)
                # Compute the metrics for every domain in candidate set
                entries = [(reference_metric_data, comparisonInfo, self.metric_indices) for comparisonInfo in self.candidates]
                results = []
                if (self.debugging is False):
                    execution_handler = ExecutionHandler(self.cores, 12 * 60)
                    results = execution_handler.parallelize(self.compute_metrics, entries)
                else:
                    for candidate in entries:
                        result = self.compute_metrics(candidate)
                        if(result is not False):
                            results.append(result)
                # Output each metric in a csv file
                if len(results) > 0:
                    results = np.array(results)
                    for metric_index in self.metric_indices:
                        data = pd.read_csv(
                            io.StringIO('\n'.join(['\t'.join(result) for result in results[:, [0, 1]+[metric_index+2]]])),
                            names=['structureId', 'domain', self._column_headers[metric_index]], sep='\t')
                        data.sort_values(by=[self._column_headers[metric_index]], ascending=[True]).to_csv(output_result_paths[metric_index], index=False)

    # Handle different types of feature data filenames (e.g. whole structure vs domain based)
    def parse_feature_filename(self, filename):
        parts = filename.split('_')
        parsed = []
        if len(parts) > 2:
            parsed = ['_'.join(parts[:-2]), parts[-2], parts[-1].replace(self.store_format, '')]
        return parsed
