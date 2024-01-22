import logging
import sys
import threading

sys.path.append('..')
import pandas as pd
import os
from src.grpc_core import jobreceiver_pb2_grpc
from src.grpc_core import jobreceiver_pb2
from src.configurationmanager import ConfigurationManager
from src.machaon import Machaon
from collections import defaultdict
import copy
import psutil
import numpy as np
import time
import shutil
import requests
import glob
import hashlib
import traceback
import uuid
import zipfile
import re


class JobReceiver(jobreceiver_pb2_grpc.JobReceiverServicer):

    def __init__(self):
        self.root_disk = ''
        self.output_path = ''
        self.temp_folder = 'temp'
        self.failed_folder = 'failed'
        self.available_cores = 8
        self.transfer_chunk_size = 1024
        self.pdb_folders = ['PDBs_vir', 'PDBs', 'AF4']
        self.cache_folder = 'PDBs_new'
        self.current_jobid = -1
        self.rcsb_url = 'https://files.rcsb.org/view/'
        self.af_url = 'https://alphafold.ebi.ac.uk/files/'
        self.esm_url = 'https://api.esmatlas.com/fetchPredictedStructure/'
        self.list_data_dirs = {'Pub_Viral_PDB_Candidate_Set': 'PDBs_vir',
                               'Human_PDB_Candidate_Set_v1': 'PDBs',
                               'Human_AF4_Candidate_Set': 'AF4'}
        self.comparison_modes = ['whole', 'domain', 'segment']
        self.alignment_levels = ['primary', 'secondary', 'hydrophobicity', 'mixed']
        self.download_delay = 4
        self.current_thread = None
        self.logger = logging.getLogger("grpc_logger")
        logging_handler = logging.StreamHandler()
        self.logger.addHandler(logging_handler)
        os.makedirs(os.path.join(os.getcwd(), 'logs'), exist_ok=True)
        logfile_handler = logging.FileHandler(os.path.join(os.getcwd(), 'logs', 'grpc.log'))
        self.logger.addHandler(logfile_handler)
        self.logger.setLevel(logging.DEBUG)
        self.measurement_delay = 10
        self.spawn_time = None
        self.maximum_process_time = 3600 * 2
        self.custom_reference = None

    def setup_directories(self, root_disk, output_path, cores):
        self.root_disk = root_disk
        self.output_path = output_path
        self.available_cores = cores
        # Create folders for custom requests
        os.makedirs(os.path.join(self.root_disk, self.cache_folder), exist_ok=True)
        os.makedirs(os.path.join(self.root_disk, self.temp_folder), exist_ok=True)
        os.makedirs(os.path.join(self.root_disk, self.failed_folder), exist_ok=True)
        for mode in self.comparison_modes:
            if mode == 'segment':
                continue
            os.makedirs(os.path.join(self.root_disk, '_'.join(['DATA', self.cache_folder, mode])), exist_ok=True)

    def setup_custom_data(self, structure_id_list, custom_folder_name, comparison_mode):
        # Create required folders
        data_path = os.path.join(self.root_disk, custom_folder_name)
        self.remove_data_dir(data_path)
        os.makedirs(data_path)
        feature_data_path = os.path.join(self.root_disk, '_'.join(['DATA', custom_folder_name, comparison_mode]))
        self.remove_data_dir(feature_data_path)
        os.makedirs(feature_data_path)
        fetched_structures, filenames = self.fetch_custom_data(structure_id_list, custom_folder_name, comparison_mode)
        return fetched_structures, filenames

    def fetch_custom_data(self, structure_id_list, custom_folder_name, comparison_mode):
        # Check existing folders if they are present soft copy to new
        # folder & data_folder else download
        fetched_structures = []
        uncached, filenames = self.populate_request_folders(structure_id_list, custom_folder_name, comparison_mode)
        if len(uncached) > 0:
            fetched_structures = self.download_structures(uncached)
            for structure_id in fetched_structures:
                destination_path = os.path.join(self.root_disk, custom_folder_name, ''.join([structure_id, '.pdb']))
                if os.path.exists(destination_path) is False:
                    os.symlink(os.path.join(self.root_disk, self.cache_folder, ''.join([structure_id, '.pdb'])),
                               destination_path)
        return fetched_structures, filenames

    def download_structures(self, structure_id_list):
        fetched_structures = []
        cache_folder_path = os.path.join(self.root_disk, self.cache_folder)
        for structure_id in structure_id_list:
            if os.path.exists(os.path.join(cache_folder_path, ''.join([structure_id, '.pdb']))) is True:
                fetched_structures.append(structure_id)
                continue
            try:
                print("G Do")
                download_url = ""
                if len(structure_id) == 4:
                    download_url = "".join([self.rcsb_url, structure_id, '.pdb'])
                elif structure_id.startswith("AF-"):
                    download_url = "".join([self.af_url, structure_id, '.pdb'])
                elif structure_id.startswith("MGYP"):
                    matched = re.findall("MGYP[0-9]{12}", structure_id)
                    if len(matched) > 0:
                        download_url = "".join([self.esm_url, matched[0], '.pdb'])
                if download_url != "":
                    print(download_url)
                    result = requests.get(download_url)
                    if result.status_code == 200:
                        with open(os.path.join(self.root_disk, self.cache_folder, ''.join([structure_id, '.pdb'])), 'wb') as structure_file:
                            structure_file.write(result.content)
                        fetched_structures.append(structure_id)
            except:
                logging.debug(' '.join(['jobid:', repr(self.current_jobid), ' | Fetching failed: ', structure_id, traceback.format_exc()]))
            time.sleep(self.download_delay)
        return fetched_structures

    def populate_request_folders(self, structure_id_list, request_folder, comparison_mode):
        remaining = [pdb_id.replace('_', '-') for pdb_id in structure_id_list]
        data_directories = self.pdb_folders.copy()
        data_directories.append(self.cache_folder)
        features_request_folder = '_'.join(['DATA', request_folder, comparison_mode])
        all_filenames = []
        for folder_name in data_directories:
            features_folder_name = '_'.join(['DATA', folder_name, comparison_mode])
            feature_list_path = os.path.join(self.root_disk, ''.join([features_folder_name, '.csv']))
            if os.path.exists(feature_list_path) is False:
                continue
            data = pd.read_csv(feature_list_path, sep='\t')
            data = data[data['structureId'].isin(remaining)]
            unique_ids = list(set(data['structureId']))
            structure_ids = list(data['structureId'])
            chain_ids = list(data['chainId'])
            domains = []
            if comparison_mode == 'domain':
                domains = list(data['domain'])
            for cached_id in unique_ids:
                cached_structure = cached_id.replace('-model-v', '-model_v') # handle AlphaFold filenames
                destination_path = os.path.join(self.root_disk, request_folder, ''.join([cached_structure, '.pdb']))
                if os.path.exists(destination_path) is False:
                    os.symlink(os.path.join(self.root_disk, folder_name, ''.join([cached_structure, '.pdb'])),
                               destination_path)
            if comparison_mode != 'segment':
                domain_mode = comparison_mode == 'domain'
                for id_index, structure_id in enumerate(structure_ids):
                    filename = [structure_id, chain_ids[id_index]]
                    if domain_mode is True:
                        filename.append(domains[id_index])
                    full_filename = '_'.join(filename)
                    destination_path = os.path.join(self.root_disk, features_request_folder, ''.join([full_filename, '.proto']))
                    if os.path.exists(destination_path) is False:
                        os.symlink(os.path.join(self.root_disk, features_folder_name, ''.join([full_filename, '.proto'])),
                                   destination_path)
                    if filename not in all_filenames:
                        all_filenames.append(filename)
            remaining = [structure_id for structure_id in remaining if structure_id not in unique_ids]
            if len(remaining) == 0:
                break
        return remaining, all_filenames

    def execute_machaon(self, entry):
        config_manager = ConfigurationManager()
        configuration = copy.deepcopy(config_manager._template_config)
        configuration['rootDisk'] = self.root_disk
        configuration['pdbDatasetPath'] = entry['data_path']
        configuration['isReferenceViral'] = True
        configuration['viralContentExists'] = True
        configuration['referencePDBID'] = entry['pdb_id']
        configuration['referenceChainID'] = entry['chain_id']
        configuration['noThirdPartyData'] = not entry['meta']
        configuration['GOSearch'] = [entry['go']]
        configuration['outputPath'] = entry['output_path']
        configuration['alignmentLevel'] = entry['alignment_level']
        configuration['alignmentBackend'] = 'parasail'
        configuration['comparisonMode'] = entry['comparison_mode']
        configuration['excludedPDBIDs'] = entry['excluded']
        if entry['comparison_mode'] == 'segment' and entry['segment_start'] > 0:
            configuration['segments'] = [[i for i in range(entry['segment_start'], entry['segment_end'] + 1)]]
        machaon_core = Machaon()
        #machaon_core.debugging = True
        machaon_core.max_cores = entry['cores']
        try:
            if(len(entry['new_pdbs']) > 0 and entry['comparison_mode'] != 'segment'):
                file_paths = ['.'.join([structure_id, 'pdb']) for structure_id in entry['new_pdbs']]
                machaon_core.process_new_pdbs([configuration], file_paths)
            machaon_core.perform_comparisons([configuration])
        except:
            self.logger.debug(traceback.format_exc())
        self.prepare_result(entry['request_hash'], entry['comparison_mode'], entry['fetched'])
        self.current_thread = None

    def clean_up(self, request_hash):
        output_path = os.path.join(self.output_path, request_hash)
        data_path = os.path.join(self.root_disk, request_hash)

        # Kill computing thread
        if self.current_thread is not None:
            self.current_thread.kill()
            self.current_thread = None

        # Reset thread spawn time
        self.spawn_time = None

        # Remove data for custom reference (for preset lists)
        print("custom ref: ", self.custom_reference)
        if self.custom_reference is not None:
            self.remove_custom_reference(request_hash)

        # Remove data folders
        for mode in self.comparison_modes:
            if mode == 'segment':
                continue
            features_dir = os.sep.join([self.root_disk, '_'.join(['DATA', request_hash, mode])])
            self.remove_data_dir(features_dir)
            if os.path.exists(''.join([features_dir, '.csv'])):
                os.remove(''.join([features_dir, '.csv']))

        self.remove_data_dir(data_path)

        # Move partial output to the 'failed' folder for inspection
        if os.path.exists(output_path):
            destination_path = os.path.join(self.root_disk, self.failed_folder, request_hash)
            if os.path.exists(destination_path) is True:
                self.remove_data_dir(destination_path)
            shutil.copytree(output_path, destination_path)
            self.remove_data_dir(output_path)

        # Handle previous remnants (e.g. server crash)
        compressed_file_path = os.path.join(self.output_path, '.'.join([request_hash, 'zip']))
        if os.path.exists(compressed_file_path):
            os.remove(compressed_file_path)

    # Clean up for custom reference when preset lists are selected
    def remove_custom_reference(self, request_hash):
        for mode in self.comparison_modes:
            features_dir = os.sep.join([self.root_disk, '_'.join(['DATA', request_hash, mode])])
            if os.path.exists(features_dir):
                for data_filepath in glob.glob(''.join([features_dir, os.sep, self.custom_reference.replace('_', '-'), '*'])):
                    if (data_filepath.endswith('.proto')):
                        os.remove(data_filepath)
        structures_path = os.path.join(self.root_disk, request_hash)
        if os.path.exists(structures_path):
            reference_structure_path = os.path.join(structures_path, '.'.join([self.custom_reference, 'pdb']))
            if os.path.exists(reference_structure_path):
                os.remove(reference_structure_path)
        self.custom_reference = None

    @staticmethod
    def remove_data_dir(data_path):
        if os.path.exists(data_path):
            if os.path.islink(data_path):
                os.remove(data_path)
            else:
                shutil.rmtree(data_path)

    @staticmethod
    def update_local_cachelist(data_folder_path, new_files, mode):
        machaon = Machaon()
        machaon.manage_file_list(data_folder_path, new_files, mode)

    def prepare_result(self, request_hash, comparison_mode, new_structures):
        features_path = os.path.join(self.root_disk, '_'.join(['DATA', request_hash, comparison_mode]))
        structures_path = os.path.join(self.root_disk, request_hash)
        output_path = os.path.join(self.output_path, request_hash)

        # Create a temp folder that will include the full
        # output folder and new pdb & feature files
        identifier = str(uuid.uuid4())
        temp_folder_path = os.path.join(self.root_disk, self.temp_folder, identifier)
        os.makedirs(temp_folder_path)
        new_pdbs_folder = os.path.join(temp_folder_path, 'pdbs')
        new_data_folder = os.path.join(temp_folder_path, comparison_mode)

        # Remove data for custom reference (for preset lists)
        print("custom ref: ", self.custom_reference)
        if self.custom_reference is not None:
            self.remove_custom_reference(request_hash)

        # Compress full output folder and move it into temp folder
        if os.path.exists(output_path):
            self.compress_dir(output_path, os.path.join(temp_folder_path, '.'.join([request_hash, 'zip'])))

        # Handle uncached structures
        print("Uncached structures: ", repr(new_structures))
        if len(new_structures) > 0:
            os.makedirs(new_pdbs_folder)
            os.makedirs(new_data_folder)

            # Copy new pdbs to temp folder's new pdbs folder
            for structure_id in new_structures:
                file_path = os.path.join(structures_path, '.'.join([structure_id, 'pdb']))
                if os.path.exists(file_path):
                    shutil.copy(file_path, new_pdbs_folder)

            # Move features from new pdbs to new features folder
            if os.path.exists(features_path):
                custom_data_path = os.path.join(self.root_disk, '_'.join(['DATA', self.cache_folder, comparison_mode]))
                new_files = []
                for structure_id in new_structures:
                    for data_filepath in glob.glob(''.join([features_path, os.sep, structure_id.replace('_', '-'), '*'])):
                        if (data_filepath.endswith('.proto')):
                            data_filename = os.path.basename(data_filepath)
                            destination_path = os.path.join(self.root_disk, custom_data_path, data_filename)
                            if os.path.exists(destination_path) is False:
                                shutil.copyfile(data_filepath, destination_path)
                                if data_filename not in new_files:
                                    new_files.append(data_filename[-1])
                                shutil.copyfile(data_filepath,  os.path.join(new_data_folder, data_filename))
                self.update_local_cachelist(custom_data_path, new_files, comparison_mode)

        # Delete request features file list
        list_path = os.path.join(self.root_disk, ''.join(['DATA', '_', request_hash, '_', comparison_mode, '.csv']))
        if os.path.exists(list_path):
            os.remove(list_path)

        # Compress the temp folder
        self.compress_dir(temp_folder_path, os.path.join(self.output_path, '.'.join([request_hash, 'zip'])))

        # Delete folders
        self.remove_data_dir(temp_folder_path)
        self.remove_data_dir(features_path)
        self.remove_data_dir(structures_path)

    def handle_preset_request(self, reference_structure, request_hash, list_data_folder, comparison_mode):
        features_folder_name = '_'.join(['DATA', list_data_folder,
                                         comparison_mode if comparison_mode != 'segment' else 'whole'])
        feature_list_path = os.path.join(self.root_disk, ''.join([features_folder_name, '.csv']))
        data = pd.read_csv(feature_list_path, sep='\t')
        data = data[data['structureId'] == reference_structure.replace('-model-v', '-model_v')]
        data_folder = request_hash
        no_reference = False
        downloaded_structures = False

        # Create virtual folders
        destination_path = os.path.join(self.root_disk, request_hash)
        if os.path.exists(destination_path) is False:
            os.symlink(os.path.join(self.root_disk, list_data_folder),
                       destination_path)
        if comparison_mode != 'segment':
            custom_features_folder = '_'.join(['DATA', request_hash, comparison_mode])
            custom_feature_list = os.path.join(self.root_disk, ''.join([custom_features_folder, '.csv']))
            destination_path = os.path.join(self.root_disk, custom_features_folder)
            if os.path.exists(destination_path) is False:
                os.symlink(os.path.join(self.root_disk, features_folder_name),
                            destination_path)
            destination_path = os.path.join(self.root_disk, custom_feature_list)
            if os.path.exists(destination_path) is False:
                shutil.copyfile(os.path.join(self.root_disk, feature_list_path),
                           os.path.join(self.root_disk, custom_feature_list))
        if (data.empty):
            # Fetch missing structure
            fetched_structures, filenames = self.fetch_custom_data([reference_structure], data_folder, comparison_mode)
            if len(fetched_structures) == 0 and len(filenames) == 0:
                no_reference = True
            downloaded_structures = len(fetched_structures) != 0

        return data_folder, data.empty, downloaded_structures, no_reference

    def create_custom_filelist(self, structure_id_list, folder_name, comparison_mode):
        feature_filename = '_'.join(['DATA', folder_name, comparison_mode])
        headers = ['structureId', 'chainId']
        if comparison_mode == 'domain':
            headers.append('domain')
        file_list = pd.DataFrame(structure_id_list, columns=headers)
        file_list.to_csv(os.path.join(self.root_disk, ''.join([feature_filename, '.csv'])), sep='\t', index=False)

    @staticmethod
    def compress_dir(source_path, output_path):
        with zipfile.ZipFile(output_path, "w", compression=zipfile.ZIP_DEFLATED) as zip_handle:
            # ziph is zipfile handle
            for root, dirs, files in os.walk(source_path):
                for file in files:
                    zip_handle.write(os.path.join(root, file),
                               os.path.relpath(os.path.join(root, file),
                                               os.path.join(source_path, '..')))

    @staticmethod
    def compute_file_hash(file_path):
        secure_hashing = hashlib.sha256()
        secure_hash = ''
        with open(file_path, 'rb') as file_handler:
            # Read and update hash string value in blocks of 4K
            for byte_block in iter(lambda: file_handler.read(4096), b""):
                secure_hashing.update(byte_block)
            secure_hash = secure_hashing.hexdigest()
        return secure_hash

    @staticmethod
    def denormalize_id(structure_id):
        # AlphaFold full names
        return structure_id.replace('-MODEL_V', '-model_v')

    def inspect_result(self, request_hash):
        valid = True
        output_folder = os.path.join(self.output_path, request_hash)
        metrics_folder = os.path.join(output_folder, 'metrics')
        candidates_folder = os.path.join(output_folder, 'candidates')
        if os.path.exists(metrics_folder):
            metric_filenames = glob.glob(''.join([metrics_folder, os.sep, '*.csv']))
            if len(metric_filenames) == 0:
                print("G M")
                valid = False
            else:
                for metric_filename in metric_filenames:
                    data = pd.read_csv(os.path.join(metrics_folder, metric_filename), sep=',')
                    if len(data.index) < 2:
                        print("G O")
                        valid = False
                        break
        else:
            print("G X")
            valid = False
        if valid is True:
            if os.path.exists(candidates_folder):
                filenames = glob.glob(''.join([candidates_folder, os.sep, '*.csv']))
                if len(filenames) == 0:
                    print("G Y")
                    valid = False
                else:
                    for filename in filenames:
                        if filename.endswith('-data.csv'):
                            continue
                        data = pd.read_csv(os.path.join(metrics_folder, filename), sep='\t')
                        if len(data.index) < 1:
                            print("G W")
                            valid = False
                            break
                if valid is True:
                    candidate_filenames = glob.glob(''.join([candidates_folder, os.sep, '*.html']))
                    if len(candidate_filenames) == 0:
                        print("G N")
                        valid = False
            else:
                print("G P")
                valid = False
        return valid


    def StartJob(self, request, context):
        response = 0
        request_hash = ''
        request_options = defaultdict()
        id_parts = self.denormalize_id(request.reference_id).split('_')
        if(self.current_thread is not None):
            if (time.time() - self.spawn_time > self.maximum_process_time):
                self.current_thread.kill()
                self.current_thread = None
                self.spawn_time = None
            else:
                response = 1
        if (np.mean(psutil.cpu_percent(interval=1, percpu=True)) > 30):
            response = 2
        if (len(id_parts) == 1 or len(id_parts) > 3):
            response = 3
        try:
            request_hash = repr(int(request.hash))
            request_id = request.request_id
            comparison_mode = self.comparison_modes[request.comparison_mode]
            alignment_level = self.alignment_levels[request.alignment_level]
            structure_ids = [self.denormalize_id(chosen_id.strip()) for chosen_id in list(request.structure_ids) if len(chosen_id.strip()) > 0]
        except:
            logging.debug(traceback.format_exc())
            response = 4
        # SEGMENT / DOMAIN NO APPLICABLE
        if response == 0:
            try:
                request_options['pdb_id'] = '_'.join(id_parts[:-1])
                request_options['chain_id'] = id_parts[-1]
                is_custom_list = len(structure_ids) > 0
                request_options['new_pdbs'] = []
                request_options['fetched'] = []
                request_options['excluded'] = []
                if(is_custom_list):
                    fetched_structures, composite_ids = self.setup_custom_data(structure_ids + [request_options['pdb_id']], request_hash, comparison_mode)
                    data_path = request_hash
                    if len(fetched_structures) == 0 and len(composite_ids) == 0:
                        response = 8
                    elif request_options['pdb_id'] not in fetched_structures:
                        reference_found = False
                        for composite_id in composite_ids:
                            if request_options['pdb_id'].replace('-model_v', '-model-v') in composite_id:
                                reference_found = True
                                break
                        if reference_found is False:
                            response = 9
                    if response == 0:
                        if len(composite_ids) > 0:
                            self.create_custom_filelist(composite_ids, request_hash, comparison_mode)
                        if len(fetched_structures) > 0:
                            request_options['new_pdbs'].extend(fetched_structures)
                        if request_options['pdb_id'] not in structure_ids:
                            request_options['excluded'].append(request_options['pdb_id'])
                        request_options['fetched'].extend(fetched_structures)
                else:
                    if request.listname in self.list_data_dirs:
                        data_path, custom_reference, fetched_reference, no_reference = self.handle_preset_request(request_options['pdb_id'],
                                                                                 request_hash,
                                                                                 self.list_data_dirs[request.listname],
                                                                                 comparison_mode)
                        if no_reference is True:
                            response = 7
                        else:
                            request_options['new_pdbs'].append(request_options['pdb_id'])

                        if custom_reference is True:
                            self.custom_reference = request_options['pdb_id']
                        else:
                            self.custom_reference = None

                        if fetched_reference is True:
                            request_options['fetched'].append(request_options['pdb_id'])

                    else:
                        response = 5
            except:
                self.logger.debug(traceback.format_exc())
                response = 6
        if response == 0:
            try:
                request_options['data_path'] = data_path
                request_options['root_disk'] = self.root_disk
                request_options['custom_list'] = is_custom_list
                request_options['cores'] = self.available_cores
                request_options['comparison_mode'] = comparison_mode
                request_options['segment_start'] = request.segment_start
                request_options['segment_end'] = request.segment_end
                request_options['meta'] = request.meta_analysis
                request_options['alignment_level'] = alignment_level
                request_options['go'] = request.go_term
                request_options['request_hash'] = request_hash
                request_options['request_id'] = request_id
                request_options['output_path'] = os.path.join(self.output_path, request_hash)
                os.makedirs(request_options['output_path'], exist_ok=True)

                # Handle previous remnants (e.g. server crash)
                compressed_file_path = os.path.join(self.output_path, '.'.join([request_hash, 'zip']))
                if os.path.exists(compressed_file_path):
                    os.remove(compressed_file_path)

                # Start job
                self.spawn_time = time.time()
                self.current_thread = threading.Thread(target=self.execute_machaon, args=([request_options]))
                self.current_thread.start()
                # self.execute_machaon(request_options)
            except:
                self.logger.debug(traceback.format_exc())
                response = 7
        return jobreceiver_pb2.JobStatus(request_id=request_id, status_code=response)

    def GetStatus(self, request, context):
        current_status = 1
        cpu_percentages = [percentage for percentage in psutil.cpu_percent(interval=1, percpu=True) if percentage > 0.0]
        if len(cpu_percentages) > 0:
            if np.mean(cpu_percentages) > 10.0:
                current_status = 2
        if self.current_thread is not None:
                current_status = 2
        return jobreceiver_pb2.ServerStatus(status_code=current_status)

    def DownloadResult(self, request, context):
        compressed_file_path = os.path.join(self.output_path, '.'.join([request.hash, 'zip']))
        status = 0
        if os.path.exists(compressed_file_path) is True:
            if self.current_thread is not None:
                self.current_thread.join()
                self.current_thread = None
            if self.inspect_result(request.hash) is True:
                print("G J")
                secure_hash = self.compute_file_hash(compressed_file_path)
                details = jobreceiver_pb2.JobDetails(request_id=request.request_id, hash=request.hash, secure_hash=secure_hash, status_code=status)
                yield jobreceiver_pb2.JobResult(file_info=details)
                with open(compressed_file_path, mode="rb") as archive:
                    while True:
                        chunk = archive.read(self.transfer_chunk_size)
                        if chunk:
                            entry_response = jobreceiver_pb2.JobResult(chunk_data=chunk)
                            yield entry_response
                        else:
                            os.remove(compressed_file_path)
                            return
            else:
                print("G R")
                self.clean_up(request.hash)
                status = -3
        else:
            try:
                status = -1
                # Detect a failed job
                if (np.mean(psutil.cpu_percent(interval=1, percpu=True)) < 10):
                    if((np.max(psutil.getloadavg()) / self.available_cores) < 0.1):
                        if self.spawn_time is not None:
                            if(time.time() - self.spawn_time > self.maximum_process_time):
                                print("G Q")
                                self.clean_up(request.hash)
                                status = -2
                        else:
                            print("G L")
                            status = -2
            except:
                print("G Z")
                self.logger.debug(traceback.format_exc())
        if status != 0:
            details = jobreceiver_pb2.JobDetails(request_id=request.request_id, hash=request.hash, secure_hash='',
                                                 status_code=status)
            yield jobreceiver_pb2.JobResult(file_info=details)

    def Synchronize(self, request, context):
        data = bytearray()
        status = 0
        identifier = str(uuid.uuid4())
        temp_folder_path = os.path.join(self.root_disk, self.temp_folder, identifier)
        if os.path.exists(temp_folder_path):
            shutil.rmtree(temp_folder_path)
        os.makedirs(temp_folder_path)
        file_path = ''.join([temp_folder_path, os.path.sep, identifier, '.zip'])
        secure_hash = None
        try:
            for request_chunk in request:
                if request_chunk.secure_hash:
                    secure_hash = request_chunk.secure_hash
                    continue
                data.extend(request_chunk.chunk_data)
            # if secure_hash is not None:
            #     raise Exception("no hash")
            with open(file_path, 'wb') as compressed_data:
                compressed_data.write(data)
            file_hash = self.compute_file_hash(file_path)
            if file_hash != secure_hash:
                raise Exception("hash mismatch")
            shutil.unpack_archive(file_path, temp_folder_path)
            for pdb_file in glob.glob(os.sep.join([temp_folder_path, 'PDBs_new', '*.pdb'])):
                shutil.copy(os.path.join(temp_folder_path, pdb_file), os.path.join(self.root_disk, self.cache_folder))
            for mode in ['whole', 'domain']:
                data_folder_path = os.path.join(self.root_disk, '_'.join(['DATA', self.cache_folder, mode]))
                new_files = []
                for data_file in glob.glob(os.sep.join([temp_folder_path, '_'.join(['DATA', self.cache_folder, mode]), '*.proto'])):
                    shutil.copy(data_file, data_folder_path)
                    if data_file not in new_files:
                        new_files.append(os.path.basename(data_file))
                self.update_local_cachelist(data_folder_path, new_files, mode)
            self.remove_data_dir(temp_folder_path)
            if os.path.exists(file_path):
                os.remove(file_path)
        except:
            self.logger.debug(traceback.format_exc())
            status = 1

        return jobreceiver_pb2.ServerStatus(status_code=status)

if __name__ == '__main__':
    reference_id = '6VXX_A'
    pdb_id = '4Z4D'
    job_receiver = JobReceiver()
    job_receiver.setup_directories('/media/panos/69aa1d93-46c8-45e1-a1e6-45cacfe63f15/structural-data',
                                   '/media/panos/69aa1d93-46c8-45e1-a1e6-45cacfe63f15/structural-data/output', 8)
    # job_receiver.handle_preset_request(pdb_id, '0000000000', job_receiver.list_data_dirs['Human_PDB_Candidate_Set_v1'], 'whole')
    # job_receiver.clean_up('0000000000', 'whole', [])
    # job_receiver.download_structures(['AF-A0A5E8GAP1-F1-model_v4.pdb', '7M0Q'])
    result_folder = os.path.join(job_receiver.output_path, '9955620089909220348.zip')
    #job_receiver.send_results('8752305630060385451', result_folder)
    job_receiver.compute_file_hash(result_folder)
