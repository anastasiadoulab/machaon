import sys

sys.path.append('..')
from src.pdbhandler import PDBHandler
from src.configurationmanager import ConfigurationManager
from src.machaon import Machaon
from src.executionhandler import ExecutionHandler
from src.evaluator import Evaluator
from src.evaluator import StructureAligner

import copy
import pandas as pd
import os
import math
import numpy as np

def calculate_accuracy_metrics(true_positive, false_positive, true_negative, false_negative):
    metrics = {}
    metrics['precision'] = true_positive / (true_positive + false_positive)
    metrics['recall'] = true_positive / (true_positive + false_negative)
    denominator = (metrics['precision'] + metrics['recall'])
    metrics['f1'] = 0 if denominator == 0 else (2 * (metrics['precision'] * metrics['recall']) / denominator)
    metrics['specificity'] = true_negative / (true_negative + false_positive)
    metrics['mcc'] = ((true_positive * true_negative) - (false_negative * false_positive))
    metrics['mcc'] = metrics['mcc'] / math.sqrt((true_positive + false_positive) * (true_positive + false_negative) *
                                                (true_negative + false_positive) * (true_negative + false_negative))
    metrics['true_positive'] = true_positive
    metrics['false_negative'] = false_negative
    return metrics

def compute_tm(entry):
    reference_id, candidate_id, pdb_dataset_path = entry
    evaluator = Evaluator(False)
    evaluator.reference_data = {'pdbId' : reference_id, 'chainIndex' : 0}
    evaluator.candidate_data = {'pdbId' : candidate_id, 'chainIndex' : 0}
    evaluator.set_pdb_path(pdb_dataset_path)
    result = evaluator.calculate_tm_score()
    return candidate_id, result[0], result[1]

def compute_sec_struct_identity(entry):
    reference_struct, candidate_id, pdb_dataset_path, stride_output_path = entry
    pdbhandler = PDBHandler()
    aligner = StructureAligner()
    aligner.pdb_dataset_path = pdb_dataset_path
    identity = -1
    chain_id = False
    parts = candidate_id.split('_')
    secondary_structure = ''
    if(len(parts) > 1):
        pdbhandler.structure_id = parts[0]
        chain_id = parts[1]
    else:
        pdbhandler.structure_id = candidate_id
        full_pdb_file_path = os.path.sep.join([pdb_dataset_path, ''.join([candidate_id, '.pdb'])])
        chain_id = pdbhandler.get_protein_chain_ids(full_pdb_file_path)
        if(chain_id is not False):
            chain_id = chain_id[0]
        else:
            print(entry)
    if (chain_id is not False):
        available_residues = pdbhandler.get_residue_range(full_pdb_file_path, chain_id)
        if(len(available_residues['fullRange']) > 1):
            _, secondary_structure = aligner.get_structures_by_stride(candidate_id, chain_id, available_residues, stride_output_path)
    if(len(secondary_structure) > 0):
        score, alignment, alignment_result = aligner.align(reference_struct, secondary_structure, '2D', {'gapChars' : ['{', ' ']})
        if(alignment is not False):
            identity, no_gap_identity, gaps = aligner.calculate_identity(alignment)
    return (candidate_id, identity)

def compute_tm_scores(reference_pdb_id, output_path, pdb_dataset_path, structure_ids):
    print('Calculating 3D similarity of protein subset...')
    current_working_directory = os.getcwd()
    os.chdir(os.path.join(current_working_directory, '../..'))
    execution_handler = ExecutionHandler(31, 12 * 60)
    input_entries = [(reference_pdb_id, x.split('_')[0], pdb_dataset_path) for x in structure_ids]
    results = execution_handler.parallelize(compute_tm, input_entries)
    os.chdir(current_working_directory)
    output = [x[0] for x in results]
    pd.DataFrame(results, columns=['candidate', '3D', '3D-cand']).to_csv(output_path, sep='\t', index=False)
    return output

def compute_sec_struct_identities(config):
    reference_id, pdb_dataset_path, stride_output_path, structure_ids = config
    aligner = StructureAligner()
    pdbhandler = PDBHandler()
    aligner.pdb_dataset_path = pdb_dataset_path
    current_working_directory = os.getcwd()
    os.chdir(os.path.join(current_working_directory, '../..'))
    pdbhandler.structure_id = reference_id
    full_pdb_file_path = os.path.sep.join([pdb_dataset_path, ''.join([reference_id.replace('_', ''), '.pdb'])])
    chain_id = pdbhandler.get_protein_chain_ids(full_pdb_file_path)[0]
    available_residues = pdbhandler.get_residue_range(full_pdb_file_path, chain_id)
    _, secondary_structure = aligner.get_structures_by_stride(reference_id, chain_id, available_residues, stride_output_path)
    print('Calculating 2D identity with the protein subset...')
    execution_handler = ExecutionHandler(31, 4 * 60)
    input_entries = [(secondary_structure, x.split('_')[0], pdb_dataset_path, stride_output_path) for x in structure_ids]
    results = execution_handler.parallelize(compute_sec_struct_identity, input_entries)

    os.chdir(current_working_directory)
    return results

if __name__ == '__main__':
    config_manager = ConfigurationManager()
    configuration = copy.deepcopy(config_manager._template_config)
    configuration['rootDisk'] = "structural-data"
    configuration['referenceGeneID'] = "0"
    configuration['referenceSequenceLength'] = 1
    configuration['pdbDatasetPath'] = "pdb70b"
    configuration['isReferenceViral'] = False
    configuration['noThirdPartyData'] = True
    pdb_dataset_path = os.path.join(configuration['rootDisk'], configuration['pdbDatasetPath'])
    stride_output_path = os.path.join(configuration['rootDisk'], 'STRIDE_DALI')
    os.makedirs(stride_output_path, exist_ok=True)

    benchmark_data_path = 'DALI-benchmark/'
    os.makedirs(os.path.join(benchmark_data_path, 'struct-identities'), exist_ok=True)
    os.makedirs(os.path.join(benchmark_data_path, '2D-identities'), exist_ok=True)
    tm_scores_dir = os.path.join(benchmark_data_path, 'benchmark/ordered_querywise/TMalign')
    scop_list = pd.read_csv(os.path.join(benchmark_data_path, 'benchmark/scope_139_targets.list'), sep='\t', names=['query', 'family', 'pdbfile'])
    benchmark_dataset = pd.read_csv(os.path.join(benchmark_data_path, 'benchmark/combinetable.pdb70'), sep='\t', names=['pdbfile', 'queries', 'family'])
    pdb_ids = [x for x in benchmark_dataset['pdbfile'].to_list()]

    query_pdb_ids = scop_list['pdbfile'].to_list()
    stride_failed_pdbs = []
    benchmark_output_dir = os.path.join(benchmark_data_path, 'output_pdb70b')
    os.makedirs(os.path.join(benchmark_data_path, 'struct-identities', 'machaon'), exist_ok=True)
    os.makedirs(os.path.join(benchmark_data_path, '2D-identities', 'machaon'), exist_ok=True)
    with open(''.join(['scop_', 'machaon', '_accuracy.log']), 'w', encoding='utf-8') as log_file:
        for query_index, query_pdb in enumerate(query_pdb_ids):
            if(query_index % 10 == 0):
                print('Current query: ', query_index)
            pdb_filename = query_pdb.strip()
            os.makedirs(os.path.join(benchmark_data_path, 'struct-identities'), exist_ok=True)
            os.makedirs(os.path.join(benchmark_data_path, '2D-identities'), exist_ok=True)
            output_filepath = os.path.join(benchmark_output_dir, pdb_filename.replace('_', ''))
            if(os.path.exists(output_filepath) is False):
                output_filepath = os.path.join(benchmark_output_dir, scop_list[(scop_list['pdbfile'] == pdb_filename) | (scop_list['pdbfile'] == query_pdb)]['query'].item())
            struct_logfile_path = os.path.join(benchmark_data_path, 'struct-identities', 'machaon', pdb_filename.replace('_', ''))
            candidates_dir = os.path.join(output_filepath, output_filepath, 'candidates')
            csv_filepath = ''
            for output_filename in os.listdir(candidates_dir):
                if('merged-notenriched.csv' in output_filename):
                    csv_filepath = os.path.join(candidates_dir, output_filename)
                    break
            benchmark_data = pd.read_csv(csv_filepath, sep='\t')
            output_pdbids = [x.split('_')[0] for x in benchmark_data['structureId'].to_list()]
            benchmark_data['candidate'] = output_pdbids
            benchmark_data = benchmark_data[benchmark_data['candidate'].isin(pdb_ids)].head(250)
            if(os.path.exists(struct_logfile_path) is False):
                identities = compute_sec_struct_identities([pdb_filename.replace('_', ''), pdb_dataset_path, stride_output_path,
                                                            benchmark_data['candidate'].to_list()])
                tm_log_path = os.path.join(''.join([benchmark_data_path, 'output_pdb70b', os.path.sep, pdb_filename.replace('_', '')]),
                                           'candidates', ''.join([pdb_filename, '-_tm_score.csv']))
                compute_tm_scores(pdb_filename.replace('_', ''), tm_log_path, pdb_dataset_path,  benchmark_data['candidate'].to_list())
                identities = pd.DataFrame(identities, columns=['candidate', '2D-identity'])
                similarities = pd.read_csv(tm_log_path, sep='\t')
                combined = pd.merge(benchmark_data, similarities, how='left', on=['candidate'])
                combined = pd.merge(combined, identities[['candidate', '2D-identity']], how='left', on=['candidate'])
                combined.to_csv(struct_logfile_path, sep='\t', index=False)
            else:
                combined = pd.read_csv(struct_logfile_path, sep='\t')
            log_file.write('\t'.join([pdb_filename, repr(combined['3D'].median()),
                                      repr(combined[combined['2D-identity'] > 0]['2D-identity'].median()),
                                      repr(len(combined[combined['2D-identity'] < 0].index)),
                                      repr(len(combined.index))]))
            log_file.write('\n')
            stride_failed_pdbs.extend(combined[combined['2D-identity'] < 0]['candidate'].to_list())
    stride_failed_pdbs = list(set(stride_failed_pdbs))
    with open(''.join(['stride_fail.log']), 'w', encoding='utf-8') as stride_fail_log:
        stride_fail_log.write('\n'.join(stride_failed_pdbs))
