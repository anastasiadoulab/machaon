import math
import sys
from src.pdbhandler import PDBHandler
from src.executionhandler import ExecutionHandler
from src.evaluator import Evaluator
from src.evaluator import StructureAligner

sys.path.append('..')
from collections import defaultdict

import copy
import pandas as pd
import os
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

def compute_tm_scores(reference_pdb_id, output_path, pdb_dataset_path, structure_ids):
    print('Calculating 3D similarity of false negatives...')
    current_working_directory = os.getcwd()
    os.chdir(os.path.join(current_working_directory, '../..'))
    execution_handler = ExecutionHandler(31, 12 * 60)
    input_entries = [(reference_pdb_id, x.split('_')[0], pdb_dataset_path) for x in structure_ids]
    # compute_tm(input_entries[0])
    results = execution_handler.parallelize(compute_tm, input_entries)
    os.chdir(current_working_directory)
    output = [x[0] for x in results]
    pd.DataFrame(list(zip(structure_ids, output, [x[1] for x in results])),
                 columns=['candidate', '3D', '3D-cand']).sort_values(by='3D', ascending=False).to_csv(output_path, sep='\t', index=False)
    return output


def compute_tm(entry):
    reference_id, candidate_id, pdb_dataset_path = entry
    evaluator = Evaluator('biopython', False)
    evaluator.reference_data = {'pdbId' : reference_id, 'chainIndex' : 0, 'fullPDBPath': os.path.join(pdb_dataset_path, ''.join([reference_id, '.pdb']))}
    evaluator.candidate_data = {'pdbId' : candidate_id, 'chainIndex' : 0, 'fullPDBPath': os.path.join(pdb_dataset_path, ''.join([candidate_id, '.pdb']))}
    evaluator.set_pdb_path(pdb_dataset_path)
    result = evaluator.calculate_tm_score()
    return result[0], result[1]


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
        _, secondary_structure = aligner.get_structures_by_stride(candidate_id, chain_id, available_residues, stride_output_path)
    if(len(secondary_structure) > 0):
        score, alignment, alignment_result = aligner.align(reference_struct, secondary_structure, '2D', {'gapChars' : ['{', ' ']})
        if(alignment is not False):
            identity, no_gap_identity, gaps = aligner.calculate_identity(alignment_result)
    return (candidate_id, identity)

def compute_sec_struct_identities(config):
    reference_id, pdb_dataset_path, stride_output_path, output_path, structure_ids = config
    pdbhandler = PDBHandler()
    aligner = StructureAligner()
    aligner.pdb_dataset_path = pdb_dataset_path
    current_working_directory = os.getcwd()
    full_pdb_file_path = os.path.sep.join([pdb_dataset_path, ''.join([reference_id, '.pdb'])])
    available_residues = pdbhandler.get_residue_range(full_pdb_file_path, 'A')
    if(len(available_residues['fullRange']) > 1):
        if (os.path.exists(output_path) is False):
            os.chdir(os.path.join(current_working_directory, '../..'))
            _, secondary_structure = aligner.get_structures_by_stride(reference_id, 'A', available_residues, stride_output_path)
            print('Calculating 2D identity with all the protein set...')
            execution_handler = ExecutionHandler(31, 4 * 60)
            input_entries = [(secondary_structure, x.split('_')[0], pdb_dataset_path, stride_output_path) for x in structure_ids]
            results = execution_handler.parallelize(compute_sec_struct_identity, input_entries)
            # results = []
            # input_entries = input_entries[1185:]
            # for entry in input_entries:
            #     results.append(compute_sec_struct_identity(entry))
            os.chdir(current_working_directory)
            results = [x[1] for x in results]
            pd.DataFrame(list(zip(structure_ids, results)),
                         columns=['candidate', '2D-identity']).sort_values(by='2D-identity', ascending=False).to_csv(output_path, sep='\t', index=False)
        else:
            results = pd.read_csv(output_path, sep='\t')
            results = results['2D-identity'].to_list()
    return results

if __name__ == '__main__':
    conformations = defaultdict(list)
    with open('SHREC2018_ref.cla', 'r', encoding='utf-8') as groundtruth_file:
        line_index = 0
        current_protein = ''
        for line in groundtruth_file:
            if (line_index > 1):
                parts = line.strip()
                if (len(parts) > 0):
                    parts = parts.split(' ')
                    if(len(parts) == 1 and len(current_protein) > 0):
                        conformations[current_protein].append(parts[0])
                    else:
                        current_protein = parts[0]
            line_index += 1

    results_path = os.path.join('../../../output_shrec_new')
    total_entries = 2269
    performance_metrics = ['precision', 'recall', 'f1', 'specificity', 'mcc', 'true_positive', 'false_negative']
    stride_output_path = 'structural-data/STRIDE_SHREC'
    pdb_dataset_path = 'structural-data/PDBs_SHREC18'
    os.makedirs(stride_output_path, exist_ok=True)
    with open('shrec_accuracy.log', 'w', encoding='utf-8') as log_file:
        log_file.write('\t'.join(['filename'] + performance_metrics + ['fp_tm_score', 'tp_tm_score', 'fn_tm_score',
                                                                       'fp_sec_identity', 'tp_sec_identity', 'fn_sec_identity']))
        log_file.write('\n')
        for conformation_index, protein in enumerate(conformations):
            if(conformation_index % 100 == 0):
                print('Entry: ', repr(conformation_index))
            if(len(conformations[protein])> 1):
                for filename in conformations[protein]:
                    file_path = os.path.join(results_path, '_'.join([filename, 'A', 'whole']), 'candidates', ''.join([filename, '_A-merged-noenriched_shrec.csv']))
                    results = pd.read_csv(file_path, sep='\t')
                    results['filename'] = results['filename'].apply(str)
                    true_positive = 0
                    false_positive = 0
                    true_negative = 0
                    false_negative = 0
                    found = []
                    fp_tm_scores = []
                    fn_tm_scores = []
                    fp_sec_identities = []
                    tp_tm_scores = []
                    tp_sec_identities = []
                    fn_sec_identities = []
                    fp_files = []
                    tp_files = []
                    for index, row in results.iterrows():
                        if row['filename'] in conformations[protein]:
                            true_positive += 1
                            found.append(row['filename'])
                            tp_files.append(row['filename'])
                            tp_tm_scores.append(row['3D'])
                        if row['filename'] not in conformations[protein]:
                            false_positive += 1
                            fp_files.append(row['filename'])
                            fp_tm_scores.append(row['3D'])
                    if (false_positive > 0):
                        sec_struct_identities_path = os.path.join(results_path, ''.join([filename, '_A_whole']),
                                                                  '2D-identities-fp.tsv')
                        fp_sec_identities = compute_sec_struct_identities( [filename, pdb_dataset_path, stride_output_path,
                                                                            sec_struct_identities_path, fp_files])
                    if (true_positive > 0):
                        sec_struct_identities_path = os.path.join(results_path, ''.join([filename, '_A_whole']),
                                                                  '2D-identities-tp.tsv')
                        tp_sec_identities = compute_sec_struct_identities( [filename, pdb_dataset_path, stride_output_path,
                                                                            sec_struct_identities_path, tp_files])
                    false_negatives = []
                    for conformation in conformations[protein]:
                        if conformation == filename:
                            continue
                        if conformation not in found:
                            false_negatives.append(conformation)
                    fn_sec_identities = compute_sec_struct_identities([filename, pdb_dataset_path, stride_output_path,
                                                                    sec_struct_identities_path, false_negatives])
                    fn_tm_log_path = os.path.join(results_path, ''.join([filename, '_A_whole']), 'candidates',
                                                  ''.join([filename, '-false_negatives_tm_score.csv']))
                    if(len(false_negatives) > 0):
                        fn_tm_scores = compute_tm_scores(filename, fn_tm_log_path, pdb_dataset_path,false_negatives)
                    fp_sec_identities = [x for x in fp_sec_identities if x >= 0]
                    tp_sec_identities = [x for x in tp_sec_identities if x >= 0]
                    fn_sec_identities = [x for x in fn_sec_identities if x >= 0]
                    fp_tm_scores = [x for x in fp_tm_scores if x >= 0]
                    tp_tm_scores = [x for x in tp_tm_scores if x >= 0]
                    fn_tm_scores = [x for x in fn_tm_scores if x >= 0]

                    true_negative = total_entries - true_positive - len(false_negatives) - false_positive
                    metrics = calculate_accuracy_metrics(true_positive, false_positive, true_negative, false_negative)
                    log_file.write('\t'.join([filename] + [repr(metrics[metric]) for metric in performance_metrics] +
                                             [repr(np.median(fp_tm_scores)) if len(fp_tm_scores) > 0 else '-1',
                                              repr(np.median(tp_tm_scores)) if len(tp_tm_scores) > 0 else '-1',
                                              repr(np.median(fn_tm_scores)) if len(fn_tm_scores) > 0 else '-1',
                                              repr(np.median(fp_sec_identities)) if len(
                                                  fp_sec_identities) > 0 else '-1',
                                              repr(np.median(tp_sec_identities)) if len(
                                                  tp_sec_identities) > 0 else '-1',
                                              repr(np.median(fn_sec_identities)) if len(
                                                  fn_sec_identities) > 0 else '-1']))
                    log_file.write('\n')


