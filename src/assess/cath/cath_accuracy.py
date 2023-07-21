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
    evaluator = Evaluator('biopython', False)
    evaluator.reference_data = {'pdbId' : reference_id, 'chainIndex' : 0}
    evaluator.candidate_data = {'pdbId' : candidate_id, 'chainIndex' : 0}
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

def compute_tm_scores(reference_pdb_id, output_path, pdb_dataset_path, structure_ids):
    print('Calculating 3D similarity of false negatives...')
    current_working_directory = os.getcwd()
    os.chdir(os.path.join(current_working_directory, '../..'))
    execution_handler = ExecutionHandler(31, 12 * 60)
    input_entries = [(reference_pdb_id, x.split('_')[0], pdb_dataset_path) for x in structure_ids]
    results = execution_handler.parallelize(compute_tm, input_entries)
    os.chdir(current_working_directory)
    output = [x[0] for x in results]
    pd.DataFrame(list(zip(structure_ids, output, [x[1] for x in results])),
                 columns=['candidate', '3D', '3D-cand']).sort_values(by='3D', ascending=False).to_csv(output_path, sep='\t', index=False)
    return output

def compute_sec_struct_identities(config):
    reference_id, pdb_dataset_path, stride_output_path, output_path, structure_ids = config
    pdbhandler = PDBHandler()
    aligner = StructureAligner()
    aligner.pdb_dataset_path = pdb_dataset_path
    current_working_directory = os.getcwd()
    reference_pdb_id, chain_id = reference_id.split('_')
    full_pdb_file_path = os.path.sep.join([pdb_dataset_path, ''.join([reference_pdb_id, '.pdb'])])
    available_residues = pdbhandler.get_residue_range(full_pdb_file_path, chain_id)
    results = []
    if(len(available_residues['fullRange']) > 1):
        if(os.path.exists(output_path) is False):
            os.chdir(os.path.join(current_working_directory, '../..'))
            _, secondary_structure = aligner.get_structures_by_stride(reference_pdb_id, chain_id, available_residues, stride_output_path)
            print('Calculating 2D identity with the protein subset...')
            execution_handler = ExecutionHandler(31, 4 * 60)
            input_entries = [(secondary_structure, x.split('_')[0], pdb_dataset_path, stride_output_path) for x in structure_ids]
            results = execution_handler.parallelize(compute_sec_struct_identity, input_entries)
            os.chdir(current_working_directory)
            results = [x[1] for x in results]
            pd.DataFrame(list(zip(structure_ids, results)), columns=['candidate', '2D-identity']).sort_values(by='2D-identity',
                                                                    ascending=False).to_csv(output_path, sep='\t', index=False)
        else:
            data = pd.read_csv(output_path, sep='\t')
            results = data['2D-identity'].to_list()
            if(len(results) > 0):
                if(isinstance(results[0], str)):
                    if( "('" in results[0]):
                        results = [float(x.split(',')[1].strip().replace(')', '')) for x in results]
                        data['2D-identity'] = results
                        data.to_csv(output_path, sep='\t', index=False)
    return results

if __name__ == '__main__':
    config_manager = ConfigurationManager()
    configuration = copy.deepcopy(config_manager._template_config)
    configuration['rootDisk'] = "structural-data"
    configuration['referenceGeneID'] = "0" 
    configuration['pdbDatasetPath'] = "PDBs_CATH"
    configuration['isReferenceViral'] = False
    configuration['noThirdPartyData'] = True
    pdb_dataset_path = os.path.join(configuration['rootDisk'], configuration['pdbDatasetPath'])
    stride_output_path = os.path.join(configuration['rootDisk'], 'STRIDE_CATH')
    os.makedirs(stride_output_path, exist_ok=True)

    cath_data = pd.read_csv('cath.tsv', sep='\t')
    cath_data['ca'] = ['.'.join(x.split('.')[:2]) for x in cath_data['cath']]
    cath_data['c'] = [x.split('.')[0] for x in cath_data['cath']]
    selected_entries = cath_data.drop_duplicates(subset='cath', keep='first').copy()
    cath_entries = [(row['filename'], row['cath'], row['ca'], row['c']) for index, row in selected_entries.iterrows()]
    total_entries = 2685

    performance_metrics = ['precision', 'recall', 'f1', 'specificity', 'mcc', 'true_positive', 'false_negative']
    with open('cath_accuracy.log', 'w', encoding='utf-8') as log_file:
        log_file.write('\t'.join(['filename'] + performance_metrics + ['fp_tm_score', 'tp_tm_score', 'fn_tm_score',
                                                                       'fp_sec_identity', 'tp_sec_identity', 'fn_sec_identity']))
        log_file.write('\n')

        for entry_index, entry in enumerate(cath_entries):
            if(entry_index % 100 == 0):
                print('Entry: ', repr(entry_index))
            filename, cath_family, ca_family, c_family  = entry
            chain_id = filename[-3]
            structure_id = ''.join([filename, '_', chain_id])
            configuration['referencePDBID'] = filename
            configuration['referenceChainID'] = chain_id
            configuration['outputPath'] =  os.path.join('../../../output_cath_new', ''.join([structure_id, '_whole']))

            cath_family_occurrences = len(cath_data[cath_data['cath'] == cath_family].index)
            merged_set = pd.read_csv(os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-merged-notenriched_cath.csv'])), dtype={'ca':'object','c':'object'}, sep='\t')
            structure_ids = [x.split('_')[0] for x in merged_set['structureId'].to_list()]
            fp_tm_scores = []
            tp_tm_scores = []
            fp_sec_identities = []
            tp_sec_identities = []
            if(cath_family_occurrences > 1):
                tm_log_path = os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-_tm_score.csv']))
                #compute_tm_scores(configuration['referencePDBID'], tm_log_path, pdb_dataset_path, cath_data['filename'].to_list())
                merged_set['ca'] = merged_set['ca'].apply(str)
                true_positives = len(merged_set[merged_set['cath'] == cath_family].index)
                false_positives = len(merged_set[merged_set['cath'] != cath_family].index)
                excluded_set = cath_data[~cath_data['filename'].isin(structure_ids)]
                assert len(merged_set.index) + len(excluded_set.index) == total_entries
                false_negative_set = excluded_set[(excluded_set['cath'] == cath_family) & (excluded_set['filename'] != configuration['referencePDBID'])]
                false_negatives = len(false_negative_set.index)
                fn_tm_scores = []
                fn_sec_identities = []
                fn_tm_log_path = os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-false_negatives_tm_score.csv']))
                if(os.path.exists(fn_tm_log_path)):
                    os.unlink(fn_tm_log_path)
                if(false_negatives > 0):
                    fn_tm_scores = compute_tm_scores(configuration['referencePDBID'], fn_tm_log_path, pdb_dataset_path, false_negative_set['filename'].to_list())
                    sec_struct_identities_path = os.path.join('../../../output_cath_new',
                                                              ''.join([structure_id, '_whole']), '2D-identities-fn.tsv')
                    fn_sec_identities = compute_sec_struct_identities(
                        [structure_id, pdb_dataset_path, stride_output_path, sec_struct_identities_path,
                         false_negative_set['filename'].to_list()])
                true_negative = total_entries - len(merged_set.index) - false_negatives
                if(false_positives > 0):
                    fp_tm_scores.extend(merged_set[merged_set['cath'] != cath_family]['3D'].to_list())
                    sec_struct_identities_path = os.path.join('../../../output_cath_new',
                                                              ''.join([structure_id, '_whole']), '2D-identities-fp.tsv')
                    fp_sec_identities.extend(compute_sec_struct_identities(
                        [structure_id, pdb_dataset_path, stride_output_path, sec_struct_identities_path,
                         merged_set[merged_set['cath'] != cath_family]['filename'].to_list()]))
                if(true_positives > 0):
                    tp_tm_scores.extend(merged_set[merged_set['cath'] == cath_family]['3D'].to_list())
                    sec_struct_identities_path = os.path.join('../../../output_cath_new',
                                                              ''.join([structure_id, '_whole']), '2D-identities-tp.tsv')
                    tp_sec_identities.extend(compute_sec_struct_identities(
                        [structure_id, pdb_dataset_path, stride_output_path, sec_struct_identities_path,
                         merged_set[merged_set['cath'] == cath_family]['filename'].to_list()]))
                if (false_positives > 0 and false_negatives > 0):
                    fp_rp = []
                    for fp_tm_score in fp_tm_scores:
                        for fn_tm_score in fn_tm_scores:
                            if(fp_tm_score > fn_tm_score):
                                fp_rp.append(fp_tm_score)
                                break
                    if(len(fp_rp) > 0):
                        print('Greater tm scores found in fp vs fn: ', len(fp_rp))
                if (false_positives > 0 and true_positives > 0):
                    fp_rp = []
                    for fp_tm_score in fp_tm_scores:
                        for tp_tm_score in tp_tm_scores:
                            if(fp_tm_score > tp_tm_score):
                                fp_rp.append(fp_tm_score)
                                break
                    if(len(fp_rp) > 0):
                        print('Greater tm scores found in fp vs tp: ', len(fp_rp))
                metrics = calculate_accuracy_metrics(true_positives, false_positives, true_negative, false_negatives)
                fp_sec_identities = [x for x in fp_sec_identities if x >= 0]
                tp_sec_identities = [x for x in tp_sec_identities if x >= 0]
                fn_sec_identities = [x for x in fn_sec_identities if x >= 0]
                fp_tm_scores = [x for x in fp_tm_scores if x >= 0]
                tp_tm_scores = [x for x in tp_tm_scores if x >= 0]
                fn_tm_scores = [x for x in fn_tm_scores if x >= 0]

                log_file.write('\t'.join([filename] + [repr(metrics[metric]) for metric in performance_metrics] + [repr(np.median(fp_tm_scores)) if len(fp_tm_scores) > 0 else '-1',
                                                                                                                   repr(np.median(tp_tm_scores)) if len(tp_tm_scores) > 0 else '-1',
                                                                                                                   repr(np.median(fn_tm_scores)) if len(fn_tm_scores) > 0 else '-1',
                                                                                                                   repr(np.median(fp_sec_identities)) if len(fp_sec_identities) > 0 else '-1',
                                                                                                                   repr(np.median(tp_sec_identities)) if len(tp_sec_identities) > 0 else '-1',
                                                                                                                   repr(np.median(fn_sec_identities)) if len(fn_sec_identities) > 0 else '-1']))
                log_file.write('\n')



