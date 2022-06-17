import sys

sys.path.append('..')
from src.pdbhandler import PDBHandler
from src.configurationmanager import ConfigurationManager
from src.machaon import Machaon
from src.executionhandler import ExecutionHandler
from src.evaluator import Evaluator

import copy
import pandas as pd
import os


def compute_tm(entry):
    reference_id, candidate_id, pdb_dataset_path = entry
    evaluator = Evaluator(False)
    evaluator.reference_data = {'pdbId' : reference_id, 'chainIndex' : 0}
    evaluator.candidate_data = {'pdbId' : candidate_id, 'chainIndex' : 0}
    evaluator.pdb_dataset_path = pdb_dataset_path
    result = evaluator.calculate_tm_score()
    return result[0], result[1]

if __name__ == '__main__':
    config_manager = ConfigurationManager()
    configuration = copy.deepcopy(config_manager._template_config)
    configuration['rootDisk'] = "structural-data"
    configuration['referenceGeneID'] = "0"
    configuration['referenceSequenceLength'] = 1
    configuration['pdbDatasetPath'] = "PDBs_CATH"
    configuration['isReferenceViral'] = False
    configuration['noThirdPartyData'] = True
    configuration['isNonRedundantSet'] = True

    machaon_core = Machaon()

    pdbhandler = PDBHandler()
    pdbhandler.root_disk = configuration['rootDisk']

    cath_data = pd.read_csv('cath.tsv', sep='\t')
    cath_data['ca'] = ['.'.join(x.split('.')[:2]) for x in cath_data['cath']]
    cath_data['c'] = [x.split('.')[0] for x in cath_data['cath']]
    selected_entries = cath_data.drop_duplicates(subset='cath', keep='first').copy()
    cath_entries = [(row['filename'], row['cath'], row['ca'], row['c']) for index, row in selected_entries.iterrows()]

    with open('cath.log', 'w', encoding='utf-8') as log_file:
        log_file.write('\t'.join(['filename', 'cath_family', 'c_found_top', 'c_found_all', 'ca_found_top', 'ca_found_all', 'cath_found_top', 'cath_found_all',
                                  'tm_found_top', 'tm_found_all', 'tm_found_top_cand', 'tm_found_all_cand']))
        log_file.write('\n')
        pdb_dataset_path = os.path.join(configuration['rootDisk'], configuration['pdbDatasetPath'])
        for entry in cath_entries:
            filename, cath_family, ca_family, c_family  = entry
            chain_id = filename[-3]
            structure_id = ''.join([filename, '_', chain_id])
            configuration['referencePDBID'] = filename
            configuration['referenceChainID'] = chain_id
            configuration['outputPath'] =  os.path.join('../../../output_cath_new', ''.join([structure_id, '_whole']))
            machaon_core.perform_comparisons([configuration])
            generated_data = pd.read_csv(os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-merged-notenriched.csv'])), sep='\t')
            generated_data['filename'] = [x.split('_')[0] for x in generated_data['structureId']]
            generated_data = generated_data[1:]
            merged_set = pd.merge(generated_data, cath_data, how='inner', on='filename')
            ca_family_occurrences = len(cath_data[cath_data['ca'] == ca_family].index)
            c_family_occurrences = len(cath_data[cath_data['c'] == c_family].index)
            cath_family_occurrences = len(cath_data[cath_data['cath'] == cath_family].index)
            print('Calculating 3D similarity...')
            current_working_directory = os.getcwd()
            os.chdir(os.path.join(current_working_directory,'../..'))
            execution_handler = ExecutionHandler(31, 12 * 60)
            input_entries = [(configuration['referencePDBID'], x.split('_')[0], pdb_dataset_path) for x in merged_set['structureId']]
            results = execution_handler.parallelize(compute_tm, input_entries)
            merged_set['3D'] = [x[0] for x in results]
            merged_set['3D-cand'] = [x[1] for x in results]
            os.chdir(current_working_directory)
            merged_set.to_csv(os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-merged-notenriched_cath.csv'])), index=False, sep='\t')
            merged_set = pd.read_csv(os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-merged-notenriched_cath.csv'])), dtype={7:'object',8:'object'}, sep='\t')
            if(cath_family_occurrences > 1):
                top_rows = merged_set.head(20).copy()
                ca_found_top = (len(top_rows[top_rows['ca'] == ca_family].index) / len(top_rows.index)) * 100
                c_found_top = (len(top_rows[top_rows['c'] == c_family].index) / len(top_rows.index)) * 100
                tm_found_top = top_rows['3D'].median()
                tm_found_top_cand = top_rows['3D-cand'].median()
                merged_set['ca'] = merged_set['ca'].apply(str)
                merged_set['c'] = merged_set['c'].apply(str)
                cath_found_top = (len(top_rows[top_rows['cath'] == cath_family].index) / len(top_rows.index)) * 100
                cath_found_all = (len(merged_set[merged_set['cath'] == cath_family].index) / cath_family_occurrences) * 100
                ca_found_all = (len(merged_set[merged_set['ca'] == ca_family].index) / ca_family_occurrences) * 100
                c_found_all = (len(merged_set[merged_set['c'] == c_family].index) / c_family_occurrences) * 100
                tm_found_all = merged_set['3D'].median()
                tm_found_all_cand = merged_set['3D-cand'].median()
                log_file.write('\t'.join([filename, cath_family, repr(c_found_top),  repr(c_found_all), repr(ca_found_top), repr(ca_found_all),
                                          repr(cath_found_top), repr(cath_found_all),  repr(tm_found_top), repr(tm_found_all),  repr(tm_found_top_cand),
                                          repr(tm_found_all_cand)]))
                log_file.write('\n')





