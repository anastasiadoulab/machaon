import sys

sys.path.append('..')
from src.pdbhandler import PDBHandler
from src.configurationmanager import ConfigurationManager
from src.machaon import Machaon
from src.executionhandler import ExecutionHandler
from src.evaluator import Evaluator
from collections import defaultdict

import copy
import pandas as pd
import os


def compute_tm(entry):
    reference_id, candidate_id, pdb_dataset_path = entry
    evaluator = Evaluator('biopython', False)
    evaluator.reference_data = {'pdb_id' : reference_id, 'chainIndex' : 0, 'fullPDBPath': os.path.join(pdb_dataset_path, ''.join([reference_id, '.pdb']))}
    evaluator.candidate_data = {'pdb_id' : candidate_id, 'chainIndex' : 0, 'fullPDBPath': os.path.join(pdb_dataset_path, ''.join([candidate_id, '.pdb']))}
    evaluator.pdb_dataset_path = pdb_dataset_path
    result = evaluator.calculate_tm_score()
    return result[0], result[1]


if __name__ == '__main__':
    config_manager = ConfigurationManager()
    configuration = copy.deepcopy(config_manager._template_config)
    configuration['rootDisk'] = "structural-data"
    configuration['referenceGeneID'] = "0" 
    configuration['pdbDatasetPath'] = "PDBs_SHREC18"
    configuration['isReferenceViral'] = False
    configuration['noThirdPartyData'] = True
    configuration['isNonRedundantSet'] = True

    comparison = Machaon()

    pdbhandler = PDBHandler()
    pdbhandler.root_disk = configuration['rootDisk']

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

    with open('shrec.log', 'w', encoding='utf-8') as log_file:
        log_file.write('\t'.join(['filename', 'chain_length', 'conformations_found']))
        log_file.write('\n')
        pdb_dataset_path = os.path.join(configuration['rootDisk'], configuration['pdbDatasetPath'])
        for protein in conformations:
            if(len(conformations[protein])> 1):
                for pdb_file in conformations[protein]:
                    chain_id = 'A'
                    structure_id = ''.join([pdb_file, '_', chain_id])
                    configuration['referencePDBID'] = pdb_file
                    configuration['referenceChainID'] = chain_id
                    configuration['outputPath'] =  os.path.join('../../../output_shrec_new', ''.join([structure_id, '_whole']))
                    comparison.perform_comparisons([configuration])
                    generated_data = pd.read_csv(os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-merged-notenriched.csv'])), sep='\t')
                    generated_data['filename'] = [x.split('_')[0] for x in generated_data['structureId']]
                    generated_data = generated_data[1:]
                    residue_range = pdbhandler.get_residue_range(os.path.join(configuration['rootDisk'], configuration['pdbDatasetPath'], '.'.join([pdb_file, 'pdb'])), chain_id)
                    length = len(residue_range['fullRange'])

                    current_working_directory = os.getcwd()
                    os.chdir(os.path.join(current_working_directory,'../..'))
                    execution_handler = ExecutionHandler(31, 12 * 60)
                    input_entries = [(configuration['referencePDBID'], filename, pdb_dataset_path) for filename in generated_data['filename']]
                    results = execution_handler.parallelize(compute_tm, input_entries)
                    generated_data['3D'] = [x[0] for x in results]
                    generated_data['3D-cand'] = [x[1] for x in results]
                    os.chdir(current_working_directory)
                    generated_data.to_csv(os.path.join(configuration['outputPath'], 'candidates', ''.join([structure_id, '-merged-noenriched_shrec.csv'])), index=False, sep='\t')
