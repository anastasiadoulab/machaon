import os
import traceback
import subprocess
import pandas as pd
from src.evaluator import Evaluator
import time
from src.executionhandler import ExecutionHandler

def evaluate(entry):
    referenceProtein, idParts, isViral = entry
    evaluator = Evaluator(isViral)
    rootDisk = 'structural-data'
    evaluator.set_root_disk(rootDisk)
    pdbDatasetPath = 'PDBs_vir'
    evaluator.set_pdb_path(os.path.sep.join([rootDisk, pdbDatasetPath]))
    evaluator.reference_data = evaluator.load_data(*referenceProtein)
    evaluator.candidate_data = evaluator.load_data(*idParts)
    result = False
    try:
        sections = ['1D', '2D', '5UTR', 'CDS', '3UTR']
        compared = evaluator.align_sequences(sections)
        compared['3D'] = evaluator.calculate_tm_score()
        output = [idParts[0], idParts[1]]
        sections.append('3D')
        for section in sections:
            if (section in compared):
                output.extend([repr(x) for x in compared[section]])
            else:
                output.extend([-1, -1, -1])  # every section has a triplet of values
        output.append(repr(evaluator.compute_chemical_similarity()))
        result = '\t'.join(output)
    except:
        print(''.join(['Failed: ', repr(idParts)]))
        print(traceback.format_exc())
    return result

referenceProteins = [ ('6VXX', 'A')]
os. chdir('..')
for referenceProtein in referenceProteins:
    metrics = ['b-phipsi', 'w-rdist', 't-alpha']
    resultPath = '6VXX_A_whole/metrics'
    for metric in metrics:
        inputEntries = []
        data = pd.read_csv(''.join([resultPath, os.path.sep, referenceProtein[0], '_', referenceProtein[1], '_', metric, '.csv']))
        results = ['\t'.join(['pdbId', 'chainId',  '1D-score', '1D-identity', '1D-identity-ng', '1D-identity-gaps', '2D-score',
                               '2D-identity', '2D-identity-ng', '2D-identity-gaps', '5UTR-score', '5UTR-identity', '5UTR-identity-ng',
                               '5UTR-identity-gaps', 'CDS-score', 'CDS-identity', 'CDS-identity-ng', 'CDS-identity-gaps', '3UTR-score',
                               '3UTR-identity','3UTR-identity-ng', '3UTR-identity-gaps', '3D-score', '3D-score-cand', 'length', 'chemSim'])]
        for index, row in data.iterrows():
            if(index % 5 != 0):
                continue
            inputEntries.append([referenceProtein, row['structureId'].split('_'), True])

        execution_handler = ExecutionHandler(18, 5 * 60)
        results.extend(execution_handler.parallelize(evaluate, inputEntries))
        with open(''.join([resultPath, os.path.sep, referenceProtein[0], '_', referenceProtein[1], '_', metric, '_eval.csv']), 'w') as outfile:
            outfile.write('\n'.join(results))
print(' '.join(['==================================[COMPLETED: ', time.ctime(), ']'])) 
