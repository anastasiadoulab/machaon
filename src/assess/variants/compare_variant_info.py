import os
from collections import defaultdict
import pandas as pd

structure_ids = ['6VXX_A', '6VXXp_A', '7T9K_A', '7V7Q_A']

output_path = '../../../output_fix'

with open('go_diff.log', 'w') as log_file:
    properties = ['molecularFunction', 'cellularComponent', 'biologicalProcess']
    for property in properties:
        virus_go_info = defaultdict(list)
        for structure_id in structure_ids:
            file_path = os.path.join(output_path, '_'.join([structure_id, 'whole']), 'go', ''.join([structure_id, '_', property, '-go-aggregate.csv']))
            virus_go_info[structure_id] = pd.read_csv(file_path, sep='\t')

        log_file.write(property)
        log_file.write('\n\n7T9K_A\n\n')
        log_file.write(virus_go_info['7T9K_A'][~virus_go_info['7T9K_A']['propertyName'].isin(virus_go_info['6VXX_A']['propertyName'])].to_string())
        log_file.write('\n\n7V7Q_A\n\n')
        log_file.write(virus_go_info['7V7Q_A'][~virus_go_info['7V7Q_A']['propertyName'].isin(virus_go_info['6VXX_A']['propertyName'])].to_string())
        log_file.write('\n\n\n')


    for structure_index, structure_id in enumerate(structure_ids):
        for struct_index in range(structure_index+1, len(structure_ids)):
            x_results = pd.read_csv(os.path.join(output_path, '_'.join([structure_id, 'whole']), 'candidates',
                                                       '-'.join([structure_id, 'merged-enriched_eval_full.csv'])), sep='\t')
            x_results['composite_id'] = x_results.apply(lambda row: '_'.join([row['pdbId'], row['chainId']]),
                                                                    axis=1)
            x_results = x_results[x_results['composite_id'] != structure_id]

            file_path = os.path.join(output_path, '_'.join([structure_ids[struct_index], 'whole']), 'candidates', '-'.join([structure_ids[struct_index], 'merged-enriched_eval_full.csv']))
            y_results = pd.read_csv(file_path, sep='\t')
            y_results['composite_id'] = y_results.apply(lambda row: '_'.join([row['pdbId'], row['chainId']]), axis=1)
            y_results = y_results[y_results['composite_id'] != structure_ids[struct_index]]

            log_file.write(''.join([structure_id,' | ', 'uniprotId', '\n\n']))
            x_ref = x_results['uniprotId'].to_list()
            y_ref = y_results['uniprotId'].to_list()
            log_file.write(''.join(['\nnot common to ', structure_ids[struct_index], '\n']))
            log_file.write('\n'.join([x for x in x_ref if x not in y_ref]))
            log_file.write(''.join(['\nnot common to ', structure_id, '\n']))
            log_file.write('\n'.join([x for x in y_ref if x not in x_ref]))
            common = len([x for x in x_ref if x in y_ref])
            log_file.write(''.join(['\ncommon: ', repr(common), '\n\n']))
            y_ref.extend(x_ref)
            superset_length = len(list(set(y_ref)))
            log_file.write(''.join(['superset length: ', repr(superset_length), '\n']))
            log_file.write(''.join(['percentage of superset that is common: ', repr(common / superset_length), '\n']))
            log_file.write('\n\n\n')