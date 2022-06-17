
import os
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({'font.size': 8})
import operator as op


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

    aggregated_values = []
    sessions = 0
    for protein in conformations:
        if (len(conformations[protein]) > 1):
            for pdb_file in conformations[protein]:
                sessions += 1
                chain_id = 'A'
                structure_id = ''.join([pdb_file, '_', chain_id])
                result_path =  os.path.join('../../../output_shrec_new', ''.join([structure_id, '_whole']))
                generated_data = pd.read_csv(os.path.join(result_path, 'candidates', ''.join([structure_id, '-merged-noenriched_shrec.csv'])), sep='\t')
                aggregated_values.extend(generated_data['3D'].to_list())
    aggregated_values.sort()
    ax = sns.violinplot(data=aggregated_values, width=0.5,  palette='colorblind')
    plt.ylabel('TM-Score (%)')
    plt.xticks(plt.xticks()[0], [''.join([repr(len(aggregated_values)), ' structures'])])
    plt.title(''.join(['3D similarity of the proteins in the results \n (', repr(sessions), ' search sessions)']))
    plt.savefig('SHREC-searches-3D-violin.png', format='png', dpi=300)
    ax.set_rasterized(True)
    plt.savefig('SHREC-searches-3D-violin.eps', format='eps', dpi=300)



