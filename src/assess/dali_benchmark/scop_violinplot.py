
import os
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({'font.size': 8})
import operator as op


if __name__ == '__main__':
    aggregated_values = {'3D' : [], '2D-identity' : []}
    sessions = 0
    result_path = 'output_path'
    for filename in os.listdir(result_path):
        sessions += 1
        generated_data = pd.read_csv(os.path.join(result_path, filename), sep='\t')
        for key in aggregated_values:
            aggregated_values[key].extend([x for x in generated_data[key].to_list() if x >= 0])

for key in aggregated_values:
    aggregated_values[key].sort()
    ax = sns.violinplot(data=aggregated_values[key], width=0.5,  palette='colorblind')
    plt.ylabel('TM-Score (%)' if key == '3D' else '2D-identity (%')
    plt.xticks(plt.xticks()[0], [''.join([repr(len(aggregated_values[key])), ' structures'])])
    plt.title(''.join(['3D similarity' if key == '3D' else '2D identity', ' of the proteins in the final cluster \n (', repr(sessions), ' search sessions)']))
    plt.savefig(''.join(['SCOP-searches-', '3D' if key == '3D' else '2D',  '-violin.png']), format='png', dpi=300)
    ax.set_rasterized(True)
    plt.savefig(''.join(['SCOP-searches-', '3D' if key == '3D' else '2D', '-violin.eps']), format='eps', dpi=300)
    plt.clf()



