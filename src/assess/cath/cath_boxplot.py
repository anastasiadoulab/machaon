
import os
from collections import defaultdict
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
plt.rcParams.update({'font.size': 8})
import operator as op


if __name__ == '__main__':
    cath_data = pd.read_csv('cath.tsv', sep='\t')
    cath_data['ca'] = ['.'.join(x.split('.')[:2]) for x in cath_data['cath']]
    cath_data['c'] = [x.split('.')[0] for x in cath_data['cath']]
    selected_entries = cath_data.drop_duplicates(subset='cath', keep='first').copy()
    cath_entries = [(row['filename'], row['cath'], row['ca'], row['c']) for index, row in selected_entries.iterrows()]

    results = defaultdict(list)
    for entry in cath_entries:
        filename, cath_family, ca_family, c_family  = entry
        chain_id = filename[-3]
        structure_id = ''.join([filename, '_', chain_id])
        result_path =  os.path.join('../../../output_cath_new', ''.join([structure_id, '_whole']))
        generated_data = pd.read_csv(os.path.join(result_path, 'candidates', ''.join([structure_id, '-merged-notenriched_cath.csv'])), sep='\t')
        generated_data['filename'] = [x.split('_')[0] for x in generated_data['structureId']]
        generated_data = generated_data[1:]
        results[structure_id] = generated_data['3D'].to_list()

    sorted_keys, sorted_values = zip(*sorted(results.items(), key=op.itemgetter(1)))
    source_data = pd.DataFrame.from_dict(results, orient='index').to_csv('CATH-searches-3D-box.csv', index=False)
    flierprops = dict(markerfacecolor='None', markersize=0.1, markeredgecolor='gray')
    ax = sns.boxplot(data=sorted_values, palette='colorblind', width=1.0, linewidth=0.3, flierprops=flierprops)
    plt.xticks(plt.xticks()[0], sorted_keys)
    plt.ylabel('TM-Score (%)')
    plt.xlabel('Reference proteins')
    plt.tick_params(axis='x', labelsize=3, labelrotation=90)
    plt.title('3D similarity of the proteins in the results')
    plt.savefig('CATH-searches-3D-box.png', format='png', dpi=300)
    ax.set_rasterized(True)
    plt.savefig('CATH-searches-3D-box.eps', format='eps', dpi=300)



