import pandas as pd
import numpy as np
import os
import time
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 8})
from scipy import stats
import seaborn as sns

def plot2d(dataValues, xlabel, ylabel, title, outputPath):
    plt.plot(xlabel, ylabel, data=dataValues, linestyle='', marker='o', markersize=0.7)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title, loc='left')


rootDisk = 'structural-data'

pdbDatasetPath = os.path.sep.join([rootDisk, 'PDBs_vir'])

start = time.time()


referenceStructureIds = ['6VXX']
referenceChains = ['A']

sections = ['1D-identity', '2D-identity', '5UTR-identity', 'CDS-identity', '3UTR-identity', '3D-score', 'chemSim']
#sections = ['2D-identity']

metrics = ['b-phipsi', 'w-rdist', 't-alpha']
metricOrder = [True, True, True]
color_palette = sns.color_palette('colorblind', len(sections))
metrics_directory = '6VXX_A_whole/metrics/'
for refIndex, referenceStructureId in enumerate(referenceStructureIds):
    referenceChain = referenceChains[refIndex]
    outputPath = ''.join( [metrics_directory, referenceStructureId, '_', str(referenceChain), '-metrics'])
    figure, axes = plt.subplots(3)
    for metricIndex, metric in enumerate(metrics):
        for sectionIndex, section in enumerate(sections):
            data = []
            evaluations = []
            print(''.join(['=======================', referenceStructureId, '_', str(referenceChain), ' ', section, ' ', metric]))
            data = pd.read_csv(''.join([metrics_directory, referenceStructureId, '_', str(referenceChain), '_', metric, '.csv']), sep=',')
            evaluations = pd.read_csv(''.join([metrics_directory, referenceStructureId, '_', str(referenceChain),'_', metric, '_eval.csv']), sep='\t')
            evaluations = evaluations[evaluations[section] >= 0]
            evaluations['composite_id'] = evaluations.apply(lambda row: '_'.join([row['pdbId'], row['chainId']]), axis=1)
            data['composite_id'] = data['structureId']
            mergedData = pd.merge(evaluations[['composite_id', section]], data[['composite_id', metrics[metricIndex]]], how='inner', on=['composite_id'])
            mergedData = mergedData.replace([np.inf, -np.inf], np.nan).dropna()
            threshold = 4 if metrics[metricIndex] != 't-alpha' else 0
            threshold = 0
            mergedData = mergedData[stats.zscore(mergedData[metrics[metricIndex]]) < threshold]
            mergedData[section] = mergedData[section] * 100
            mergedData[section] = mergedData[section].ewm(span=400).mean().values
            lineLabel = section.replace('-identity', '').replace('-score', '').replace('chemSim', 'chemical')
            axes[metricIndex].plot(mergedData[[metrics[metricIndex]]], mergedData[[section]], color=color_palette[sectionIndex], label=lineLabel)
            axes[metricIndex].set_xlabel(metrics[metricIndex])

    plt.subplots_adjust(hspace=0.5)
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc='lower right', fancybox=True, framealpha=1.0, bbox_to_anchor=(0.92, 0.05))
    figure.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.ylabel('Percentage (%)')
    plt.savefig(''.join([outputPath, '.png']), format='png', dpi=300, bbox_inches='tight')
    plt.savefig(''.join([outputPath, '.eps']), format='eps', dpi=300, bbox_inches='tight')
