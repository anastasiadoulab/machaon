import pandas as pd
import numpy as np
import os
import time
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size' : 16})
from scipy import stats

def plot2d(dataValues, title, outputPath):
    plt.clf()
    plt.plot(dataValues, linestyle='', marker='o', markersize=0.7)
    plt.xlabel('metric')
    plt.ylabel('occurences')
    plt.title(title, loc='left')
    plt.savefig(outputPath)

def drawHistogram(dataValues, title, outputPath, logvalues=True):
    plt.clf()
    _ = plt.hist(dataValues, bins='fd')
    plt.title(title)
    if(logvalues is True):
        plt.yscale('log')
        plt.ylabel('occurences (log)')
    else:
        plt.ylabel('occurences')
    plt.xlabel('values')
    plt.tight_layout()
    plt.savefig(outputPath, format='png', dpi=300)
    plt.savefig(outputPath.replace('.png', '.eps'), format='eps', dpi=300)

rootDisk = 'structural-data'


pdbDatasetPath = os.path.sep.join([rootDisk, 'PDBs_vir'])

start = time.time()




referenceStructureIds = ['6VXX']
referenceChains = ['A']

sections = ['1D-score', '1D-identity', '1D-identity-ng', '2D-score',
           '2D-identity', '2D-identity-ng', '5UTR-score', '5UTR-identity', '5UTR-identity-ng',
           'CDS-score', 'CDS-identity', 'CDS-identity-ng', '3UTR-score', '3UTR-identity',
           '3UTR-identity-ng', '3D-score']

metrics = ['b-phipsi', 'w-rdist', 't-alpha']
metricOrder = [True, True, True]
metricValues = {}

for refIndex, referenceStructureId in enumerate(referenceStructureIds):
        results = []
        referenceChain = referenceChains[refIndex]
        data = []
        evaluations = []
        for metricIndex, metric in enumerate(metrics):
            metricValues[metrics[metricIndex]] = {}
            print(''.join(['=======================', referenceStructureId, '_', str(referenceChain), ' ',  metric]))
            data = pd.read_csv(''.join(['6VXX_A_whole/metrics/', referenceStructureId, '_', str(referenceChain), '_', metric, '.csv']), sep=',')
            outputPath = ''.join(['6VXX_A_whole/metrics/', referenceStructureId, '_', str(referenceChain), '_', metric, '.png'])
            data = data.replace([np.inf, -np.inf], np.nan).dropna()

            for index, row in data.iterrows():
                metricValues[metrics[metricIndex]][row['structureId']] = row[metrics[metricIndex]]
            toPlot = data[metrics[metricIndex]].values
            zscores = stats.zscore(toPlot)
            indices = np.where(zscores < 3.5)
            print(' '.join([metrics[metricIndex], 'min', repr(np.min(toPlot)), 'max', repr(np.max(toPlot)), 'median', repr(np.median(toPlot)), 'std', repr(np.std(toPlot))]))
            #toPlot = toPlot[indices]
            drawHistogram(toPlot, metrics[metricIndex], outputPath)

fullOccurrences = []
for i in range(len(metrics)):
    fullOccurrences.append([])

for structureid in metricValues[metrics[0]]:
    computed = True
    for i in range(1, len(metrics)):
        computed = computed and structureid in metricValues[metrics[i]]
    if(computed):
        for metricIndex in range(len(metrics)):
            fullOccurrences[metricIndex].append(metricValues[metrics[metricIndex]][structureid])

for i in range(len(metrics)):
    fullOccurrences[i] = np.array(fullOccurrences[i])
metricValues = np.array(fullOccurrences)
covMatrix = np.cov(metricValues)
print(covMatrix)
print(np.linalg.eig(covMatrix))


df = pd.DataFrame(metricValues.reshape(-1, len(metrics)), columns=metrics)
print(df.corr())
