import pandas as pd
root = 'output'

#filenames = ['6VXX_A_segment/candidates/6VXX_A_site0-metrics-merged-enriched_eval.csv', '6VXX_A_segment/candidates/6VXX_A_site1-metrics-merged-enriched_eval.csv',
             # '6VXX_A_segment/candidates/6VXX_A_site2-metrics-merged-enriched_eval.csv', '6VXX_A_segment/candidates/6VXX_A_site3-metrics-merged-enriched_eval.csv',
             # '6VXX_A_segment/candidates/6VXX_A_site4-metrics-merged-enriched_eval.csv']
filenames = ['6VXX_A_domain/candidates/6VXX_A_BetaCoV-S1-CTD-merged-enriched_eval.csv', '6VXX_A_domain/candidates/6VXX_A_BetaCoV-S1-NTD-merged-enriched_eval.csv',
            '6VXX_A_domain/candidates/6VXX_A_RBD-merged-enriched_eval.csv']
# filenames = ['6VXX_A_whole/candidates/6VXX_A-merged-enriched_eval.csv']
columnHeaders = ['b-phipsi', 'w-rdist', 't-alpha']
metricOrder = [True, True, True]
aggregatedData = None

evalMetrics = ['b-phipsi', 'w-rdist', 't-alpha','1D-identity',
'1D-identity-ng', '1DPDB-identity',
'1DPDB-identity-ng', '2D-identity', '2D-identity-ng',
'5UTR-identity', '5UTR-identity-ng',
'CDS-identity', 'CDS-identity-ng',
'3UTR-identity', '3UTR-identity-ng', '3D-score', 'chemSim', 'length']


single = True

for fileName in filenames:
    if(aggregatedData is None):
        aggregatedData = pd.read_csv(''.join([root, fileName]), sep='\t')
    else:
        aggregatedData.append(pd.read_csv(''.join([root, fileName]), sep='\t'))

aggregatedData = aggregatedData.drop_duplicates()
evalMetricData = aggregatedData[evalMetrics]
evalMetricData = evalMetricData[evalMetricData >= 0]
evalMetricData = evalMetricData.dropna()
print(round(evalMetricData.median(), 4))
print(round(evalMetricData.max(), 4))
