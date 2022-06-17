import pandas as pd

results = pd.read_csv('6VXX_A-merged-enriched.csv', sep='\t')

results = results[results['pdbId'] != '6VXX']
results = results[results['geneId'] == '43740568']

pdb_chains = ['_'.join([row['pdbId'], row['chainId']]) for rowIndex, row in results.iterrows()]

with open('pdb_chains.log', 'w') as log_file:
    log_file.write('\n'.join(pdb_chains))

print('Total structures: ', len(results.index))

for metric in ['b-phipsi', 'w-rdist', 't-alpha']:
    print(metric, repr(results[metric].median()))

