import pandas as pd

raw_reference_results = pd.read_csv('6VXX_A-merged-enriched_eval_full.csv', sep='\t')
preprocessed_reference_results = pd.read_csv('6VXXp_A-merged-enriched_eval_full.csv', sep='\t')

raw_reference_results['composite_id'] = raw_reference_results.apply(lambda row: '_'.join([row['pdbId'], row['chainId']]), axis=1)
preprocessed_reference_results['composite_id'] = preprocessed_reference_results.apply(lambda row: '_'.join([row['pdbId'], row['chainId']]), axis=1)

for key in ['uniprotId']:
    print(key)
    raw_ref = raw_reference_results[key].to_list()[1:]
    preproc_ref = preprocessed_reference_results[key].to_list()[1:]
    common = len([x for x in raw_ref if x in preproc_ref])
    raw_unique = [x for x in raw_ref if x not in preproc_ref]
    raw_unique_positions = [index + 1 for index, x in enumerate(raw_ref) if x not in preproc_ref]
    print('raw_unique_positions:', repr(raw_unique_positions))
    preproc_unique = [x for x in preproc_ref if x not in raw_ref]
    preproc_unique_positions = [index + 1 for index, x in enumerate(preproc_ref) if x not in raw_ref]
    print('preproc_unique_positions:', repr(preproc_unique_positions))
    print('common: ', repr(common))
    preproc_ref.extend(raw_ref)
    superset_length = len(list(set(preproc_ref)))
    print('superset length: ', repr(superset_length))
    print('percentage of superset: ', repr(common / superset_length))
    print('\n')
preprocessed_reference_results[preprocessed_reference_results['uniprotId'].isin(preproc_unique)].to_csv('preproc_unique.tsv', sep='\t')
raw_reference_results[raw_reference_results['uniprotId'].isin(raw_unique)].to_csv('raw_unique.tsv', sep='\t')

