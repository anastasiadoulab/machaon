import os
import pandas as pd

results_folder = 'output_pdb70'
queries_list = 'scope_140_targets.list'
chain_ids_list = 'scop.pdbid.chains.list'
output_folder = 'test_input'
os.makedirs(output_folder, exist_ok=True)

excluded_id = 'd2zjq51'
scop_queries = pd.read_csv(queries_list, sep='\t', header=None)
query_chains = pd.read_csv(chain_ids_list, sep='\t', names=['sid', 'cid'])
structure_ids = scop_queries[2].to_list()
chain_ids =  query_chains.set_index('sid').to_dict()['cid']
for structure_id in structure_ids:
    if structure_id in excluded_id:
        continue
    results_path = os.path.join(results_folder, structure_id.replace('_',''), 'candidates', ''.join([structure_id.replace('_',''), '_',
                                 chain_ids[structure_id.replace('_','')], '-merged-cluster.csv']))
    query_pairs = pd.read_csv(results_path, sep='\t')
    total_entries = len(query_pairs.index)
    query_pairs['structureId'] = query_pairs['structureId'].str.split(r'\_').str.get(0)
    query_pairs['queryid'] = [structure_id for x in reversed(range(total_entries))]
    query_pairs['rank'] = [float(x) for x in reversed(range(total_entries))]
    tupple_pairs = query_pairs[['structureId', 'queryid', 'rank']]
    tupple_pairs.to_csv(os.path.join(output_folder, structure_id), index=False, header=False, sep='\t')
