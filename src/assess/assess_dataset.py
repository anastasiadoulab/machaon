import os
from multiprocessing import Pool
from tqdm import tqdm
from enricher import Enricher
from pdbhandler import PDBHandler
import pandas as pd

rootDisk = 'structural-data'
pdbDatasetPath = os.path.sep.join([rootDisk, 'PDBs_vir'])
pdbFilenames = [filename for filename in os.listdir(pdbDatasetPath) if os.path.isfile(os.path.join(pdbDatasetPath, filename))]
results = []
enricher = Enricher()
enricher.pdb_dataset_path = 'PDBs_vir'
enricher.set_root_disk(rootDisk)

def getOrganisms(PDBfilename):
    info = []
    pdbhandler = PDBHandler()
    pdbhandler.root_disk = rootDisk
    chains = pdbhandler.fetch_pdb_peptidic_chains(os.path.join(pdbDatasetPath, PDBfilename))
    if(chains is not False):
        for chainId in chains:
            entry = PDBfilename.split('.')[0], chainId
            info.append([*enricher.access_pdb_info(entry)])
    return info

with Pool(20) as pool:
    for result in tqdm(pool.imap(getOrganisms, pdbFilenames), total=len(pdbFilenames)):
        if(len(result) > 0):
            results.extend(result)
    pool.close()
    pool.join()

dataframe = pd.DataFrame(results, columns=['pdb_id', 'residues', 'resolution', 'organisms', 'genes'])
dataframe.to_csv(os.path.sep.join([rootDisk, 'PDBs_vir_info.csv']), sep='\t')