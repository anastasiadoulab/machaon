import sys
sys.path.append('../..')
import os
from src.pdbhandler import PDBHandler


pdbhandler = PDBHandler()
pdbhandler.camera_params_path = 'camera-params/'

pdbhandler.root_disk = 'structural-data'
pdb_dataset_path = 'structural-data'

pdbhandler.load_points(os.path.sep.join([pdb_dataset_path, 'PDBs_vir', '6VXX.pdb']), 'A', normalized=False, alpha_carbons_only=True)
pdbhandler.structure_id = '6VXX_A'
pdbhandler.view_point_cloud(pdbhandler.points, 'output/point-cloud.png')


pdbhandler.load_points(os.path.sep.join([pdb_dataset_path, 'PDBs_vir', '6VXX.pdb']), 'A', normalized=False, alpha_carbons_only=True)
pdbhandler.structure_id = '6VXX_A'
pdbhandler.view_point_cloud_distances(pdbhandler.points, 'output/lines.png')


pdbhandler.load_points(os.path.sep.join([pdb_dataset_path, 'PDBs_vir', '6VXX.pdb']), 'A')
pdbhandler.structure_id = '6VXX_A'
pdbhandler.view_alpha_shape(pdbhandler.points, 'output/alpha-shape.png')
