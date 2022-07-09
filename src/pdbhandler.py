import os
from math import fmod
import numpy as np
import Bio.PDB
from Bio.PDB import is_aa
from Bio.PDB.PDBParser import PDBParser
import warnings
from rdkit import Chem
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import rdmolops
import requests
import time
import json
from sklearn import mixture
from src.silhouette import SilhouetteAnalysis
import pandas as pd
import random
import open3d as o3d
from sklearn.preprocessing import MinMaxScaler
import subprocess
import traceback
from Bio import SeqIO
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

class PDBHandler:

    def __init__(self, pdb_id=''):
        self.parser = PDBParser(PERMISSIVE=1)
        self.structure_id = pdb_id
        self.chain_angles = defaultdict()
        self.uniprot_accession_number = ''
        self.residue_selection = []
        self.domains = []
        self.points = []
        self.residue_distances = []
        self.root_disk = ''
        self.verbose = True
        self.residue_alternative_names = ['CYX', 'CYM', 'HID', 'HIE', 'HIP', 'ARN', 'ASH', 'GLH', 'LYN', 'TYM']
        self.camera_params_path = os.path.join('draw','camera-params')
        warnings.filterwarnings("ignore")  # supress PDBConstructionWarning for chains with missing residues

    def get_protein_chain_ids(self, pdb_path):
        pdb_path = self.handle_alphafold_pdbid(pdb_path)
        chains = []
        try:
            structure = self.parser.get_structure(self.structure_id, pdb_path)
            for model in structure:
                for chain in model:
                    if (self.is_protein(chain)):
                        chains.append(chain._id)
                break
        except:
            if (self.verbose is True):
                print(traceback.format_exc())
        return chains

    # Read PDB and calculate Phi and Psi angles
    def fetch_angles(self, pdb_path, chain_id=''):
        path_parts = pdb_path.split(os.path.sep)
        [self.structure_id, _] = path_parts[-1].split('.')

        try:
            structure = self.parser.get_structure(self.structure_id, pdb_path)  # making the structure of the .pdb file
            for model in structure:
                for chain in model:
                    if (chain_id != '' and chain_id != chain._id):
                        continue
                    if (chain_id == chain._id or self.is_protein(chain)):
                        self.chain_angles[chain._id] = defaultdict()
                        phi_angles = []
                        psi_angles = []
                        phi_psi_pairs = []

                        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                        poly.get_phi_psi_list()  # taking the torsion angles of the back bone atoms
                        angles = []

                        for residue in poly:
                            if (len(self.residue_selection) > 0 and residue.get_id()[1] not in self.residue_selection):
                                continue
                            phi = residue.xtra['PHI']
                            psi = residue.xtra['PSI']
                            angles.append((phi, psi))
                            if (phi is not None):
                                phi = self.process_angle(phi)
                                phi_angles.append(phi)

                            if (psi is not None):
                                psi = self.process_angle(psi)
                                psi_angles.append(psi)

                            if (phi is not None and psi is not None):
                                phi_psi_pairs.append([phi, psi])

                        if(len(self.residue_selection) > 0):
                            assert len(self.residue_selection) >= len(angles), ' '.join([pdb_path, chain_id])

                        self.chain_angles[chain._id]['angles'] = angles
                        self.chain_angles[chain._id]['phi'] = phi_angles
                        self.chain_angles[chain._id]['psi'] = psi_angles
                        self.chain_angles[chain._id]['phiPsiPairs'] = phi_psi_pairs

                # single model process only
                break
        except:
            if (self.verbose is True):
                print(traceback.format_exc())

    def handle_alphafold_pdbid(self, pdb_path):
        # AlphaFold PDBs support
        if ('AF-' in self.structure_id):
            parts = self.structure_id.split('-')
            if (len(parts) > 4):
                self.structure_id = '_'.join(['-'.join(parts[:-1]), parts[-1]])
            path_parts = pdb_path.split(os.path.sep)
            structure_path =  ''.join([self.structure_id, '.pdb'])
            pdb_path = ''.join([os.path.sep.join(path_parts[:-1]), os.path.sep, structure_path])
        return pdb_path

    def get_residue_range(self, pdb_path, chain_id):
        # Get all residue positions covered by the specified PDB chain
        full_range = []
        unusual_start_numbering = 0
        path_parts = pdb_path.split(os.path.sep)
        [self.structure_id, _] = path_parts[-1].split('.')
        try:
            #AlphaFold PDBs support
            pdb_path = self.handle_alphafold_pdbid(pdb_path)

            structure = self.parser.get_structure(self.structure_id, pdb_path)  # making the structure of the .pdb file
            for model in structure:
                for chain in model:
                    if (chain_id != chain._id):
                        continue
                    if (chain_id == chain._id):
                        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                        for residue in poly:
                            if ((is_aa(residue) is False and residue.resname not in self.residue_alternative_names) or residue.get_id()[0] != ' '):
                                continue
                            residue_id = residue.get_id()[1]
                            # https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
                            if (residue_id > 0):
                                full_range.append(residue_id)
                            else:
                                unusual_start_numbering += 1
                break
            full_range.sort()
        except:
            if (self.verbose is True):
                print(traceback.format_exc())

        return {'fullRange' : full_range, 'discardedStart' : unusual_start_numbering}

    def get_residues_names(self, pdb_path, chain_id):
        # Get the protein sequence of the PDB file
        sequence = []
        path_parts = pdb_path.split(os.path.sep)
        [self.structure_id, _] = path_parts[-1].split('.')
        try:
            #AlphaFold PDBs support
            pdb_path = self.handle_alphafold_pdbid(pdb_path)

            structure = self.parser.get_structure(self.structure_id, pdb_path)  # making the structure of the .pdb file
            for model in structure:
                for chain in model:
                    if (chain_id != chain._id):
                        continue
                    if (chain_id == chain._id):
                        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                        for residue in poly:
                            if ((is_aa(residue) is False and residue.resname not in self.residue_alternative_names) or residue.get_id()[0] != ' '):
                                continue
                            residue_id = residue.get_id()[1]
                            # https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
                            if (residue_id > 0):
                                sequence.append(residue.resname)
                break
        except:
            if (self.verbose is True):
                print(traceback.format_exc())

        return sequence

    # Determine if a chain of in the current PDB file is a peptidic one
    def is_protein(self, chain):
        for residue in chain.get_residues():
            if(residue.full_id[3][0] == ' '): #ignore HETATM entries
                return is_aa(residue) or residue.resname in self.residue_alternative_names

    # Convert angles from radiants to degrees
    def process_angle(self, angle, all_positive=False):
        angle = np.degrees(angle)

        if (all_positive is True):
            angle = fmod(angle, 360)
            if (angle < 0):
                angle += 360

        return angle

    # Retrieve all peptidic chains in the PDB file located in the given path
    def fetch_pdb_peptidic_chains(self, pdb_path):
        # AlphaFold PDBs support
        pdb_path = self.handle_alphafold_pdbid(pdb_path)

        molecule = Chem.MolFromPDBFile(pdb_path)
        if (molecule is not None):
            chains = rdmolops.SplitMolByPDBChainId(molecule)
            chain_ids = []
            for chain_id in chains:
                for atom in chains[chain_id].GetAtoms():
                    # check if it is protein
                    if (atom.GetPDBResidueInfo().GetResidueName().strip() not in ['A', "G", "C", "U"]):
                        chain_ids.append(chain_id)
                    break
            return {filtered_chain_id: chains[filtered_chain_id] for filtered_chain_id in chain_ids}
        else:
            if (self.verbose is True):
                print(pdb_path, ' - validation failed.')
            return False

    # Filter a molecule according to a residue selection
    def filter_molecule(self, molecule):
        molecule = Chem.EditableMol(molecule)
        atom_ids = []
        for atom in molecule.GetMol().GetAtoms():
            if (atom.GetPDBResidueInfo().GetResidueNumber() not in self.residue_selection):
                atom_ids.append(atom.GetIdx())
        atom_ids.reverse()
        for atom_id in atom_ids:
            molecule.RemoveAtom(atom_id)
        return molecule.GetMol()

    # Calculate a molecular fingerprint
    def calculate_fingerprint(self, molecule):
        return FingerprintMols.FingerprintMol(molecule)

    # Parse a PDB ID from a file path
    def parse_pdbid(self, pdb_path):
        parts = pdb_path.split(os.path.sep)
        return parts[-1].split('.')[0]

    def get_uniprot_accession_by_alphafold_pdbid(self):
        # Support for AlphaFold PDBs
        parts = self.structure_id.split('-')
        if(len(parts) > 1):
            if(parts[0] == 'AF'):
                self.uniprot_accession_number = parts[1]

    # Retrieve a UniProt accession number connected with a PDB ID
    def get_uniprot_accession_number(self, chain_id, pdb_id='', pdb_path=''):
        if (pdb_id != ''):
            self.structure_id = pdb_id
        self.uniprot_accession_number = ''
        uniprot_map_path = os.path.join(self.root_disk, 'uniprotmap')
        if (os.path.exists(uniprot_map_path) is False):
            os.makedirs(uniprot_map_path)
        # Check if there was a previously cached retrieval of this identifier
        data_path = os.path.join(uniprot_map_path, ''.join([self.structure_id, '.', chain_id, '.log']))
        if(self.uniprot_accession_number == ''):
            if (os.path.exists(data_path)):
                with open(data_path) as dataFile:
                    self.uniprot_accession_number = dataFile.readline().strip()
            else:
                self.get_uniprot_accession_by_alphafold_pdbid()
                if (self.uniprot_accession_number == ''):
                    # Extract the UniProt number from the PDB file
                    if (pdb_id != '' and pdb_path != ''):
                        self.uniprot_accession_number = self.get_uniprot_accession_by_seqres(pdb_path, pdb_id, chain_id)

                    # Search in the local mapping resources
                    if (self.uniprot_accession_number == ''):
                        self.uniprot_accession_number = self.read_uniprot_accession_number(chain_id)

                    # Query RCSB GraphQL service for the UniProt number
                    if (self.uniprot_accession_number == ''):
                        time.sleep(random.randint(10, 25))
                        request_url = ''.join(['https://data.rcsb.org/graphql'])
                        query = '''{
                                  polymer_entity_instances(instance_ids: ["PDB_CHAIN"])
                                  {
                                    polymer_entity
                                    {
                                      rcsb_polymer_entity_container_identifiers
                                      {
                                        uniprot_ids
                                      }
                                    }
                                  }
                                }'''.replace('PDB_CHAIN', ''.join([self.structure_id, '.', chain_id]))
                        if (self.verbose is True):
                            print(query)
                        response = requests.post(request_url, json={'query': query})
                        if (response.status_code == 200):
                            json_data = json.loads(response.text)
                            if (json_data['data'] is not None
                                    and len(json_data['data']) > 0
                                    and len(json_data['data']['polymer_entity_instances']) > 0
                                    and json_data['data']['polymer_entity_instances'][0]['polymer_entity'][
                                        "rcsb_polymer_entity_container_identifiers"]['uniprot_ids'] is not None
                                    and len(json_data['data']['polymer_entity_instances'][0]['polymer_entity'][
                                                "rcsb_polymer_entity_container_identifiers"]['uniprot_ids']) > 0):
                                self.uniprot_accession_number = \
                                    json_data['data']['polymer_entity_instances'][0]['polymer_entity'][
                                        "rcsb_polymer_entity_container_identifiers"]['uniprot_ids'][0]
                            else:
                                self.uniprot_accession_number = ''
                        if (self.verbose is True):
                            print(response.text)
        if (self.uniprot_accession_number != ''):
            with open(data_path, 'w') as outputFile:
                outputFile.write(self.uniprot_accession_number)
        return self.uniprot_accession_number

    # Seek a UniProt-PDB mapping in ID mapping resource
    def read_uniprot_accession_number(self, chain_id):
        uniprot_dataset_path = os.path.join(self.root_disk, 'idmapping_selected.tab.gz')
        pdb_id = ''.join([self.structure_id, ':', chain_id, ';'])

        uniprot_accession_number = ''

        try:
            output = subprocess.run(['zgrep', '-i', pdb_id, '-m', '1', uniprot_dataset_path], capture_output=True,
                                    timeout=30)
            parts = output.stdout.decode("utf-8").split('\t')

            # 'parts' includes the following information:
            # 1. UniProtKB-AC
            # 2. UniProtKB-ID
            # 3. GeneID (EntrezGene)
            # 4. RefSeq
            # 5. GI
            # 6. PDB
            # 7. GO
            # 8. UniRef100
            # 9. UniRef90
            # 10. UniRef50
            # 11. UniParc
            # 12. PIR
            # 13. NCBI-taxon
            # 14. MIM
            # 15. UniGene
            # 16. PubMed
            # 17. EMBL
            # 18. EMBL-CDS
            # 19. Ensembl
            # 20. Ensembl_TRS
            # 21. Ensembl_PRO
            # 22. Additional PubMed

            uniprot_accession_number = parts[0]
        except:
            if (self.verbose is True):
                print(''.join([pdb_id, ' not found in idmapping']))
                print(traceback.format_exc())
        return uniprot_accession_number

    # Get associated GO information for a PDB chain
    def get_gene_ontology_information(self, chain_id):
        time.sleep(random.randint(5, 15))
        request_url = ''.join(['https://data.rcsb.org/graphql'])
        query = '''{
				  polymer_entity_instances(instance_ids: ["PDB_CHAIN"])
				  {
					polymer_entity
					{
					  rcsb_polymer_entity_annotation {
						annotation_id
					  }
					}
				  }
				}'''.replace('PDB_CHAIN', ''.join([self.structure_id, '.', chain_id]))

        print(query)
        response = requests.post(request_url, json={'query': query})
        if (response.status_code == 200):
            json_data = json.loads(response.text)
            if (json_data['data'] is not None
                    and len(json_data['data']) > 0
                    and len(json_data['data']['polymer_entity_instances']) > 0
                    and json_data['data']['polymer_entity_instances'][0]['polymer_entity'][
                        "rcsb_polymer_entity_container_identifiers"]['uniprot_ids'] is not None
                    and len(json_data['data']['polymer_entity_instances'][0]['polymer_entity'][
                                "rcsb_polymer_entity_container_identifiers"]['uniprot_ids']) > 0):
                self.uniprot_accession_number = json_data['data']['polymer_entity_instances'][0]['polymer_entity'][
                    "rcsb_polymer_entity_container_identifiers"]['uniprot_ids'][0]
            else:
                self.uniprot_accession_number = ''
        else:
            if (self.verbose is True):
                print(response.text)

    # Extract metadata from a PDB file for a specified PDB chain
    def get_pdb_info(self, pdb_path, chain_id):
        path_parts = pdb_path.split(os.path.sep)
        [self.structure_id, _] = path_parts[-1].split('.')

        organisms = []
        gene_names = []
        residues = 0.0
        structure = self.parser.get_structure(self.structure_id, pdb_path)  # making the structure of the .pdb file
        resolution = structure.header["resolution"]

        # Count the residues of the chain
        for model in structure:
            for chain_index, chain in enumerate(model):
                if (chain._id != chain_id):
                    continue
                for r in chain.get_residues():
                    if r.id[0] == ' ':
                        residues += 1

        # Collect molecule IDs of the chains in the PDB file
        chain_molecule_ids = []
        if ('compound' in structure.header):
            for molecule_id in structure.header['compound']:
                if ('chain' in structure.header['compound'][molecule_id]):
                    chain_ids = structure.header['compound'][molecule_id]['chain'].split(',')
                    for molecule_chain_id in chain_ids:
                        if (molecule_chain_id.strip() == chain_id.lower()):
                            chain_molecule_ids.append(molecule_id)

        # Get the associated gene and organism of the chain
        if (len(chain_molecule_ids) > 0):
            if ('source' in structure.header):
                for chain_molecule_id in chain_molecule_ids:
                    if (chain_molecule_id in structure.header['source']):
                        if ('organism_scientific' in structure.header['source'][chain_molecule_id]):
                            organisms.append(structure.header['source'][chain_molecule_id]['organism_scientific'].strip())
                        if ('gene' in structure.header['source'][chain_molecule_id]):
                            gene_names.append(structure.header['source'][chain_molecule_id]['gene'].strip())

        return residues, resolution if resolution is not None else -1.0, ','.join(organisms) if len(
            organisms) > 0 else '-', \
               ','.join(gene_names) if len(gene_names) > 0 else '-'

    # Extract SEQRES entry of the PDB file located in the given path
    def get_seqres_entry(self, pdb_path):
        seqres = {record.id: record.seq for record in SeqIO.parse(pdb_path, 'pdb-seqres')}
        return seqres

    # Extract UniProt accession number from a SEQRES entry for a specified PDB chain
    def get_uniprot_accession_by_seqres(self, pdb_path, pdb_id, chain_id):
        self.uniprot_accession_number = ''
        for record in SeqIO.parse(pdb_path, 'pdb-seqres'):
            if (record.id == ':'.join([pdb_id, chain_id])):
                if (len(record.dbxrefs) > 0):
                    for db_reference in record.dbxrefs:
                        if ('UNP:' in db_reference):
                            self.uniprot_accession_number = db_reference.split(':')[1].strip()
                            break
            if (self.uniprot_accession_number != ''):
                break
        return self.uniprot_accession_number

    # Retrieve information on domains for a UniProt accession number
    def get_domain_information(self):
        data_path = os.path.join(self.root_disk, 'domains', ''.join([self.uniprot_accession_number, '.json']))
        # Check if there was a previously cached retrieval of this information
        if (os.path.exists(data_path)):
            with open(data_path) as jsonFile:
                parsed = json.load(jsonFile)
        else:
            # Attempt to retrieve domain information from UniProt
            while True:
                time.sleep(random.randint(5, 15))
                request_url = ''.join(['https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=',
                                      self.uniprot_accession_number, '&categories=DOMAINS_AND_SITES&types=DOMAIN'])
                print(request_url)
                response = requests.get(request_url)
                if (response.status_code == 200 and '<!doctype html>' not in response.text):
                    parsed = json.loads(response.text)
                    with open(data_path, 'w') as outputFile:
                        outputFile.write(response.text)
                    break
                else:
                    if (self.verbose is True):
                        print(response.text)
        self.domains = []
        if (len(parsed) > 0):
            for domain in parsed[0]['features']:
                self.domains.append([domain['description'], domain['begin'], domain['end']])

    # Remove spatial outliers by Interquartile Range (IQR)
    def remove_spatial_outliers(self, residues):
        # IQR
        data_frame = pd.DataFrame(residues, columns=['residue'])
        first_quantile = data_frame['residue'].quantile(0.25)
        third_quantile = data_frame['residue'].quantile(0.75)
        iqr = third_quantile - first_quantile
        data_frame = data_frame[
            (data_frame['residue'] >= first_quantile - 1.5 * iqr) & (data_frame['residue'] <= third_quantile + 1.5 * iqr)]
        interval = data_frame.max() - data_frame.min()

        return interval, data_frame

    # Remove spatial outliers and identify contiguous areas in the segment
    def process_segment(self, segment_residues, available_residues):
        parts = defaultdict()
        existing = []

        for residue_id in segment_residues:
            if residue_id in available_residues:
                existing.append(residue_id)

        if (len(existing) > 0):
            interval, data_frame = self.remove_spatial_outliers(existing)

            if (interval['residue'] > 100):
                # Gaussian Mixture Models clustering of the residues in the segment
                x = np.array(existing)
                x = x.reshape(-1, 1)
                max_score, max_score_number = SilhouetteAnalysis.perform_analysis(x)
                gmm = mixture.GaussianMixture(n_components=max_score_number, init_params='kmeans', random_state=42)
                gmm.fit(x)

                labels = gmm.predict(x)
                for index, label in enumerate(labels):
                    if (label not in parts):
                        parts[label] = []
                    parts[label].append(existing[index])
            else:
                parts[0] = data_frame['residue'].to_list()

            for label in parts:
                interval, data_frame = self.remove_spatial_outliers(parts[label])
                parts[label] = data_frame['residue'].to_list()

        return parts

    # Euclidean distance of alpha carbons in two residues
    def calculate_residue_distance(self, residue_one, residue_two):
        if ('CA' not in residue_one or 'CA' not in residue_two):
            return False
        diff_vector = residue_one['CA'].coord - residue_two['CA'].coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

    # Read PDB and calculate distances
    def calculate_residue_distances(self, pdb_path, chain_id=''):
        self.residue_distances = []
        path_parts = pdb_path.split(os.path.sep)
        [self.structure_id, _] = path_parts[-1].split('.')
        try:
            structure = self.parser.get_structure(self.structure_id, pdb_path)  # making the structure of the .pdb file
            for model in structure:
                for chain in model:
                    if (chain_id != '' and chain_id != chain._id):
                        continue
                    if (chain_id == chain._id or self.is_protein(chain)):
                        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                        residues = [residue for residue in poly]
                        for residue_one_index, residue_one in enumerate(residues):
                            for residue_two_index in range(residue_one_index + 1, len(residues)):
                                residue_two = residues[residue_two_index]
                                if (len(self.residue_selection) > 0 and
                                    residue_one.get_id()[1] not in self.residue_selection):
                                    continue
                                if (len(self.residue_selection) > 0 and
                                    residue_two.get_id()[1] not in self.residue_selection):
                                    continue
                                if (residue_one.get_id()[1] == residue_two.get_id()[1]):
                                    continue
                                distance = self.calculate_residue_distance(residue_one, residue_two)
                                if (distance is not False):
                                    self.residue_distances.append(distance)
                break
            self.residue_distances = np.array(self.residue_distances, dtype='float')
        except:
            if (self.verbose is True):
                print(traceback.format_exc())
        return self.residue_distances

    # Construct and store a point cloud Open3D visualization
    def view_point_cloud(self, points, output_path):
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(points)
        pcd.colors = o3d.utility.Vector3dVector(sns.color_palette('colorblind', len(points)))
        self.save_visualization(pcd, output_path, 'pcd')

    # Construct and store a set of lines between points in the clous (Open3D visualization)
    def view_point_cloud_distances(self, points, output_path):
        line_set = o3d.geometry.LineSet()
        line_set.points = o3d.utility.Vector3dVector(points)
        lines = []
        for i in range(len(points)):
            for j in range(i, len(points)):
                lines.append([i,j])
        line_set.lines = o3d.utility.Vector2iVector(lines)
        line_set.colors = o3d.utility.Vector3dVector(sns.color_palette('colorblind', len(lines)))
        self.save_visualization(line_set, output_path, 'lineset')

    # Save an Open3D visualization
    def save_visualization(self, generated_geometry, output_path, geometry_type=''):
        vis = o3d.visualization.Visualizer()
        vis.create_window()
        if(geometry_type == 'pcd'):
            vis.get_render_option().point_size = 5.0
        elif(geometry_type == 'alpha'):
            vis.get_render_option().mesh_show_back_face = True
            vis.get_render_option().mesh_show_wireframe = True

        vis.add_geometry(generated_geometry)

        camera_params_path = os.path.join(self.camera_params_path, ''.join([ self.structure_id, '-', geometry_type,'-camera-params.json']))
        if(os.path.exists(camera_params_path)):
            # Read camera params (press 'p' first for saving a json with camera details)
            ctr = vis.get_view_control()
            parameters = o3d.io.read_pinhole_camera_parameters(camera_params_path)
            ctr.convert_from_pinhole_camera_parameters(parameters)
        else:
            vis.run()

        # Capture image
        vis.update_geometry(generated_geometry)
        vis.poll_events()
        vis.update_renderer()
        time.sleep(1)
        image = vis.capture_screen_float_buffer(True)
        image = np.asarray(image)

        #https://gist.github.com/zhou13/b4ee8e815aee83e88df5b865896aaf5a - imshow without margins
        plt.gca().set_axis_off()
        plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        plt.margins(0, 0)
        plt.gca().xaxis.set_major_locator(plt.NullLocator())
        plt.gca().yaxis.set_major_locator(plt.NullLocator())

        # Remove margins from output
        output_margins_path = os.path.join(self.camera_params_path, ''.join([self.structure_id, '-img-trunc-margins.json']))
        configuration = False
        if(os.path.exists(output_margins_path)):
            with open(output_margins_path, 'r') as margins_file:
                configuration = margins_file.read()
            configuration = json.loads(configuration)
        if (configuration is not False):
            plt.imshow(image[int(image.shape[0] * configuration[geometry_type]['height-left']):int(image.shape[0] * configuration[geometry_type]['height-right']),
                       int(image.shape[1] * configuration[geometry_type]['width-left']):int(image.shape[1] * configuration[geometry_type]['width-right'])])
        else:
            plt.imshow(image)

        plt.savefig(output_path, format='png', dpi=300)
        plt.savefig(''.join([output_path[:-4], '.eps']), format='eps', dpi=300)
        plt.close()
        plt.cla()
        plt.clf()
        vis.destroy_window()

    # Construct and store an alpha shape Open3D visualization
    def view_alpha_shape(self, points, output_path):
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(points)
        alpha = 0.085
        print(f"alpha={alpha:.3f}")
        mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
        mesh.compute_vertex_normals()
        # o3d.visualization.draw_geometries([mesh], mesh_show_back_face=True, mesh_show_wireframe=True)
        self.save_visualization(mesh, output_path, 'alpha')

    # Create a point cloud based on given coordinates
    def load_point_cloud(self, points):
        pcd = o3d.geometry.PointCloud()
        try:
            pcd.points = o3d.utility.Vector3dVector(points)
        except:
            print(traceback.format_exc())
            pcd = None
        return pcd

    # Construct a triangle mesh based on a given point cloud
    def get_mesh_triangles(self, points):
        alpha = 0.085
        pcd = self.load_point_cloud(points)
        result = None
        if (pcd is not None):
            try:
                mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_alpha_shape(pcd, alpha)
                result = mesh.triangles
            except:
                result = None
        return result

    # Load coordinates from a PDB file
    def load_points(self, pdb_path, chain_id='', normalized=True, alpha_carbons_only=False):
        path_parts = pdb_path.split(os.path.sep)
        [self.structure_id, _] = path_parts[-1].split('.')
        self.points = []
        try:
            structure = self.parser.get_structure(self.structure_id, pdb_path)  # making the structure of the .pdb file
            for model in structure:
                for chain in model:
                    if (chain_id != '' and chain_id != chain._id):
                        continue
                    if (chain_id == chain._id or self.is_protein(chain)):
                        poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                        for residue in poly:
                            if (len(self.residue_selection) > 0 and residue.get_id()[1] not in self.residue_selection):
                                continue
                            for atom in residue:
                                if(alpha_carbons_only is True and 'CA' != atom.name):
                                    continue
                                self.points.append(atom.get_coord())
                break
        except:
            if (self.verbose is True):
                print(traceback.format_exc())

        self.points = np.array(self.points)
        if (normalized is True and len(self.points) > 0):
            scaler = MinMaxScaler()
            self.points = scaler.fit_transform(self.points)

    # Get an ascending order number for a chain in a PDB file at the given path
    def get_chain_index(self, full_pdb_file_path, chain_id):
        chain_index = -1
        chains = self.fetch_pdb_peptidic_chains(full_pdb_file_path)
        if (chains is not False):
            indexing = [index for index, key in enumerate(list(chains.keys())) if key == chain_id]
            if (len(indexing) > 0):
                chain_index = indexing[0]
        return chain_index

    # Discard potential illegal characters for the host file system
    def sanitize_domain_name(self, domain_name):
        return domain_name.replace(os.path.sep, '-').replace(' ', '-').replace('_', '-').replace('?', '').replace(
            '!', '').replace('\\', '-').replace('@', '-').replace('*', '-').replace('%', '-').replace(
            '\'', '').replace('"', '').replace(';', '').replace(';', '').replace(',', '').replace(':', '').replace(
            '+', 'PLUS').replace('(', '').replace(')', '')

    def filter_pdb_by_whitelist(self, full_pdb_file_path, chain_id, residue_whitelist):
        with open(full_pdb_file_path) as pdbFile:
            with open("".join([full_pdb_file_path.replace(".pdb", ""), "c.pdb"]), 'w') as outputFile:
                for line in pdbFile.readlines():
                    if (line[0:4] == 'ATOM' or line[0:6] == 'HETATM'):
                        residue_number = int(line[22:26])
                        if(residue_number not in residue_whitelist and line[21] == chain_id):
                            continue
                    outputFile.write(line)
