import sys
sys.path.append('..')
import os
from src.structurealigner import StructureAligner
from src.pdbhandler import PDBHandler
import numpy as np
import pandas as pd
from src.machaon import Machaon
from src.configurationmanager import ConfigurationManager
import copy
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D
import functools
import operator
plt.rcParams.update({'font.size': 8})

def plot_alignments(mask, reference_residue_range, residue_range, labels, figure_filename):
    ax = plt.axes()
    plt.figure(figsize=(6, 2))
    ax.set_yticks([1])
    plt.tick_params(labelleft=False)
    color_palette = sns.color_palette('colorblind', 1)

    # Plot alignments
    plt.plot(range(1, len(mask) + 1), mask, color=color_palette[0], linewidth=0.5)

    # Plot unavailable residues
    unavailable_residues = [x for x in range(1, len(mask) + 1) if x not in reference_residue_range]
    plt.scatter(unavailable_residues, [0 for x in range(len(unavailable_residues))], color='black', marker='|')
    unavailable_residues = [x for x in range(1, len(mask) + 1) if x not in residue_range]
    plt.scatter(unavailable_residues, [-0.1 for x in range(len(unavailable_residues))], color='gray', marker='|')

    # Annotate plot
    plt.xlabel('Residues')
    plt.ylabel('2D Matches')

    # Create legends for the plot
    legend_elements = []
    legend_elements.append(plt.Rectangle((0, 0), 1, 1, label=labels[0], fc='black'))
    legend_elements.append(plt.Rectangle((0, 0), 1, 1, label=labels[1], fc='gray'))
    plt.legend(handles=legend_elements, facecolor='white', framealpha=0.5, title='Missing', loc='best', handlelength=1, handleheight=1.125)

    # Store plot
    plt.tight_layout()
    ax = plt.gca()
    ax.set_rasterized(True)
    plt.savefig(''.join([figure_filename, '.png']), format='png', dpi=300)
    plt.savefig(''.join([figure_filename, '.eps']), format='eps', dpi=300)
    plt.clf()
    plt.figure(figsize=plt.rcParams.get('figure.figsize'))

# Initializations
aligner = StructureAligner()
aligner.set_root_disk('structural-data')
aligner.pdb_dataset_path = 'structural-data/PDBs_vir'

pdbhandler = PDBHandler()
pdbhandler.root_disk = 'structural-data'
pdb_dataset_path = 'structural-data/PDBs_vir'

sequence_length = 1273

# Get residue ranges, sequences
residue_range_N = pdbhandler.get_residue_range(os.path.sep.join([pdb_dataset_path, '6VXX.pdb']), 'A')
sec_N, pri_N = aligner.get_pdb_sequences('6VXX', 'A', residue_range_N)

residue_range_O = pdbhandler.get_residue_range(os.path.sep.join([pdb_dataset_path, '7T9K.pdb']), 'A')
sec_O, pri_O = aligner.get_pdb_sequences('7T9K', 'A', residue_range_O)

residue_range_D = pdbhandler.get_residue_range(os.path.sep.join([pdb_dataset_path, '7V7Q.pdb']), 'A')
sec_D, pri_D = aligner.get_pdb_sequences('7V7Q', 'A', residue_range_D)

# Intersect residue positions
residue_range_intersection = []
residue_range_intersection.extend(residue_range_N['fullRange'])
residue_range_intersection.extend(residue_range_O['fullRange'])
residue_range_intersection.extend(residue_range_D['fullRange'])
residue_range_intersection = list(set(residue_range_N['fullRange']) & set(residue_range_O['fullRange']) & set(residue_range_D['fullRange']))
print('Intersected available residue positions: ', repr(residue_range_intersection))
print('Total residues: ', len(residue_range_intersection))

# Align secondary structures and produce zero-padded masks to protein sequence length
result_NO = aligner.align(sec_N, sec_O, 'OTHER', {'gapChars': ['{', ' ']})
result_NO_mask = aligner.get_alignment_mask(result_NO[1])
mask_start = min(residue_range_N['fullRange']) - 1
result_NO_mask = [0 if (index + 1) not in residue_range_N['fullRange'] else result_NO_mask[index - mask_start] for index in range(sequence_length)]

result_ND = aligner.align(sec_N, sec_D, 'OTHER', {'gapChars': ['{', ' ']})
result_ND_mask = aligner.get_alignment_mask(result_ND[1])
mask_start = min(residue_range_N['fullRange']) - 1
result_ND_mask = [0 if (index + 1) not in residue_range_N['fullRange'] else result_ND_mask[index - mask_start] for index in range(sequence_length)]

result_DO = aligner.align(sec_D, sec_O, 'OTHER', {'gapChars': ['{', ' ']})
result_DO_mask = aligner.get_alignment_mask(result_DO[1])
mask_start = min(residue_range_D['fullRange']) - 1
result_DO_mask = [0 if (index + 1) not in residue_range_D['fullRange'] else result_DO_mask[index - mask_start] for index in range(sequence_length)]

# Plot alignment results
plot_alignments(result_NO_mask, residue_range_N['fullRange'], residue_range_O['fullRange'], ['Native', 'Omicron'], 'alignment_NO')
plot_alignments(result_ND_mask, residue_range_N['fullRange'], residue_range_D['fullRange'], ['Native', 'Delta'], 'alignment_ND')
plot_alignments(result_DO_mask, residue_range_D['fullRange'], residue_range_O['fullRange'], ['Delta', 'Omicron'], 'alignment_DO')

# Create PDBs with the same residue positions
pdbhandler.filter_pdb_by_whitelist(os.path.sep.join([pdb_dataset_path, '6VXX.pdb']), 'A', residue_range_intersection)
pdbhandler.filter_pdb_by_whitelist(os.path.sep.join([pdb_dataset_path, '7T9K.pdb']), 'A', residue_range_intersection)
pdbhandler.filter_pdb_by_whitelist(os.path.sep.join([pdb_dataset_path, '7V7Q.pdb']), 'A', residue_range_intersection)

# Create configurations for each PDB
current_working_directory = os.getcwd()
os.chdir(os.path.join(current_working_directory,'../..'))
go_target_properties  = ['virion attachment to host cell', 'receptor-mediated virion attachment to host cell']
job_configs = []

config_manager = ConfigurationManager()
configuration = copy.deepcopy(config_manager._template_config)
configuration['rootDisk'] = "structural-data"
configuration['referenceChainID'] = 'A'
configuration['referenceSequenceLength'] = sequence_length
configuration['referenceGeneID'] = 43740568
configuration['excludedOrganisms'] = ['severe acute respiratory syndrome coronavirus2',
                                      'severe acute respiratory syndrome coronavirus 2',
                                      'covid-19 virus']
configuration['pdbDatasetPath'] = "PDBs_vir"
configuration['isReferenceViral'] = True
configuration['GOTargetProperties'] = ['virion attachment to host cell', 'receptor-mediated virion attachment to host cell']
configuration['GOProperty'] = 'biologicalProcess'

segments, known = config_manager.parse_segment_config()

for pdb_id in ['6VXXc', '7V7Qc', '7T9Kc']:
    job_configs.append(copy.deepcopy(configuration))
    job_configs[-1]['referencePDBID'] = pdb_id
    job_configs[-1]['overridePDBID'] = pdb_id.replace('c', '')
    structure_id = '.'.join([job_configs[-1]['overridePDBID'], job_configs[-1]['referenceChainID']])
    if(structure_id in known):
        job_configs[-1]['known'] = known[structure_id]
    job_configs[-1]['outputPath'] = ''.join(['../output_fix/', pdb_id, '_A_whole'])

#  include custom PDBs in the dataset and run comparisons

machaon_core = Machaon()
machaon_core.debugging = False
machaon_core.max_cores = 31
machaon_core.process_new_pdbs([job_configs[0]])
machaon_core.perform_comparisons(job_configs)
os.chdir(current_working_directory)