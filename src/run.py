import argparse
import sys
sys.path.append('..')
from src.configurationmanager import ConfigurationManager
from src.machaon import Machaon

# Console arguments for running Machaon
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Machaon v1.0 : Find common folds between proteins.')
    parser.add_argument('cores', type=int, help='The number of cores available for Machaon.')
    parser.add_argument('rootDisk', type=str, help='The directory for source data and cache.')
    parser.add_argument('pdbDatasetPath', type=str, help='The directory for the target dataset of PDB files (reference PDB included).')
    parser.add_argument('outputPath', type=str, help='The directory for the outputs of searches')
    parser.add_argument('referencePDBChain', type=str, help='Provide PDB ID and Chain ID separated by a dot: <PDB ID>.<CHAIN ID> e.g. ""6VXX.A".',)
    parser.add_argument('referenceSequenceLength', type=int, help='Provide protein sequence length of the reference protein.')
    parser.add_argument('isReferenceViral', type=int, help='1 if reference protein is a viral one else 0.', choices=[0,1])
    parser.add_argument('--segmentResidues', type=str, help='Provide residue positions (optional).', default='')
    parser.add_argument('--segmentResidueRanges', type=str, help='Provide residue ranges (optional).', default='')
    parser.add_argument('--referenceGeneID', type=str, help='Provide NCBI Entrez Gene ID to exclude all related structures in the final result (optional).', default='')
    parser.add_argument('--comparisonMode', type=str, help='Choose between "whole", "domain" or "segment" options. Default is "whole".', default='whole', choices=['whole', 'domain', 'segment'])
    parser.add_argument('--alignmentLevel', type=str, help='This is used for "segment" mode. Choose between "1D" (protein sequence), "2D" (protein secondary structure) "hydrophobicity", "mixed". Default is "mixed".', default='mixed', choices=['1D', '2D', 'hydrophobicity', 'mixed'])
    parser.add_argument('--excludedOrganisms', type=str, help='Exclude organisms (as shown in the metadata of PDB files) from the search. Separate names by commas. (optional)', default='')
    parser.add_argument('--excludedGeneNames', type=str, help='Exclude genes (as shown in the metadata of PDB files) from the search. Separate names by commas. (optional)', default='')
    parser.add_argument('--excludedPDBIDs', type=str, help='Exclude PDB IDs from the search. Separate IDs by commas.', default='')
    parser.add_argument('--GOProperty', type=str, help='Choose a property type for analysis: "biologicalProcess", "molecularFunction", "cellularComponent". (optional)', default='', choices=['biologicalProcess', 'molecularFunction', 'cellularComponent'])
    parser.add_argument('--GOTargetProperties', type=str, help='Choose properties for analysis. Separate names by commas. (optional)', default='')
    parser.add_argument('--GOSearch', type=str, help='Choose a term to be searched in all available GO Terms belonging to the results. (optional)', default='')
    parser.add_argument('--GOAlignmentLevel', type=str, help='Choose between "1D" (protein sequence), "2D" (protein secondary structure), "hydrophobicity", "mixed". Default is "2D". (optional)', default='2D', choices=['1D', '2D', 'hydrophobicity', 'mixed'])
    parser.add_argument('--goAnalysisOnly', help='Perform GO meta-analysis only (for already completed searches).', default=False, action='store_true')
    parser.add_argument('--noThirdPartyData', help='Do not use any local or online third-party data', default=False, action='store_true')
    parser.add_argument('--pdbValidation', help='Validate PDBs', default=True, action='store_true')
    parser.add_argument('--overridePDBID', type=str, help='Override the reference PDBID for Uniprot ID retrieval (for renamed reference PDB files, e.g. 6VXX_processed.pdb)', default='')

    args = vars(parser.parse_args())

    print('Machaon\'s arguments:\n', args)

    config_manager = ConfigurationManager()
    configuration = config_manager.parse_console_arguments(args)

    comparison = Machaon()
    comparison.max_cores = args['cores']
    comparison.perform_comparisons([configuration])