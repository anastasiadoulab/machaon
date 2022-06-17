import copy
import os
from collections import defaultdict
import yaml
import traceback

class ConfigurationManager:

    def __init__(self):
        self._config_file_path = '../config/config.yaml'
        self._segment_config_path = '../config/segments.yaml'

        #A job configuration template
        self._template_config = defaultdict()

        # The directory where the structures that will be compared are stored.
        # This will also be the caching location for the extracted features.
        self._template_config['rootDisk'] = ''

        # Choose the reference PDB IDs (1 search per reference)
        self._template_config['referencePDBID'] = ''

        # Override the reference PDBID for Uniprot ID retrieval (for renamed reference PDB files, e.g. 6VXX_processed.pdb)
        self._template_config['overridePDBID'] = ''

        #  Choose the chain of the reference PDB
        self._template_config['referenceChainID'] = ''

        #  Provide the gene id (Entrez) of the reference PDB
        self._template_config['referenceGeneID'] = ''

        #  Provide the protein sequence length of the reference protein
        self._template_config['referenceSequenceLength'] = 0

        #  Choose 'whole', 'domain' or 'segment'
        self._template_config['comparisonMode'] = 'whole'

        # Choose 'primary', 'secondary', 'mixed', 'hydrophobicity'. Default is 'mixed'. (Only from segment scans)
        self._template_config['alignmentLevel'] = 'mixed'

        # Relative path for PDB data folder
        self._template_config['pdbDatasetPath'] = ''

        # The location of the outputs (can be relative or full path)
        self._template_config['outputPath'] = ''

        # Filtering out structures originating from the same organism as the reference one
        self._template_config['excludedOrganisms'] = []

        # Filtering out structures originating from the same gene as the reference one
        self._template_config['excludedGeneNames'] = []

        # Exclude PDB IDs
        self._template_config['excludedPDBIDs'] = []

        # Meta-analysis skips the search in viral genome data for the reference,
        # if it is not a viral protein
        self._template_config['isReferenceViral'] = None

        # Choose a term to be searched in all available GO Terms belonging to the results e.g. 'ubiquit' (could be a stem of a word)
        self._template_config['GOSearch'] = ''

        # Choose a property type for analysis: 'biologicalProcess', 'molecularFunction', 'cellularComponent'
        self._template_config['GOProperty'] = ''

        # Choose properties for joint analysis
        self._template_config['GOTargetProperties'] = []

        # Choose target alignment level : ['primary', 'secondary', 'mixed', 'hydrophobicity'. Default is 'secondary']
        self._template_config['GOAlignmentLevel'] = 'secondary'

        # Do not use external local or online resources. PDB data only.
        self._template_config['noThirdPartyData'] = False

        # Do not use external local or online resources. PDB data only.
        self._template_config['isNonRedundantSet'] = False

        # Perform only GO Meta-analysis (for completed searches).
        self._template_config['GOAnalysisOnly'] = False

        # Segment presets ('segment' constrained mode only)
        self._template_config['segments'] = []

        # Denotes whether a segment is a characterized site or not
        self._template_config['known'] = []

        # Validation for PDB files (strict)
        self._template_config['pdbValidation'] = False

    def assign_default_root_folder(self, path):
        self._template_config['rootDisk'] = path

    def parse_console_arguments(self, console_arguments):
        configuration = copy.deepcopy(self._template_config)
        residues = []
        try:
            # Build configuration by each console argument
            for key in console_arguments:
                if(key == 'isReferenceViral'):
                    configuration[key] = console_arguments[key] == 1
                elif(key == 'referencePDBChain'):
                    # Split structure id to get pdb and chain ids
                    structure_id = console_arguments[key].strip().split('.')
                    if(len(structure_id) != 2):
                        raise ValueError(''.join(['wrong input for ', key]))
                    else:
                        configuration['referencePDBID'] = structure_id[0]
                        configuration['referenceChainID'] = structure_id[1]
                elif key == 'segmentResidueRanges':
                    # Collect residues in the specified ranges that belong to a segment
                    if(len(console_arguments[key]) > 0):
                        residues.extend(self.expand_residue_ranges(console_arguments[key]))
                elif key == 'segmentResidues':
                    # Collect residues that belong to a segment
                    if (len(console_arguments[key]) > 0):
                        residues.extend([int(x) for x in console_arguments[key].split(',')])
                elif key in ['excludedOrganisms', 'excludedGeneNames', 'excludedPDBIDs', 'GOTargetProperties']:
                    # Set exclusion criteria for entries in candidate set
                    if(len(console_arguments[key].strip()) < 1):
                        continue
                    configuration[key] = console_arguments[key].strip().split(',')
                else:
                    configuration[key] = console_arguments[key]
        except ValueError as error_object:
            print('There is an error in your argument:\n "', str(error_object), '". \nPlease run the script with the -h option only ',
                   'to review the valid arguments of the program.')
            print(traceback.format_exc())
            exit(0)
        # Set collected residue statements for the constrained reference
        if(len(residues) > 0):
            residues = list(set(residues))
            residues.sort()
            configuration['segments'] = {console_arguments['referencePDBChain'] : [residues]}
        else:
            print('Note: No residue positions are specified. Ignore this message for whole structure comparisons.\n')
        return configuration

    # Parse configuration file in YAML format
    def parse_config(self):
        with open(self._config_file_path) as config_file:
            contents = yaml.safe_load(config_file)
        collected = []
        single_chain_config = True
        single_reference_config = True
        segments, known = self.parse_segment_config()
        for config_index, configuration in enumerate(contents['configurations']):
            if('ignore' in configuration and configuration['ignore'] is True):
                continue
            # Create a copy of the default configuration template and update it with the current settings
            default_preset = copy.deepcopy(self._template_config)
            for key in configuration:
                default_preset[key] = configuration[key]
            if(isinstance(default_preset['referencePDBID'], list) is False):
                default_preset['referencePDBID'] = [default_preset['referencePDBID']]
            else:
                # If multiple references are set, create configurations for individual sessions
                single_reference_config = False
                if(isinstance(default_preset['referenceChainID'], list) is True):
                    assert len(default_preset['referenceChainID']) == len(default_preset['referencePDBID']), \
                        ''.join(['The number of PDB IDs and chain IDs in the configuration are not equal. [Job No. ', config_index, ' in yaml list]'])
                    single_chain_config = False
            for reference_index, reference_pdb_id in enumerate(default_preset['referencePDBID']):
                preset = copy.deepcopy(default_preset)
                preset['referencePDBID'] = reference_pdb_id
                if(single_reference_config is False):
                    os.makedirs(preset['outputPath'] , exist_ok=True)
                    preset['outputPath'] = os.path.join(default_preset['outputPath'], reference_pdb_id)
                if(single_chain_config is False):
                    preset['referenceChainID'] = default_preset['referenceChainID'][reference_index]
                self.verify_config(preset, config_index)
                structure_id = '.'.join([preset['referencePDBID'],preset['referenceChainID']])
                # Collect information about preset segments
                if(structure_id in segments):
                    preset['segments'] = segments[structure_id]
                    preset['known'] = known[structure_id]
                # include segments for the alternative PDBID
                if(len(preset['overridePDBID']) > 0):
                    structure_id = '.'.join([preset['overridePDBID'], preset['referenceChainID']])
                    if (structure_id in segments):
                        preset['segments'] = segments[structure_id]
                        preset['known'] = known[structure_id]
                collected.append(preset)
        return collected

    def parse_segment_config(self):
        with open(self._segment_config_path) as config_file:
            contents = yaml.safe_load(config_file)
        collected = defaultdict(list)
        known = defaultdict(list)
        for preset in contents['segments']:
            # Ignore this entry if specified
            if ('ignore' in preset and preset['ignore'] is True):
                continue
            # Collect participating residues from presets
            residues = []
            for key in preset:
                if(key == 'residueRanges'):
                    residues.extend(self.expand_residue_ranges(preset[key]))
                elif(key == 'residues'):
                    residues.extend([int(x) for x in preset[key]])
            if(len(residues) > 0):
                residues = list(set(residues))
                residues.sort()
                collected[preset['referencePDBChain']].append(residues)
                # Statement for a segment of known function or other characteristic
                if ('known' in preset and  preset['known'] is True):
                    known[preset['referencePDBChain']].append(residues)
        return collected, known

    # Convert residue ranges configuration to individual residue position statements
    def expand_residue_ranges(self, residue_ranges):
        residues = []
        try:
            residue_ranges = residue_ranges.split(',')
            for residue_range in residue_ranges:
                residue_range = residue_range.split('-')
                if (len(residue_range) != 2):
                    raise ValueError
                start = int(residue_range[0])
                end = int(residue_range[1])
                if(start > end):
                    raise ValueError
                residues.extend([x for x in range(start, end + 1)])
        except ValueError:
            print('There is an error in your segment preset: ', residue_range)
            exit(0)
        return residues

    # Verification for the presence of minimum required settings
    def verify_config(self, configuration, config_index):
        for key in configuration:
            condition = key in ['rootDisk', 'referencePDBID', 'referenceChainID', 'pdbDatasetPath', 'outputPath',] and len(configuration[key]) == 0
            condition = condition or (key == 'referenceSequenceLength' and configuration[key] < 1)
            condition = condition or (key == 'isReferenceViral' and configuration[key] is None)
            if(condition):
                print('A required field is missing in your job configuration: "', key, '" [Job No. ', config_index, ' in yaml list]')
                exit(0)