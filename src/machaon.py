from src.segmentscanner import SegmentScanner
from src.domainscanner import DomainScanner
from src.candidatesampler import CandidateSampler
from src.metaanalysis import MetaAnalysis
from src.presenter import Presenter
from src.scanner import Scanner
import os
import pandas as pd

class Machaon:

    def __init__(self):
        self.debugging = False
        self.max_cores = 2
        self.check_external()
        self.meta_metrics = ['1D-identity', '2D-identity', '5UTR-identity', 'CDS-identity', '3UTR-identity', '3D-score',
                             'chemSim', 'molecularFunctionSim', 'cellularComponentSim', 'biologicalProcessSim']
        os.makedirs('logs', exist_ok=True)

    def get_scanner(self, mode, alignment_backend, alignment_level=None):
        # Choose the appropriate Scanner
        if(mode == 'segment'):
            scanner = SegmentScanner()
            scanner.selected_alignment_level = alignment_level
            scanner.alignment_backend = alignment_backend
        elif(mode == 'domain'):
            scanner = DomainScanner()
        else:
            scanner = Scanner()
        return scanner

    def check_external(self):
        # Perform checks on external executables
        if(os.path.exists("./mkdssp") is False):
            print("Warning: DSSP executable is not found.\n")
        if(os.path.exists("./TMalign") is False):
            print("Warning: TM-Align executable is not found.\n")
        if(os.path.exists("./stride") is False):
            print("Warning: STRIDE executable is not found.\n")

    def perform_comparisons(self, job_config=None):
        if(job_config is None):
            print("There is no assignment for Machaon. Please provide a configuration.")
        else:
            # For every reference in configuration
            for reference_index, configuration in enumerate(job_config):
                # Collect features
                pdb_id = configuration['referencePDBID']
                print(' '.join(['--- Reference PDB ID: ', pdb_id, ' chain: ', configuration['referenceChainID'],
                                ' mode:', configuration['comparisonMode']]))
                scanner = self.get_scanner(configuration['comparisonMode'], configuration['alignmentBackend'],
                                           configuration['alignmentLevel'])
                scanner.cores = self.max_cores
                scanner.debugging = self.debugging
                scanner.verbose = self.debugging
                scanner.reference_pdb_id = pdb_id
                scanner.override_pdb_id = configuration['overridePDBID']
                scanner.reference_chain_id = configuration['referenceChainID']
                scanner.root_disk = configuration['rootDisk']
                scanner.pdb_dataset_path = configuration['pdbDatasetPath']
                scanner.pdb_validation = configuration['pdbValidation']
                scanner.features_path = '_'.join(['DATA', configuration['pdbDatasetPath'],
                                                  configuration['comparisonMode']])
                scanner.features_file_list = ''.join([scanner.features_path, '.csv'])
                scanner.metrics_output_path = os.path.join(configuration['outputPath'], 'metrics')
                if(configuration['GOAnalysisOnly'] is False):
                    if(configuration['comparisonMode'] != 'segment'):
                        print("*** Extracting features")
                        scanner.collect_features()
                    else:
                        scanner.segments = configuration['segments']
                    # Compute metrics for each protein in the candidate set
                    scanner.set_candidates()
                    print("*** Computing metrics")
                    scanner.scan_candidates()

                # Select candidates
                print("*** Selecting candidates")
                candidatesampler = CandidateSampler()
                candidatesampler.cores = self.max_cores
                candidatesampler.debugging = self.debugging
                candidatesampler.root_disk = scanner.root_disk
                candidatesampler.metrics_path = scanner.metrics_output_path
                candidatesampler.output_path = os.path.join(configuration['outputPath'], 'candidates')
                candidatesampler.plots_output_path = os.path.join(configuration['outputPath'], 'plots')
                candidatesampler.pdb_dataset_path = scanner.pdb_dataset_path
                candidatesampler.reference_pdb_id = scanner.reference_pdb_id
                candidatesampler.reference_chain_id = scanner.reference_chain_id
                candidatesampler.reference_gene_id = configuration['referenceGeneID']
                candidatesampler.excluded_organisms = [x.lower() for x in configuration['excludedOrganisms']]
                candidatesampler.selected_organisms = [x.lower() for x in configuration['selectedOrganisms']]
                candidatesampler.excluded_gene_names = configuration['excludedGeneNames']
                candidatesampler.excluded_pdb_ids = configuration['excludedPDBIDs']
                candidatesampler.no_enrichment = configuration['noThirdPartyData']
                candidatesampler.override_pdb_id = configuration['overridePDBID']
                candidatesampler.vector_format_output = configuration['epsOutput']
                candidatesampler.tiff_format_output = configuration['tiffOutput']
                candidatesampler.max_result_size = 250 if configuration['isNonRedundantSet'] is True else 800
                if (configuration['comparisonMode'] != 'whole'):
                    candidatesampler.set_segments(configuration['comparisonMode'], configuration['segments'])
                else:
                    candidatesampler.reference_segments = ['']
                if (configuration['GOAnalysisOnly'] is False):
                    candidatesampler.select_candidates()

                # Meta-analysis: compute metrics for transcript, chemical levels
                if (configuration['noThirdPartyData'] is False):
                    print("*** Meta-analysis")
                    metaanalysis = MetaAnalysis()
                    metaanalysis.cores = min(4, self.max_cores) # RAM usage per core is supposed to be high
                    metaanalysis.root_disk = scanner.root_disk
                    metaanalysis.debugging = self.debugging
                    metaanalysis.pdb_dataset_path = scanner.pdb_dataset_path
                    metaanalysis.reference_pdbid = scanner.reference_pdb_id
                    metaanalysis.reference_chain_id = scanner.reference_chain_id
                    metaanalysis.known_areas =  configuration['known']
                    metaanalysis.reference_segments = candidatesampler.reference_segments if \
                                                      len(candidatesampler.reference_segments) > 0 else ['']
                    metaanalysis.is_reference_viral = configuration['isReferenceViral']
                    metaanalysis.viral_content = configuration['viralContentExists']
                    metaanalysis.alignment_backend = configuration['alignmentBackend']
                    metaanalysis.output_path = candidatesampler.output_path
                    metaanalysis.go_output_path = os.path.join(configuration['outputPath'], 'go')
                    metaanalysis.go_alignment_level = configuration['GOAlignmentLevel']
                    metaanalysis.override_pdb_id = configuration['overridePDBID']
                    metaanalysis.vector_format_output = configuration['epsOutput']
                    metaanalysis.tiff_format_output = configuration['tiffOutput']

                    # Create meta-analysis data
                    if (configuration['GOAnalysisOnly'] is False):
                        print("Evaluating results...")
                        metaanalysis.evaluate_results()
                        print("Aggregating GO information...")
                        metaanalysis.aggregate_go_information()
                        print("Counting the 2D matches...")
                        for reference_segment in candidatesampler.reference_segments:
                            if (reference_segment != ''):
                                print(''.join(['Analysis for: ', reference_segment]))
                            naming_extension = reference_segment
                            if ('_' in naming_extension):
                                naming_extension = ''.join([naming_extension, '-metrics'])
                            entry = pdb_id, scanner.reference_chain_id, \
                                    naming_extension,candidatesampler.plots_output_path
                            metaanalysis.find_batch_structure_matches(entry)

                # Presenting the results
                if (configuration['GOAnalysisOnly'] is False):
                    print("*** Presenting results")
                    for reference_segment in candidatesampler.reference_segments:
                        naming_extension = reference_segment
                        if ('_' in naming_extension):
                            naming_extension = ''.join([naming_extension, '-metrics'])
                        presenter = Presenter()
                        presenter.root_disk = scanner.root_disk
                        presenter.verbose = self.debugging
                        presenter.debugging = self.debugging
                        presenter.cores = self.max_cores
                        presenter.comparison_mode = configuration['comparisonMode']
                        presenter.vector_format_output = configuration['epsOutput']
                        presenter.tiff_format_output = configuration['tiffOutput']
                        presenter.data_path = candidatesampler.output_path
                        presenter.plots_output_path = candidatesampler.plots_output_path
                        # create report for final cluster
                        print('Creating report for the final cluster...')
                        presenter.create_report([pdb_id, scanner.reference_chain_id, naming_extension],
                                                reduced_data=True)
                        if (configuration['noThirdPartyData'] is False):
                            # Create report for top 100
                            print('Creating report for top 100...')
                            presenter.create_report([pdb_id, scanner.reference_chain_id, naming_extension])
                            # Create report for top 100 human entries if exist
                            human_report_path =  ''.join([configuration['outputPath'], os.path.sep, 'candidates',
                                                          os.path.sep, pdb_id, '_', configuration['referenceChainID'],
                                                          naming_extension, '-merged-h-enriched.csv'])
                            if os.path.exists(human_report_path) is True:
                                print('Creating report for top 100 of human proteins...')
                                presenter.create_report([pdb_id, scanner.reference_chain_id, naming_extension],
                                                        file_name=human_report_path,  reduced_data=True)

                            # Create various data visualizations
                            print('Creating taxonomy trees..')
                            presenter.display_taxonomy([pdb_id, scanner.reference_chain_id, naming_extension])
                            print('Creating word clouds..')
                            presenter.create_word_cloud([pdb_id, scanner.reference_chain_id, naming_extension])
                            print('Creating polar plots..')
                            for meta_metric in self.meta_metrics:
                                presenter.plot_results_polar([pdb_id, scanner.reference_chain_id, naming_extension,
                                                              meta_metric])
                            print('Creating 2D coverage plots..')
                            presenter.plot_2D_coverage([pdb_id, scanner.reference_chain_id, naming_extension])

                if (configuration['noThirdPartyData'] is False):
                    # Gene Ontology analysis
                    print("*** Performing Gene Ontology analysis.")
                    for reference_segment in metaanalysis.reference_segments:
                        if(reference_segment != ''):
                            print(''.join(['Analysis for: ', reference_segment]))
                        naming_extension = reference_segment
                        if ('_' in naming_extension):
                            naming_extension = ''.join([naming_extension, '-metrics'])
                        entry = pdb_id, scanner.reference_chain_id , naming_extension
                        metaanalysis.get_available_go_information(entry)
                        # If there is a GO search term or specified GO terms in the configurations for investigation
                        # perform GO meta-analysis
                        if((len(configuration['GOProperty']) > 0
                            and len(configuration['GOTargetProperties']) > 0)
                                or len(configuration['GOSearch']) > 0):
                            entry = pdb_id, scanner.reference_chain_id , naming_extension, \
                                    configuration['GOProperty'], configuration['GOTargetProperties'], \
                                    configuration['GOSearch']
                            metaanalysis.perform_go_meta_analysis(entry)

                print('Completed.\n\n')

    def process_new_pdbs(self, job_config=None, new_file_paths=[]):
        print("*** Processing new PDBs")
        if len(new_file_paths):
            print(repr(new_file_paths))
        if(job_config is None):
            print("There is no assignment for Machaon. Please provide a configuration.")
        else:
            # For every reference in configuration
            for reference_index, configuration in enumerate(job_config):
                if configuration['comparisonMode'] == 'segment':
                    continue
                # Collect features from new PDBs
                pdb_id = configuration['referencePDBID']
                print(' '.join(['--- Reference PDB ID: ', pdb_id, ' chain: ', configuration['referenceChainID'], ' mode:', configuration['comparisonMode']]))
                scanner = self.get_scanner(configuration['comparisonMode'], configuration['alignmentBackend'],
                                           configuration['alignmentLevel'])
                scanner.cores = self.max_cores
                scanner.debugging = self.debugging
                scanner.verbose = self.debugging
                scanner.reference_pdb_id = pdb_id
                scanner.reference_chain_id = configuration['referenceChainID']
                scanner.root_disk = configuration['rootDisk']
                scanner.pdb_dataset_path = configuration['pdbDatasetPath']
                scanner.pdb_validation = configuration['pdbValidation']
                scanner.features_path = '_'.join(['DATA', configuration['pdbDatasetPath'], configuration['comparisonMode']])
                scanner.features_file_list = ''.join([scanner.features_path, '.csv'])
                scanner.metrics_output_path = os.path.join(configuration['outputPath'], 'metrics')
                scanner.extend_feature_data = True
                scanner.collect_features(new_file_paths)

                print('Completed.\n\n')

    # Update an existing filelist or create a new instance of a modified filelist
    def manage_file_list(self, data_folder_path, new_files, mode):
        list_path = ''.join([data_folder_path, '.csv'])
        if mode != 'segment':
            file_exists = os.path.exists(list_path)
            scanner = self.get_scanner(mode, None)
            if file_exists is False:
                scanner.create_feature_filelist(data_folder_path)
            else:
                scanner.update_feature_filelist(new_files, data_folder_path)
