from collections import defaultdict
import pandas as pd
import numpy as np
import os
import hdbscan
from src.enricher import Enricher
from umap import UMAP
from src.pdbhandler import PDBHandler
from tqdm import tqdm
from multiprocessing import Pool
import seaborn as sns
import ranky as rk
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys

plt.rcParams.update({'font.size': 8})


class CandidateSampler:

    def __init__(self):
        self.root_disk = ''
        self.debugging = False
        self.cores = 1
        self.reference_pdb_id = ''
        self.reference_chain_id = ''
        self.reference_gene_id = ''
        self.reference_segments = []
        # e.g. domain names:  self.referenceSegments = [['BetaCoV-S1-CTD', 'BetaCoV-S1-NTD']]
        # e.g. binding sites: self.referenceSegments = [[''.join(['_site', repr(x)]) for x in range(0, 5)]]
        self.no_enrichment = False
        self.excluded_organisms = []
        self.excluded_gene_names = []
        self.excluded_pdb_ids = []
        self._suffixed_filename = ''
        self.plots_output_path = 'plots'
        self.metrics_path = ''
        self.output_path = ''
        self.pdb_dataset_path = ''
        self._min_dataset_size = 300
        self.max_result_size = 250
        self._enricher = Enricher()
        self._column_headers = ['b-phipsi', 'w-rdist', 't-alpha']
        self._metric_order = [True, True, True]
        self._min_clustering_size = 100
        self._maximum_candidates_limit = 100
        self._max_clustering_attempts = 3
        self._extra_merging_columns = []
        self.keep_duplicates = False
        self.maximum_unsampled_dataset = 10000
        self.min_clustering_probability = 0.1
        self.override_pdb_id = ''

    # Filter an enriched entry by user specified criteria
    def filter_enriched_entry(self, entry):
        pdb_id, residues, resolution, organisms, genes = entry
        condition = 0
        if (organisms != '-'):
            organisms = [x.strip() for x in organisms.lower().split(',')]
            condition = min(len([x for x in self.excluded_organisms if x.lower() in organisms]), 1)
        if (condition == 0 and genes != '-'):
            genes = [x.strip() for x in genes.lower().split(',')]
            condition = min(len([x for x in self.excluded_gene_names if x.lower() in genes]), 1)
        return condition

    # Exclude all structures in the candidate set that do not match the user's preferences
    def filter_out(self, samples):
        print('Filtering structures based on user-specified organism, gene names and PDB IDs')
        samples['compositeId'] = [x.split('_') for x in samples['structureId'].to_list()]
        filtered = False
        # Exclude PDB IDs except for the reference PDB ID (if it was included in exclusions)
        if(len(self.excluded_pdb_ids) > 0):
            samples['excluded'] = [1 if (x[0] in self.excluded_pdb_ids and x[0] not in self.reference_pdb_id) else 0 for x in samples['compositeId'].to_list()]
            samples = samples[samples['excluded'] < 1].copy()
            del samples['excluded']
            filtered = True
        # Exclude organisms or genes
        if (len(self.excluded_organisms) > 0 or len(self.excluded_gene_names) > 0):
            results = []
            input_set = samples['compositeId'].to_list()
            if (self.debugging is False):
                with Pool(self.cores) as pool:
                    for result in tqdm(pool.imap(self._enricher.access_pdb_info, input_set), total=len(input_set), file=sys.stdout):
                        results.append(self.filter_enriched_entry(result))
                    pool.close()
                    pool.join()
            else:
                for compositeID in input_set:
                    results.append(self.filter_enriched_entry(self._enricher.access_pdb_info(compositeID)))
            samples['reject'] = results
            samples = samples[samples['reject'] < 1].copy()
            del samples['reject']
            filtered = True
        # Rank remaining samples
        if(filtered is True):
            samples = self.rank_samples(samples.copy(), ['structureId'])
        del samples['compositeId']
        samples.to_csv(''.join([self.output_path, os.path.sep, self._suffixed_filename, '-merged-filtered.csv']),
                       index=False, sep='\t')
        return samples

    # Determine the area of the reference protein where the features will be extracted (constrained mode)
    def set_segments(self, segment_type, target_segments=[]):
        # Retrieve the domains
        if (segment_type == 'domain'):
            domains = []
            pdbhandler = PDBHandler()
            pdbhandler.root_disk = self.root_disk
            pdbhandler.verbose = self.debugging
            pdbhandler.structure_id = self.reference_pdb_id
            pdb_path = ''.join([self.root_disk, os.path.sep, self.pdb_dataset_path, self.reference_pdb_id, '.pdb'])
            pdbhandler.get_uniprot_accession_number(self.reference_chain_id, self.reference_pdb_id, pdb_path)
            if (pdbhandler.uniprot_accession_number != ''):
                pdbhandler.get_domain_information()
                for domainInfo in pdbhandler.domains:
                    name, start, end = domainInfo
                    domain = pdbhandler.sanitize_domain_name(name)
                    domains.append(domain)
                self.reference_segments = domains
        # Determine the user-specified segments
        elif (segment_type == 'segment'):
            sites = []
            if(len(target_segments) > 0 and len(target_segments[0]) > 0):
                for site_index, segment in enumerate(target_segments):
                    sites.append(''.join(['_site', repr(site_index)]))
            self.reference_segments = sites
        if (len(self.reference_segments) == 0):
            raise Exception(''.join(['No segments were specified for ', self.reference_pdb_id]))

    def get_top_samples(self, data, percentage=0.01):
        filtered = None
        top_samples = int(data.index.shape[0] * percentage)
        # Collect top samples per metric and merge the collections, dropping the duplicates (intersection)
        for column_header in self._column_headers:
            if filtered is None:
                filtered = data.sort_values(by=column_header, ascending=True).head(top_samples).copy()
            else:
                filtered = filtered.append(
                    data.sort_values(by=column_header, ascending=column_header).head(top_samples).copy())
        filtered = self.rank_samples(filtered, ['structureId'])
        filtered.drop_duplicates(subset=['structureId'], keep='first', inplace=True)
        return filtered

    def rank_samples(self, data, merging_columns):
        data['rank'] = rk.borda(data[self._column_headers], reverse=True)
        # Sort the entries by their aggregated rank
        data = data.sort_values(by='rank', ascending=True)
        data.reset_index(drop=True, inplace=True)
        del data['rank']
        return data

    def select_candidates(self):
        # Initializations / path creations
        if (len(self.reference_segments) == 0):
            raise Exception('No segments are specified')
        if (os.path.exists(self.plots_output_path) is False):
            os.makedirs(self.plots_output_path)
        if (os.path.exists(self.output_path) is False):
            os.makedirs(self.output_path)

        self._enricher.set_root_disk(self.root_disk)
        self._enricher.pdb_dataset_path = self.pdb_dataset_path
        self._enricher.verbose = self.debugging
        self._enricher.debugging = self.debugging
        self._enricher.cores = self.cores

        # Select candidates for each reference segment
        for segment_index, reference_segment in enumerate(self.reference_segments):
            data = []
            suffix = ''
            reference_id = ''.join([self.reference_pdb_id, '_', self.reference_chain_id])
            data_headers = self._column_headers.copy()
            # filename suffixes for constrained mode (user-specified segments or domains)
            if reference_segment != '':
                if '_site' in reference_segment:
                    suffix = ''.join([reference_segment, '-metrics'])
                    data_headers.extend(['alignedRange', 'referenceSegmentIndex', 'referenceRange', 'alignedSequences'])
                else:
                    suffix = ''.join(['_', reference_segment])
                    data_headers.extend(['domain'])
                    self._extra_merging_columns = ['domain']
            self._suffixed_filename = ''.join([self.reference_pdb_id, '_', self.reference_chain_id, suffix])
            final_output_path = ''.join([self.output_path, os.path.sep, self._suffixed_filename, '-merged-enriched.csv'])
            notenriched_output_path = ''.join([self.output_path, os.path.sep, self._suffixed_filename, '-merged-notenriched.csv'])
            # if (os.path.exists(final_output_path)):
            #     os.unlink(final_output_path)
            # if (os.path.exists(notenriched_output_path)):
            #     os.unlink(notenriched_output_path)
            if (os.path.exists(final_output_path) or os.path.exists(notenriched_output_path)):
                print(''.join(['Candidates are already selected. If you wish to repeat the process, delete the following folder and execute the method again:\n',
                                self.output_path]))
            else:
                if '_site' not in reference_segment:
                    # Aggregate metrics to a single dataframe
                    for metric in self._column_headers:
                        data_filename = ''.join(
                            [self.metrics_path, os.path.sep, self._suffixed_filename, '_', metric, '.csv'])
                        data.append(pd.read_csv(data_filename, sep=','))
                    merging_columns = ['structureId']
                    merging_columns.extend(self._extra_merging_columns)
                    merged_data = data[0]
                    for metric_index in range(1, len(self._column_headers)):
                        slice_columns = ['structureId', self._column_headers[metric_index]]
                        slice_columns.extend((self._extra_merging_columns))
                        merged_data = pd.merge(merged_data, data[metric_index][slice_columns],
                                               how='inner', on=merging_columns)
                else:
                    data_filename = ''.join([self.metrics_path, os.path.sep, self._suffixed_filename, '.csv'])
                    merged_data = pd.read_csv(data_filename, sep=',')
                print(''.join([repr(merged_data.index.shape[0]), ' proteins for sampling']))

                # Remove entries with missing values
                merged_data = merged_data.replace([np.inf, -np.inf], np.nan).dropna()

                # Remove other domains from reference during domain scanning
                if ('domain' in data_headers):
                    indices = merged_data[
                        (merged_data['structureId'] == reference_id) & (merged_data['domain'] != reference_segment)].index
                    merged_data.drop(indices, inplace=True)

                # Sort entries by rank
                merged_data = self.rank_samples(merged_data, ['structureId'])
                subset = ['structureId']
                subset.extend(self._extra_merging_columns)

                # Remove duplicates
                merged_data.drop_duplicates(subset=subset, keep='first', inplace=True)

                # Exclude reference from candidate set
                if (reference_segment != '' and '_site' not in reference_segment):
                    reference_entry = merged_data[
                        (merged_data['structureId'] == reference_id) & (merged_data['domain'] == reference_segment)].copy()
                    merged_data = merged_data[
                        (merged_data['structureId'] != reference_id) & (merged_data['domain'] != reference_segment)]
                else:
                    reference_entry = merged_data[merged_data['structureId'] == reference_id].copy()
                    merged_data = merged_data[merged_data['structureId'] != reference_id]

                # Filter entries
                if (len(self.excluded_organisms) > 0 or len(self.excluded_gene_names) > 0 or len(
                        self.excluded_pdb_ids) > 0):
                    merged_data = self.filter_out(merged_data)

                # Get top samples per metric
                if (merged_data.index.shape[0] > self.maximum_unsampled_dataset):
                    samples = self.get_top_samples(merged_data)
                else:
                    samples = merged_data

                # Cluster entries by considering the resulting cluster's size
                # Repeat the process with a different size threshold if the
                # cluster is smaller than desired.
                available_attempts = 3
                min_cluster_size = self._min_clustering_size
                metric_values = None
                embedded = None
                total_samples = 0
                if (samples.index.shape[0] > self._min_dataset_size):
                    while (available_attempts > 0):
                        samples.reset_index(drop=True, inplace=True)
                        # Prune clustering space
                        samples = self.rank_samples(samples, ['structureId'])
                        # Visualize space of candidates
                        if (available_attempts == self._max_clustering_attempts):
                            metric_values = samples[self._column_headers]
                            embedded = self.visualize_candidates(metric_values)
                            total_samples = samples.index.shape[0]
                        else:
                            print('Data were not clustered. Retrying...')
                        clustering_info = self.cluster_data(metric_values, min_cluster_size)
                        # If there is at least one cluster
                        # ('0' for first cluster and not '-1' only labels that refer to noise)
                        if (0 in clustering_info['labels']):
                            # Store clustering result
                            previous_samples = samples.copy()
                            samples['cluster'] = clustering_info['labels'].tolist()
                            samples['originalLabels'] = clustering_info['originalLabels'].tolist()
                            samples['outlierScores'] = clustering_info['outlierScores'].tolist()
                            samples['probabilities'] = clustering_info['probabilities'].tolist()
                            samples.to_csv(notenriched_output_path.replace('notenriched', 'allclusters'), index=False, sep='\t')
                            # Visualize clusters
                            self.plot_clusters(samples, embedded, clustering_info, '-initial')

                            # Treat entries clustered with low probability as noise
                            labels = clustering_info['labels'].tolist()
                            clustering_info['labels'] = np.array(
                                [-1 if clustering_info['probabilities'][
                                           index] < self.min_clustering_probability and x > 0 else x for index, x in
                                 enumerate(labels)])

                            # Traverse the candidates starting with the one of the highest order (sorted list)
                            # and check their cluster labels. Pick the first positive cluster label (not noise)
                            # and keep also the previous entries.
                            selected_cluster = -1
                            traversed = 0
                            for i in range(clustering_info['minClusterSize']):
                                selected_cluster = samples['cluster'].iloc[i]
                                traversed = i
                                if (selected_cluster != -1):
                                    break
                            print('Samples top 15% count: ',repr(int(len(samples.index) * 0.15)), ' / traversed looking for cluster: ', repr(traversed))
                            top_rows_count = max(int(len(samples.index) * 0.15), traversed)
                            for i in range(top_rows_count):
                                if(samples.loc[i, 'cluster'] != selected_cluster):
                                    samples.loc[i, 'cluster'] = selected_cluster
                                    clustering_info['labels'][i] = selected_cluster
                            samples = samples[samples['cluster'] == selected_cluster].copy()
                            clustering_info['selectedCluster'] = selected_cluster
                            if (samples.index.shape[0] < self._min_dataset_size):
                                samples = previous_samples
                            else:
                                # Visualize final result
                                self.plot_clusters(samples, embedded, clustering_info)
                            break
                        if (samples.index.shape[0] == total_samples):
                            available_attempts -= 1
                            if (available_attempts <= 0):
                                print('Clustering attempts are ended.')
                                break
                            elif (available_attempts > 0):
                                min_cluster_size -= 20
                        else:
                            break
                    # Store the remaining entries
                    samples.to_csv(notenriched_output_path.replace('notenriched', 'cluster'), index=False, sep='\t')
                    samples = self.rank_samples(samples, ['structureId'])
                samples = samples.iloc[:self.max_result_size].copy()
                checkpoint = pd.concat([reference_entry, samples]).copy()
                checkpoint.to_csv(notenriched_output_path, index=False, sep='\t')
            # Enrich the current set of candidate entries
            if (self.no_enrichment is False and os.path.exists(final_output_path) is False):
                print('Enriching candidate selections...')
                if (len(data) == 0):
                    if(os.path.exists(notenriched_output_path)):
                        samples = pd.read_csv(notenriched_output_path, sep='\t')
                        reference_entry = samples[samples.index==0]
                        samples = samples[samples.index>0]
                        samples.reset_index(drop=True, inplace=True)
                    else:
                        samples = None
                if(samples is not None):
                    self.enrich_samples(samples, reference_entry, data_headers)

    # Visualize candidate space using UMAP
    def visualize_candidates(self, metric_values):
        result_filename = ''.join([self.plots_output_path, os.path.sep, self._suffixed_filename, '-UMAP.npy'])
        umap_2d = UMAP(n_components=2, init='random', random_state=42, min_dist=0.25)
        embedding = umap_2d.fit_transform(metric_values)
        with open(result_filename, 'wb') as resultFile:
            np.save(resultFile, embedding)

        # Save plot in png and eps formats
        plt.clf()
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))
        plt.scatter(embedding[:, 0], embedding[:, 1])
        plt.savefig(''.join([self.plots_output_path, os.path.sep, self._suffixed_filename, '-UMAP.png']),
                    format='png',
                    dpi=300, bbox_inches='tight')
        plt.savefig(''.join([self.plots_output_path, os.path.sep, self._suffixed_filename, '-UMAP.eps']),
                    format='eps',
                    dpi=300, bbox_inches='tight')
        plt.clf()
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))
        return embedding

    def cluster_data(self, metric_values, min_cluster_size=100):
        total = metric_values.index.shape[0]
        print(''.join(['Attempting to cluster data: total samples: ', repr(total), ' / min cluster size: ',
                       repr(min_cluster_size)]))
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=5)
        clusterer.fit(metric_values)
        clustering_info = defaultdict()
        clustering_info['totalClusters'] = clusterer.labels_.max() + 1
        clustering_info['probabilities'] = clusterer.probabilities_
        clustering_info['outlierScores'] = clusterer.outlier_scores_
        clustering_info['minClusterSize'] = min_cluster_size
        clustering_info['originalLabels'] = clusterer.labels_
        clustering_info['labels'] = clusterer.labels_

        return clustering_info

    def plot_clusters(self, samples, embedding, clustering_info, file_suffix=''):
        print(''.join(['Cluster size: ', repr(len(samples.index)), ' total clusters: ',
                       repr(clustering_info['totalClusters'])]))

        # Assign a color to each cluster
        color_palette = sns.color_palette('colorblind', clustering_info['totalClusters'])
        cluster_colors = []
        coloring = defaultdict()
        for label in clustering_info['labels']:
            sample_color = color_palette[label] if label >= 0 else (0.3, 0.3, 0.3)
            cluster_colors.append(sample_color)
            if label not in coloring:
                coloring[label] = sample_color
        cluster_colors = np.array(cluster_colors)

        # Change opacity of the plotted data points according to their clustering probability
        bins = [x for x in np.arange(0.0, 1.1, 0.1)]
        binned_indices = np.digitize(clustering_info['probabilities'], bins)
        plt.clf()
        for bin_index in range(len(bins)):
            bin_indices = np.argwhere(binned_indices == bin_index + 1)
            bin_indices = bin_indices.reshape(len(bin_indices))
            plt.scatter(*embedding[bin_indices].T, s=50, linewidth=0, c=cluster_colors[bin_indices],
                        alpha=(max(bins[bin_index], 0.01) if file_suffix=='' else 1.0))

        # Create legends for the plot
        legend_elements = []
        sorted_labels = sorted(list(coloring.keys()))
        for label in sorted_labels:
            legend_label = 'Noise'
            if (int(label) >= 0):
                if('selectedCluster' in clustering_info):
                    legend_label = ''.join(
                        ['Cluster ', repr(label), '' if label != clustering_info['selectedCluster'] else ' (Top)'])
                else:
                    legend_label = ''.join(['Cluster ', repr(label)])
            legend_elements.append(
                Line2D([0], [0], marker='o', color='w', label=legend_label, markerfacecolor=coloring[label],
                       markersize=10))
        plt.legend(handles=legend_elements, facecolor='white', framealpha=1)

        # Save plot in png and eps formats
        plt.savefig(''.join([self.plots_output_path, os.path.sep, self._suffixed_filename, '-clusters', file_suffix, '.png']),
                    format='png', dpi=300)
        ax = plt.gca()
        ax.set_rasterized(True)
        plt.savefig(''.join([self.plots_output_path, os.path.sep, self._suffixed_filename, '-clusters', file_suffix, '.eps']),
                    format='eps', dpi=300)
        plt.clf()
        plt.figure(figsize=plt.rcParams.get('figure.figsize'))

    def enrich_samples(self, samples, reference_entry, column_headers):
        # Enrich candidate entries
        self._enricher.load_go_cache()
        self._enricher.file_name = ''.join([self.output_path, os.path.sep, self._suffixed_filename, '-merged.csv'])
        self._enricher.cores = self.cores
        if(len(self.override_pdb_id) > 0):
            # Use a different PDB ID for enrichment (custom PDBs)
            overriden_reference = reference_entry.copy()
            overriden_reference['structureId'].iloc[0] = '_'.join([self.override_pdb_id, self.reference_chain_id])
            samples = pd.concat([overriden_reference, samples])
        else:
            samples = pd.concat([reference_entry, samples])
        self._enricher.fetch_enrichment(samples, column_headers)

        enriched = pd.read_csv(self._enricher.file_name.replace('.csv', '-enriched.csv'), sep='\t')

        # Get ennriched reference entry (restore reference PDB ID if an alternative was used for enrichment)
        reference_entry = enriched[(enriched['pdbId'] == (self.override_pdb_id if(len(self.override_pdb_id) > 0) else self.reference_pdb_id)) &
                                   (enriched['chainId'] == self.reference_chain_id)].copy()
        if (len(self.override_pdb_id) > 0):
            pd.options.mode.chained_assignment = None
            reference_entry['pdbId'].iloc[0] = self.reference_pdb_id
            pd.options.mode.chained_assignment = 'warn'

        if(self.keep_duplicates is False):
            # Exclude reference gene and remove multiple rows for each gene (keep the best one)
            # Filtering ignores the entries without gene ids
            if(self.reference_gene_id != ''):
                enriched = enriched[(enriched['geneId'] != self.reference_gene_id)]
                enriched = enriched[(enriched['geneId'] != int(self.reference_gene_id))]
            enriched = pd.concat([enriched[enriched['geneId'] == '#'], enriched[enriched['geneId'] != '#'].drop_duplicates(subset=['geneId'], keep='first')])

        # Remove duplicate reference entry if an alternative PDB ID was used during enrichment
        if (len(self.override_pdb_id) > 0):
            enriched.drop(enriched[enriched['pdbId'] == self.override_pdb_id].index, inplace=True)

        # Rank and prune entries. Include reference entry in the results as the first entry.
        enriched = self.rank_samples(enriched, ['pdbId', 'chainId'])
        enriched = enriched.head(self._maximum_candidates_limit)
        if(len(enriched[(enriched['pdbId'] == self.reference_pdb_id) & (enriched['chainId'] == self.reference_chain_id)].index) == 0):
            enriched = pd.concat([reference_entry, enriched])

        enriched.to_csv(self._enricher.file_name.replace('.csv', '-enriched.csv'), index=False, sep='\t')
