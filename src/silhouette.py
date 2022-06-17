from collections import Counter
from sklearn import mixture
from sklearn.metrics import silhouette_score


# From: https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html

class SilhouetteAnalysis:

    @staticmethod
    def perform_analysis(x, cluster_number_range=[3, 4, 5, 6, 7]):
        clustering_scores = []
        max_score = -1
        max_score_number = -1

        for clusters_number in cluster_number_range:

            clusterer = mixture.GaussianMixture(n_components=clusters_number, covariance_type='full',
                                                init_params='kmeans', random_state=42)
            clusterer.fit(x)
            cluster_labels = clusterer.predict(x)

            if (len(set(cluster_labels)) < 2):
                print("Single cluster found")
                continue

            # The silhouette_score gives the average value for all the samples.
            # This gives a perspective into the density and separation of the formed
            # clusters

            silhouette_average = silhouette_score(x, cluster_labels)
            print(''.join(
                ['Clusters(N=', repr(clusters_number), '): ', repr(Counter(cluster_labels)), ' | Silhouette score: ',
                 repr(silhouette_average)]))

            clustering_scores.append(silhouette_average)

            if (max_score < silhouette_average):
                max_score = silhouette_average
                max_score_number = clusters_number

        return max_score, max_score_number
