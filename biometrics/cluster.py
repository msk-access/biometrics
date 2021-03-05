import pandas as pd
import networkx as nx

from biometrics.utils import get_logger

logger = get_logger()


class Cluster:

    def __init__(self, discordance_threshold=0.05):
        self.discordance_threshold = discordance_threshold

    def cluster(self, comparisons):

        assert comparisons is not None, "There is no fingerprint comparison data available."

        if len(comparisons) < 1:
            logger.warning('There are not enough comparisons to cluster.')
            return None

        sample2group = dict(zip(
            comparisons['ReferenceSample'], comparisons['ReferenceSampleGroup']))
        sample2group.update(dict(zip(
            comparisons['QuerySample'], comparisons['QuerySampleGroup'])))

        comparisons['is_same_group'] = comparisons['DiscordanceRate'].map(
            lambda x: 1 if x <= self.discordance_threshold else 0)

        graph = nx.from_pandas_edgelist(
            comparisons[comparisons['is_same_group']==1], 'ReferenceSample', 'QuerySample')

        clusters = []

        for cluster_idx, group in enumerate(nx.connected_components(graph)):
            samples_group = list(group)

            for i, sample in enumerate(samples_group):

                comparisons_sample = comparisons[
                    (comparisons['ReferenceSample']==sample) &
                    (comparisons['QuerySample']!=sample)]
                comparisons_cluster = comparisons_sample[
                    comparisons_sample['QuerySample'].isin(samples_group)]

                sample_status_counts = comparisons_sample['Status'].value_counts()
                cluster_status_counts = comparisons_cluster['Status'].value_counts()

                mean_discordance = 'NA'
                if len(comparisons_cluster) > 0:
                    mean_discordance = comparisons_cluster['DiscordanceRate'].mean()

                row = {
                    'sample_name': sample,
                    'expected_sample_group': sample2group[sample],
                    'cluster_index': cluster_idx,
                    'cluster_size': len(samples_group),
                    'avg_discordance': mean_discordance,
                    'count_expected_matches': cluster_status_counts.get('Expected Match', 0),
                    'count_unexpected_matches': cluster_status_counts.get('Unexpected Match', 0),
                    'count_expected_mismatches': sample_status_counts.get('Expected Mismatch', 0),
                    'count_unexpected_mismatches': sample_status_counts.get('Unexpected Mismatch', 0)
                }
                clusters.append(row)

        clusters = pd.DataFrame(clusters)

        logger.info(
            'Clustering finished. Grouped {} samples into {} clusters. Expected {} clusters.'.format(
            len(sample2group), cluster_idx + 1, len(set(clusters['expected_sample_group']))))

        return clusters
