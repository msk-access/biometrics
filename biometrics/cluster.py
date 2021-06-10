from collections import Counter

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
            lambda x: 1 if ((x <= self.discordance_threshold) & (~pd.isna(x))) else 0)

        graph = nx.from_pandas_edgelist(
            comparisons[comparisons['is_same_group']==1], 'ReferenceSample', 'QuerySample')

        clusters = []
        cluster_idx = 0

        for cluster_idx, group in enumerate(nx.connected_components(graph)):
            samples_in_group = list(group)
            sample_groups = list(set(
                comparisons[comparisons['ReferenceSample'].isin(samples_in_group)]['ReferenceSampleGroup']))

            occurences = Counter(sample_groups)
            max_occurence = occurences.most_common()
            max_occurence_val = max_occurence[0][1]
            most_common_group = list(filter(lambda x: x[1] == max_occurence_val, max_occurence))
            most_common_group = ','.join([i[0] for i in most_common_group])

            for i, sample in enumerate(samples_in_group):

                comparisons_sample = comparisons[
                    (comparisons['ReferenceSample']==sample) &
                    (comparisons['QuerySample']!=sample)]
                comparisons_cluster = comparisons_sample[
                    comparisons_sample['QuerySample'].isin(samples_in_group)]

                sample_status_counts = comparisons_sample['Status'].value_counts()
                cluster_status_counts = comparisons_cluster['Status'].value_counts()

                mean_discordance = 'NA'
                if len(comparisons_cluster) > 0:
                    mean_discordance = comparisons_cluster['DiscordanceRate'].mean()

                row = {
                    'sample_name': sample,
                    'expected_sample_group': sample2group[sample],
                    'predicted_sample_group': most_common_group,
                    'cluster_index': cluster_idx,
                    'cluster_size': len(samples_in_group),
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
