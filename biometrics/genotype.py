import os
from multiprocessing import Pool

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import networkx as nx

from biometrics.utils import exit_error

EPSILON = 1e-9


class Genotyper:

    def __init__(self, no_db_compare, discordance_threshold=0.05, threads=1, zmin=None, zmax=None):
        self.no_db_compare = no_db_compare
        self.discordance_threshold = discordance_threshold
        self.threads = threads
        self.zmax = zmax
        self.zmin = zmin
        self.sample_type_ratio = 1
        self.comparisons = None
        self.clusters = None

    def are_samples_same_group(self, sample1, sample2):

        if sample1.sample_group is None or sample2.sample_group is None:
            return np.nan

        if sample1.sample_group == sample2.sample_group:
            return True
        else:
            return False

    def _plot_heatmap(self, data, outdir, name, title="Discordance calculations between samples", size_ratio=None):

        width = None
        height = None

        if size_ratio is not None and size_ratio != 1:
            width = 1400
            height = (width * size_ratio)/4

        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=data['ReferenceSample'],
                y=data['QuerySample'],
                z=data['DiscordanceRate'],
                legendgroup="Discordance",
                name='Discordance',
                customdata=data.to_numpy(),
                hovertemplate='<b>Reference sample:</b> %{customdata[0]}' +
                              '<br><b>Query sample:</b> %{customdata[1]}' +
                              '<br><b>Homozygous count in reference:</b> %{customdata[3]}' +
                              '<br><b>Total match count:</b> %{customdata[4]}' +
                              '<br><b>Homozygous match count:</b> %{customdata[5]}' +
                              '<br><b>Heterozygous match count:</b> %{customdata[6]}' +
                              '<br><b>Homozygous mismatch count:</b> %{customdata[7]}' +
                              '<br><b>Heterozygous mismatch count:</b> %{customdata[8]}' +
                              '<br><b>Discordance rate:</b> %{customdata[9]}' +
                              '<br><b>Status:</b> %{customdata[12]}' +
                              '<extra></extra>',
                zmin=self.zmin,
                zmax=self.zmax,
                colorscale='Blues_r'
            ))

        # add red dots to sample pairs that are unexpected match/mismatch

        data_sub = data[(data['Status']=='Unexpected Match') | (data['Status']=='Unexpected Mismatch')].copy()

        if len(data_sub) > 0:
            fig.add_trace(
                go.Scatter(
                    mode="markers",
                    x=data_sub['ReferenceSample'],
                    y=data_sub['QuerySample'],
                    marker_symbol=[17],
                    marker_color="red",
                    marker_line_width=0,
                    marker_size=10,
                    customdata=data_sub.to_numpy(),
                    hovertemplate='%{customdata[12]}<extra></extra>'))

        fig.update_layout(
            yaxis_title="Query samples",
            xaxis_title="Reference samples",
            legend_title_text="Discordance",
            title_text=title,
            width=width, height=height)
        fig.write_html(os.path.join(outdir, name))

        data = data[[
            'ReferenceSample', 'QuerySample', 'DatabaseComparison',
            'HomozygousInRef', 'TotalMatch', 'HomozygousMatch',
            'HeterozygousMatch', 'HomozygousMismatch', 'HeterozygousMismatch',
            'DiscordanceRate', 'Matched', 'ExpectedMatch', 'Status']]

    def plot(self, data, outdir):

        # make plot for comparing input samples with each other

        data_sub = data[~data['DatabaseComparison']].copy()
        data_sub['DiscordanceRate'] = data_sub['DiscordanceRate'].map(
            lambda x: round(x, 4))

        if data_sub.shape[0] > 1:
            self._plot_heatmap(
                data_sub, outdir, name='genotype_comparison_input_only.html',
                title="Discordance calculations between input samples")

        # make plot for comparing input samples with database samples

        data_sub = data[data['DatabaseComparison']].copy()
        data_sub['DiscordanceRate'] = data_sub['DiscordanceRate'].map(
            lambda x: round(x, 4))

        if data_sub.shape[0] > 1:
            self._plot_heatmap(
                data_sub, outdir, name='genotype_comparison_database.html',
                title="Discordance calculations between input samples and database samples",
                size_ratio=self.sample_type_ratio)

    def _compute_discordance(self, reference_sample, query_sample):
        """
        Compute discordance between two samples
        """

        pileup_ref = reference_sample.pileup
        pileup_query = query_sample.pileup

        common_covered = ~pd.isna(pileup_ref['genotype_class']) & ~pd.isna(pileup_query['genotype_class'])
        pileup_ref = pileup_ref[common_covered]
        pileup_query = pileup_query[common_covered]

        row = {
            'ReferenceSample': reference_sample.sample_name,
            'ReferenceSampleGroup': reference_sample.sample_group,
            'QuerySample': query_sample.sample_name,
            'QuerySampleGroup': query_sample.sample_group}

        if len(pileup_query) > 0:
            row['HomozygousInRef'] = sum(pileup_ref['genotype_class'] == 'Hom')
            row['TotalMatch'] = sum(pileup_ref['genotype_class'] == pileup_query['genotype_class'])
            row['HomozygousMatch'] = sum((pileup_ref['genotype_class'] == pileup_query['genotype_class']) & (pileup_ref['genotype_class'] == 'Hom'))
            row['HeterozygousMatch'] = sum((pileup_ref['genotype_class'] == pileup_query['genotype_class']) & (pileup_ref['genotype_class'] == 'Het'))
            row['HomozygousMismatch'] = sum((pileup_ref['genotype'] != pileup_query['genotype']) & ((pileup_ref['genotype_class'] == 'Hom') & (pileup_query['genotype_class'] == 'Hom')))
            row['HeterozygousMismatch'] = sum((pileup_ref['genotype_class'] != pileup_query['genotype_class']) & ((pileup_ref['genotype_class'] == 'Het') | (pileup_query['genotype_class'] == 'Het')))
            row['CountOfCommonSites'] = len(pileup_query)
        else:
            # if there are no regions with enough coverage

            row['HomozygousInRef'] = np.nan
            row['TotalMatch'] = np.nan
            row['HomozygousMatch'] = np.nan
            row['HeterozygousMatch'] = np.nan
            row['HomozygousMismatch'] = np.nan
            row['HeterozygousMismatch'] = np.nan
            row['CountOfCommonSites'] = 0

        return row

    def _compute_discordance_batch_job(self, batch_sample_pairs):
        """
        batch job to take a batch of pairs of samples to compute discordance
        """
        rows = []

        for reference_sample, query_sample in batch_sample_pairs:
            rows.append(self._compute_discordance(
                reference_sample, query_sample))

        return rows

    def _compare_sample_lists(self, sample_set1, sample_set2, samples):
        """
        Compare two lists of samples. Does so by first grouping the sample
        comparisons into batches for parallel processing
        """

        jobs = []
        thread_pool = Pool(self.threads)
        total_comparisons = len(sample_set1) * len(sample_set2)
        parallel_batch_size = max(int(total_comparisons / self.threads), 1)

        current_batch = []

        for i, sample_name1 in enumerate(sample_set1):
            for j, sample_name2 in enumerate(sample_set2):
                current_batch.append([samples[sample_name1], samples[sample_name2]])

                if len(current_batch) >= parallel_batch_size:
                    jobs.append(current_batch)
                    current_batch = []

        # add any remaining jobs

        if len(current_batch) > 0:
            jobs.append(current_batch)

        # analyze, collect results, and flatten the lists

        results = thread_pool.map(self._compute_discordance_batch_job, jobs)
        results = [item for sublist in results for item in sublist]

        return results


    def cluster_samples(self, samples):

        samples_input = dict(filter(
            lambda x: not x[1].query_group, samples.items()))

        if len(samples_input) == 0 or self.comparisons is None:
            return

        samples_input_names = [i.sample_name for i in samples_input.values()]

        comparisons = self.comparisons[self.comparisons['ReferenceSample'].isin(samples_input_names)]
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
                    'expected_sample_group': samples_input[sample].sample_group,
                    'cluster_index': cluster_idx,
                    'cluster_size': len(samples_group),
                    'avg_discordance': mean_discordance,
                    'count_expected_matches': cluster_status_counts.get('Expected Match', 0),
                    'count_unexpected_matches': cluster_status_counts.get('Unexpected Match', 0),
                    'count_expected_mismatches': sample_status_counts.get('Expected Mismatch', 0),
                    'count_unexpected_mismatches': sample_status_counts.get('Unexpected Mismatch', 0)
                }
                clusters.append(row)

        self.clusters = pd.DataFrame(clusters)
        return self.clusters

    def compare_samples(self, samples):

        comparisons = []
        samples_db = dict(filter(lambda x: x[1].query_group, samples.items()))
        samples_input = dict(filter(
            lambda x: not x[1].query_group, samples.items()))

        # get the number of each type of sample and compute a ratio
        # this is used to plot the heatmap when comparing with database
        # samples

        sample_n_db = len(samples_db)
        sample_n_input = len(samples_input)

        if sample_n_db > 0 and sample_n_input > 0 and sample_n_db > sample_n_input:
            self.sample_type_ratio = sample_n_db / sample_n_input

        # check to see if there are appropriate number of samples to
        # do the analysis

        if self.no_db_compare:
            if len(samples_input) <= 1:
                exit_error("You need to specify 2 or more samples in order to compare genotypes.")
        else:
            if len(samples_input) <= 1 and len(samples_db) < 1:
                exit_error("There are no samples in the database to compare with")

        # compare all the input samples to each other

        if sample_n_input > 1:
            results = self._compare_sample_lists(
                samples_input, samples_input, samples)

            for i in range(len(results)):
                results[i]['DatabaseComparison'] = True
            comparisons += results

        # for each input sample, compare with all the samples in the db

        if not self.no_db_compare and sample_n_db > 0:
            results = self._compare_sample_lists(
                samples_input, samples_db, samples)

            for i in range(len(results)):
                results[i]['DatabaseComparison'] = True
            comparisons += results

        comparisons = pd.DataFrame(comparisons)

        # compute discordance rate

        comparisons['DiscordanceRate'] = comparisons['HomozygousMismatch'] / (comparisons['HomozygousInRef'] + EPSILON)
        # data['DiscordanceRate'] = data['DiscordanceRate'].map(lambda x: round(x, 6))
        comparisons.loc[comparisons['HomozygousInRef'] < 10, 'DiscordanceRate'] = np.nan

        # for each comparison, indicate if the match/mismatch is expected
        # or not expected

        comparisons['Matched'] = comparisons['DiscordanceRate'] < self.discordance_threshold
        comparisons['ExpectedMatch'] = comparisons.apply(
            lambda x: self.are_samples_same_group(
                samples[x['ReferenceSample']],
                samples[x['QuerySample']]), axis=1)

        comparisons['Status'] = ''
        comparisons.loc[comparisons.Matched & comparisons.ExpectedMatch, 'Status'] = "Expected Match"
        comparisons.loc[comparisons.Matched & ~comparisons.ExpectedMatch, 'Status'] = "Unexpected Match"
        comparisons.loc[~comparisons.Matched & comparisons.ExpectedMatch, 'Status'] = "Unexpected Mismatch"
        comparisons.loc[~comparisons.Matched & ~comparisons.ExpectedMatch, 'Status'] = "Expected Mismatch"

        self.comparisons = comparisons[[
            'ReferenceSample', 'ReferenceSampleGroup', 'QuerySample', 'QuerySampleGroup',
            'CountOfCommonSites', 'HomozygousInRef', 'TotalMatch', 'HomozygousMatch',
            'HeterozygousMatch', 'HomozygousMismatch',
            'HeterozygousMismatch', 'DiscordanceRate', 'Matched',
            'ExpectedMatch', 'Status']]

        return self.comparisons
