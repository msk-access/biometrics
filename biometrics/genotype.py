import os
from multiprocessing import Pool

import pandas as pd
import numpy as np
import plotly.graph_objects as go

from utils import exit_error

EPSILON = 1e-9


class Genotyper:

    def __init__(self, no_db_compare, discordance_threshold=0.05, threads=1, zmin=None, zmax=None):
        self.no_db_compare = no_db_compare
        self.discordance_threshold = discordance_threshold
        self.threads = threads
        self.zmax = zmax
        self.zmin = zmin
        self.sample_type_ratio = 1

    def are_samples_same_group(self, sample1, sample2):

        if sample1.group is None or sample2.group is None:
            return np.nan

        if sample1.group == sample2.group:
            return True
        else:
            return False

    def _plot_heatmap(self, data, outdir, name, size_ratio=None):

        width = None
        height = None

        if size_ratio is not None and size_ratio != 1:
            width = 2000
            height = (width * size_ratio)/2

        fig = go.Figure()
        fig.add_trace(
            go.Heatmap(
                x=data['ReferenceSample'],
                y=data['QuerySample'],
                z=data['DiscordanceRate'],
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
        fig.update_layout(
            yaxis_title="Query samples",
            xaxis_title="Reference samples",
            legend_title_text="Discordance",
            title_text="Discordance calculations between samples",
            width=width, height=height)
        fig.write_html(os.path.join(outdir, name))

        data = data[[
            'ReferenceSample', 'QuerySample', 'DatabaseComparison',
            'HomozygousInRef', 'TotalMatch', 'HomozygousMatch',
            'HeterozygousMatch', 'HomozygousMismatch', 'HeterozygousMismatch',
            'DiscordanceRate', 'Matched', 'ExpectedMatch', 'Status']]

    def plot(self, data, outdir):
        data_sub = data[~data['DatabaseComparison']].copy()

        if data_sub.shape[0] > 1:
            self._plot_heatmap(
                data_sub, outdir, 'genotype_comparison_input_only.html')

        data_sub = data[data['DatabaseComparison']].copy()

        if data_sub.shape[0] > 1:
            self._plot_heatmap(
                data_sub, outdir, 'genotype_comparison_database.html',
                self.sample_type_ratio)

    def _compute_discordance(self, samples):

        reference_sample, query_sample = samples

        row = {
            'ReferenceSample': reference_sample.name,
            'QuerySample': query_sample.name}

        row['HomozygousInRef'] = sum(reference_sample.pileup['genotype_class'] == 'Hom')
        row['TotalMatch'] = sum(reference_sample.pileup['genotype_class'] == query_sample.pileup['genotype_class'])
        row['HomozygousMatch'] = sum((reference_sample.pileup['genotype_class'] == query_sample.pileup['genotype_class']) & (reference_sample.pileup['genotype_class'] == 'Hom'))
        row['HeterozygousMatch'] = sum((reference_sample.pileup['genotype_class'] == query_sample.pileup['genotype_class']) & (reference_sample.pileup['genotype_class'] == 'Het'))
        row['HomozygousMismatch'] = sum((reference_sample.pileup['genotype'] != query_sample.pileup['genotype']) & ((reference_sample.pileup['genotype_class'] == 'Hom') & (query_sample.pileup['genotype_class'] == 'Hom')))
        row['HeterozygousMismatch'] = sum((reference_sample.pileup['genotype_class'] != query_sample.pileup['genotype_class']) & ((reference_sample.pileup['genotype_class'] == 'Het') | (query_sample.pileup['genotype_class'] == 'Het')))

        return row

    def genotype(self, samples):

        data = []
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

        thread_pool = Pool(self.threads)

        if sample_n_input > 1:
            # compare all the input samples to each other

            jobs = []

            for i, sample_name1 in enumerate(samples_input):
                for j, sample_name2 in enumerate(samples_input):
                    jobs.append([samples[sample_name1], samples[sample_name2]])

            results = thread_pool.map(self._compute_discordance, jobs)

            for i in range(len(results)):
                results[i]['DatabaseComparison'] = False
            data += results

        # for each input sample, compare with all the samples in the db

        if not self.no_db_compare and sample_n_db > 0:

            jobs = []

            for i, sample_name1 in enumerate(samples_input):
                for j, sample_name2 in enumerate(samples_db):
                    jobs.append([samples[sample_name1], samples[sample_name2]])

            results = thread_pool.map(self._compute_discordance, jobs)
            for i in range(len(results)):
                results[i]['DatabaseComparison'] = True
            data += results

        data = pd.DataFrame(data)

        # compute discordance rate

        data['DiscordanceRate'] = data['HomozygousMismatch'] / (data['HomozygousInRef'] + EPSILON)
        # data['DiscordanceRate'] = data['DiscordanceRate'].map(lambda x: round(x, 6))
        data.loc[data['HomozygousInRef'] < 10, 'DiscordanceRate'] = np.nan

        # for each comparison, indicate if the match/mismatch is expected
        # or not expected

        data['Matched'] = data['DiscordanceRate'] < self.discordance_threshold
        data['ExpectedMatch'] = data.apply(
            lambda x: self.are_samples_same_group(
                samples[x['ReferenceSample']],
                samples[x['QuerySample']]), axis=1)

        data['Status'] = ''
        data.loc[data.Matched & data.ExpectedMatch, 'Status'] = "Expected Match"
        data.loc[data.Matched & ~data.ExpectedMatch, 'Status'] = "Unexpected Match"
        data.loc[~data.Matched & data.ExpectedMatch, 'Status'] = "Unexpected Mismatch"
        data.loc[~data.Matched & ~data.ExpectedMatch, 'Status'] = "Expected Mismatch"

        data = data[[
            'ReferenceSample', 'QuerySample', 'DatabaseComparison',
            'HomozygousInRef', 'TotalMatch', 'HomozygousMatch',
            'HeterozygousMatch', 'HomozygousMismatch', 'HeterozygousMismatch',
            'DiscordanceRate', 'Matched', 'ExpectedMatch', 'Status']]

        return data
