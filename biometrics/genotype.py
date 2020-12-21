import os

import pandas as pd
import numpy as np
import plotly.graph_objects as go

from utils import exit_error

EPSILON = 1e-9


class Genotyper:

    def __init__(self, no_db_compare, discordance_threshold=0.05):
        self.no_db_compare = no_db_compare
        self.discordance_threshold = discordance_threshold

    def are_samples_same_group(self, sample1, sample2):

        if sample1.group is None or sample2.group is None:
            return np.nan

        if sample1.group == sample2.group:
            return True
        else:
            return False

    def _plot_heatmap(self, data, outdir, name):
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
                zmin=0,
                zmax=1,
                colorscale='Blues_r'
            ))
        fig.update_layout(
            yaxis_title="Query samples",
            xaxis_title="Reference samples",
            legend_title_text="Discordance",
            title_text="Discordance calculations between samples")
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
                data_sub, outdir, 'genotype_comparison_database.html')

    def compute_discordance(self, row, reference, query):

        row['HomozygousInRef'] = sum(reference.pileup['genotype_class'] == 'Hom')
        row['TotalMatch'] = sum(reference.pileup['genotype_class'] == query.pileup['genotype_class'])
        row['HomozygousMatch'] = sum((reference.pileup['genotype_class'] == query.pileup['genotype_class']) & (reference.pileup['genotype_class'] == 'Hom'))
        row['HeterozygousMatch'] = sum((reference.pileup['genotype_class'] == query.pileup['genotype_class']) & (reference.pileup['genotype_class'] == 'Het'))
        row['HomozygousMismatch'] = sum((reference.pileup['genotype_class'] != query.pileup['genotype_class']) & ((reference.pileup['genotype_class'] == 'Hom') | (query.pileup['genotype_class'] == 'Hom')))
        row['HeterozygousMismatch'] = sum((reference.pileup['genotype_class'] != query.pileup['genotype_class']) & ((reference.pileup['genotype_class'] == 'Het') | (query.pileup['genotype_class'] == 'Het')))

        return row

    def genotype(self, samples):

        data = []
        samples_db = dict(filter(lambda x: x[1].query_group, samples.items()))
        samples_input = dict(filter(
            lambda x: not x[1].query_group, samples.items()))

        if self.no_db_compare:
            if len(samples_input) <= 1:
                exit_error("You need to specify 2 or more samples in order to compare genotypes.")
        else:
            if len(samples_input) <= 1 and len(samples_db) < 1:
                exit_error("There are no samples in the database to compare with")

        if len(samples_input) > 1:
            # compare all the input samples to each other

            for i, sample_name1 in enumerate(samples_input):
                for j, sample_name2 in enumerate(samples_input):

                    row = {
                        'ReferenceSample': sample_name1,
                        'QuerySample': sample_name2,
                        'DatabaseComparison': False}
                    row = self.compute_discordance(
                        row, samples[sample_name1], samples[sample_name2])
                    data.append(row)

        # for each input sample, compare with all the samples in the db

        if not self.no_db_compare and len(samples_db) > 0:
            for i, sample_name1 in enumerate(samples_input):
                for j, sample_name2 in enumerate(samples_db):

                    row = {
                        'ReferenceSample': sample_name1,
                        'QuerySample': sample_name2,
                        'DatabaseComparison': True}
                    row = self.compute_discordance(
                        row, samples[sample_name1], samples[sample_name2])
                    data.append(row)

        data = pd.DataFrame(data)

        # compute discordance rate

        data['DiscordanceRate'] = data['HomozygousMismatch'] / (data['HomozygousInRef'] + EPSILON)
        data['DiscordanceRate'] = data['DiscordanceRate'].map(lambda x: round(x, 6))
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
