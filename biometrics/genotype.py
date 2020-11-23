import pandas as pd
import numpy as np

from utils import exit_error

EPSILON = 1e-9
DISCORDANCE_THRES = 0.05


class Genotyper:

    def __init__(self, args):
        self.no_db_comparison = args.no_db_comparison

    def are_samples_same_group(self, sample1, sample2):

        if sample1.group is None or sample2.group is None:
            return np.nan

        if sample1.group == sample2.group:
            return True
        else:
            return False

    def get_sample_groups(self, samples):

        groups = {}

        for sample_name, sample in samples.items():
            if sample.group is None:
                continue

            if sample.group not in groups:
                groups[sample.group] = set(sample_name)
            else:
                groups[sample.group].add(sample_name)

        return groups

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
        samples_db = dict(filter(lambda x: x[1].is_in_db, samples.items()))
        samples_input = dict(filter(
            lambda x: not x[1].is_in_db, samples.items()))

        if not self.no_db_comparison:
            if len(samples_input) <= 1:
                exit_error("You need to specify 2 or more samples in order to compare genotypes.")
        else:
            if len(samples_db) <= 1:
                exit_error("There are no samples in the database to compare with")

        # compare all the input samples to each other

        for i, sample_name1 in enumerate(samples_input):
            for j, sample_name2 in enumerate(samples_input):

                if i == j:
                    continue

                row = {
                    'ReferenceSample': sample_name1,
                    'QuerySample': sample_name2}
                row = self.compute_discordance(
                    row, samples[sample_name1], samples[sample_name2])
                data.append(row)

        # for each input sample, compare with all the samples in the db

        if not self.no_db_comparison and len(samples_db) > 0:
            for i, sample_name1 in enumerate(samples_input):
                for j, sample_name2 in enumerate(samples_db):

                    if i == j:
                        continue

                    row = {
                        'ReferenceSample': sample_name1,
                        'QuerySample': sample_name2}
                    row = self.compute_discordance(
                        row, samples[sample_name1], samples[sample_name2])
                    data.append(row)

        data = pd.DataFrame(data)

        # compute discordance rate

        data['DiscordanceRate'] = data['HomozygousMismatch'] / (data['HomozygousInRef'] + EPSILON)
        data.loc[data['HomozygousInRef'] < 10, 'DiscordanceRate'] = np.nan

        # for each comparison, indicate if the match/mismatch is expected
        # or not expected

        data['Matched'] = data['DiscordanceRate'] < DISCORDANCE_THRES
        data['ExpectedMatch'] = data.apply(
            lambda x: self.are_samples_same_group(
                samples[x['ReferenceSample']],
                samples[x['QuerySample']]), axis=1)

        data['Status'] = ''
        data.loc[data.Matched & data.ExpectedMatch, 'Status'] = "Expected Match"
        data.loc[data.Matched & ~data.ExpectedMatch, 'Status'] = "Unexpected Match"
        data.loc[~data.Matched & data.ExpectedMatch, 'Status'] = "Unexpected Mismatch"
        data.loc[~data.Matched & ~data.ExpectedMatch, 'Status'] = "Expected Mismatch"

        return data
