import pandas as pd
import numpy as np

EPSILON = 1e-9


class Genotyper:

    def __init__(self, args):
        self.no_db_comparison = args.no_db_comparison

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
        samples_db = list(filter(lambda x: x.is_in_db, samples))
        samples_input = list(filter(lambda x: not x.is_in_db, samples))

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

        data['DiscordanceRate'] = data['HomozygousMismatch'] / (data['HomozygousInRef'] + EPSILON)
        data.loc[data['HomozygousInRef'] < 10, 'DiscordanceRate'] = np.nan

        return data
