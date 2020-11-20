import pandas as pd
import numpy as np

EPSILON = 1e-9


class Genotyper:

    def __init__(self, args):
        pass

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

        for i, sample1 in enumerate(samples):
            for j, sample2 in enumerate(samples):

                if i == j:
                    continue

                row = {'ReferenceSample': sample1.name, 'QuerySample': sample2.name}
                row = self.compute_discordance(row, sample1, sample2)
                data.append(row)

        data = pd.DataFrame(data)
        data['DiscordanceRate'] = data['HomozygousMismatch'] / (data['HomozygousInRef'] + EPSILON)
        data.loc[data['HomozygousInRef'] < 10, 'DiscordanceRate'] = np.nan

        import pdb; pdb.set_trace()

        return data
