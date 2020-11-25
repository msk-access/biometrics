import pandas as pd
import numpy as np


class MinorContamination():
    """
    Minor contamination.
    """

    def __init__(self, args):
        pass

    def to_dataframe(self, samples):

        data = pd.DataFrame(
            columns=['sample', 'sample_group', 'sample_sex', 'sample_type',
                     'minor_contamination', 'total_homozygous_sites'])

        for sample_name, sample in samples.items():

            row = {
                'sample': sample_name,
                'sample_group': sample.group,
                'sample_sex': sample.sex,
                'sample_type': sample.sample_type,
                'minor_contamination': sample.metrics['minor_contamination'],
                'total_homozygous_sites': sample.metrics['total_homozygous_sites']
            }

            data = data.append(row, ignore_index=True)

        return data

    def estimate(self, samples):

        for sample_name, sample in samples.items():

            sites = samples[sample_name].pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            hom_sites = sites_notna[sites_notna['genotype_class'] == 'Hom']

            sample.metrics = {
                'total_homozygous_sites': hom_sites.shape[0]
            }

            if hom_sites.shape[0] == 0:
                sample.metrics['minor_contamination'] = np.nan
            else:
                sample.metrics['minor_contamination'] = \
                    round(hom_sites['minor_allele_freq'].mean(), 4)

        return samples
