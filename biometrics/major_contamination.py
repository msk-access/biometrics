import pandas as pd
import numpy as np


class MajorContamination():
    """
    Major contamination.
    """

    def __init__(self, args):
        pass

    def to_dataframe(self, samples):

        data = pd.DataFrame(
            columns=['sample', 'sample_group', 'sample_sex', 'sample_type',
                     'major_contamination', 'total_sites',
                     'total_heterozygous_sites'])

        for sample_name, sample in samples.items():

            row = {
                'sample': sample_name,
                'sample_group': sample.group,
                'sample_sex': sample.sex,
                'sample_type': sample.sample_type,
                'major_contamination': sample.metrics['major_contamination'],
                'total_sites': sample.metrics['total_sites'],
                'total_heterozygous_sites': sample.metrics['total_heterozygous_sites']
            }

            data = data.append(row, ignore_index=True)

        return data

    def estimate(self, samples):

        for sample_name, sample in samples.items():

            sites = sample.pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            het_sites = sites_notna[sites_notna['genotype_class'] == 'Het']

            sample.metrics = {
                'total_sites': sites_notna.shape[0],
                'total_heterozygous_sites': het_sites.shape[0]
            }

            if sites_notna.shape[0] == 0:
                sample.metrics['major_contamination'] = np.nan
            else:
                sample.metrics['major_contamination'] = \
                    round(het_sites.shape[0] / sites_notna.shape[0], 4)

        return samples
