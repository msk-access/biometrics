import pandas as pd
import numpy as np


class MinorContamination:
    """
    Minor contamination.
    """

    def __init__(self, args):
        pass

    def estimate(self, samples):

        for i, sample_name in enumerate(samples):

            sites = samples[sample_name].pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            hom_sites = sites_notna[sites_notna['genotype_class'] == 'Hom']

            if hom_sites.shape[0] == 0:
                samples[sample_name].metrics['minor_contamination'] = np.nan
            else:
                samples[sample_name].metrics['minor_contamination'] = \
                    hom_sites['minor_allele_freq'].mean()

        return samples
