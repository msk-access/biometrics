import pandas as pd
import numpy as np


class MajorContamination:
    """
    Major contamination.
    """

    def __init__(self, args):
        pass

    def estimate(self, samples):

        for i, sample in enumerate(samples):

            sites = sample.pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            het_sites = sites_notna[sites_notna['genotype_class']=='Het']

            if sites_notna.shape[0] == 0:
                samples[i].metrics['major_contamination'] = np.nan
            else:
                samples[i].metrics['major_contamination'] = \
                    het_sites.shape[0] / sites_notna.shape[0]

        return samples
