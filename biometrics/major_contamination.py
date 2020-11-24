import pandas as pd
import numpy as np

from contamination import BaseContamination


class MajorContamination(BaseContamination):
    """
    Major contamination.
    """

    def __init__(self, args):
        pass

    def estimate(self, samples):

        for i, sample_name in enumerate(samples):

            sites = samples[sample_name].pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            het_sites = sites_notna[sites_notna['genotype_class'] == 'Het']

            if sites_notna.shape[0] == 0:
                samples[sample_name].metrics['major_contamination'] = np.nan
            else:
                samples[sample_name].metrics['major_contamination'] = \
                    round(het_sites.shape[0] / sites_notna.shape[0], 4)

        return samples
