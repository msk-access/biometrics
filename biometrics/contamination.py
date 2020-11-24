import pandas as pd
import numpy as np


class BaseContamination:

    def __init__(self, args):
        pass

    def to_dataframe(self, samples, contamination_type):

        data = pd.DataFrame(
            columns=['sample', 'sample_group', 'sample_sex', 'sample_type',
                     contamination_type])

        for sample_name, sample in samples.items():

            row = {
                'sample': sample_name,
                'sample_group': sample.group,
                'sample_sex': sample.sex,
                'sample_type': sample.sample_type,
                contamination_type: sample.metrics.get(contamination_type, np.nan)
            }

            data = data.append(row, ignore_index=True)

        return data
