import pandas as pd
import numpy as np


class SexMismatch:
    """
    Class to detect sex mismatch
    """

    def __init__(self, threshold):
        self.threshold = threshold

    def predict_sex(self, sample):

        if sample.region_counts is None:
            return np.nan

        total_count = sample.region_counts['count'].sum()

        predicted_sex = 'M' if total_count > self.threshold else 'F'

        return predicted_sex

    def detect_mismatch(self, samples):

        results = []

        for i, sample_name in enumerate(samples):

            sample = samples[sample_name]

            predicted_sex = self.predict_sex(sample)

            results.append({
                'sample': sample_name,
                'expected_sex': sample.sample_sex,
                'predicted_sex': predicted_sex
            })

        results = pd.DataFrame(results)

        results['sex_mismatch'] = \
            (results['expected_sex'] != results['predicted_sex']).astype(str)
        results.loc[pd.isna(results['predicted_sex']), 'sex_mismatch'] = np.nan

        return results
