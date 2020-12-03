import pickle
import os

import pandas as pd
from utils import exit_error


class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, name=None, alignment_file=None, group=None, sex=None, sample_type=None, db=None, is_in_db=False):
        self.alignment_file = alignment_file
        self.name = name
        self.sex = sex
        self.group = group
        self.sample_type = sample_type
        self.pileup = None
        self.region_counts = None
        self.extraction_file = None
        self.is_in_db = is_in_db
        self.metrics = {}

        if db is not None and self.name is not None:
            self.extraction_file = os.path.join(db, self.name + '.pk')

    def save_to_file(self):

        pileup_data = self.pileup.to_dict("records")
        region_counts = self.region_counts.to_dict('records')
        sample_data = {
            'alignment_file': self.alignment_file,
            'name': self.name,
            'sex': self.sex,
            'group': self.group,
            'sample_type': self.sample_type,
            'pileup_data': pileup_data,
            'region_counts': region_counts
        }

        pickle.dump(sample_data, open(self.extraction_file, "wb"))

    def load_from_file(self, extraction_file=None):

        if extraction_file is not None:
            self.extraction_file = extraction_file

        if self.extraction_file is None:
            exit_error('Extraction file path is None.')

        if not os.path.exists(self.extraction_file):
            exit_error('Extraction file does not exist: {}'.format(
                self.extraction_file))

        sample_data = pickle.load(open(self.extraction_file, "rb"))

        self.pileup = pd.DataFrame(sample_data['pileup_data'])
        self.alignment_file = sample_data['alignment_file']
        self.name = sample_data['name']
        self.sex = sample_data['sex']
        self.group = sample_data['group']
        self.sample_type = sample_data['sample_type']
        self.region_counts = pd.DataFrame(
            sample_data['region_counts'], dtype=object)
