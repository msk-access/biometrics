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
        self.is_in_db = is_in_db
        self.metrics = {}

        if db is not None:
            self.extraction_file = os.path.join(db, self.name + '.pk')
        else:
            self.extraction_file = os.path.join(os.getcwd(), self.name, '.pk')

    def save_to_file(self):

        pileup_data = self.pileup.to_dict()
        sample_data = {
            'alignment_file': self.alignment_file,
            'name': self.name,
            'sex': self.sex,
            'group': self.group,
            'sample_type': self.sample_type,
            'pileup_data': pileup_data
        }

        pickle.dump(sample_data, open(self.extraction_file, "wb"))

    def load_from_file(self, sample):

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
