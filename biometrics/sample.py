import pickle
import os

import pandas as pd


class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, sample_name=None, sample_bam=None, sample_group=None,
                 sample_sex=None, sample_type=None, db=None, query_group=False):
        self.sample_bam = sample_bam
        self.sample_name = sample_name
        self.sample_sex = sample_sex
        self.sample_type = sample_type

        if sample_group is None:
            self.sample_group = sample_name
        else:
            self.sample_group = sample_group

        self.pileup = None
        self.region_counts = None
        self.extraction_file = None
        self.query_group = query_group
        self.metrics = {}

        if self.sample_name is not None:
            if db is not None:
                self.extraction_file = os.path.join(db, self.sample_name + '.pk')
            else:
                self.extraction_file = self.sample_name + '.pk'

    def save_to_file(self):

        pileup_data = self.pileup.to_dict("records")

        if self.region_counts is not None:
            region_counts = self.region_counts.to_dict('records')
        else:
            region_counts = None

        sample_data = {
            'sample_bam': self.sample_bam,
            'sample_name': self.sample_name,
            'sample_sex': self.sample_sex,
            'sample_group': self.sample_group,
            'sample_type': self.sample_type,
            'pileup_data': pileup_data,
            'region_counts': region_counts
        }

        pickle.dump(sample_data, open(self.extraction_file, "wb"))

    def load_from_file(self, extraction_file=None):

        if extraction_file is not None:
            self.extraction_file = extraction_file

        assert self.extraction_file is not None, 'Extraction file path is None.'
        assert os.path.exists(self.extraction_file), 'Extraction file does not exist: {}'.format(
            self.extraction_file)

        sample_data = pickle.load(open(self.extraction_file, "rb"))

        region_counts = None
        if sample_data.get('region_counts') is not None:
            region_counts = pd.DataFrame(
                sample_data['region_counts'], dtype=object)

        self.pileup = pd.DataFrame(sample_data['pileup_data'])
        self.sample_bam = sample_data['sample_bam']
        self.sample_name = sample_data['sample_name']
        self.sample_sex = sample_data['sample_sex']
        self.sample_group = sample_data['sample_group'] if sample_data['sample_group'] is not None else sample_data['sample_name']
        self.sample_type = sample_data['sample_type']
        self.region_counts = region_counts
