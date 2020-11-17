import glob
import os

from utils import exit_error


class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, name, alignment_file=None, group=None, sex=None, sample_type=None, db=None):
        self.alignment_file = alignment_file
        self.name = name
        self.sex = sex
        self.group = group
        self.sample_type = sample_type
        self.pileup = None
        self.metrics = {}

        if db is not None:
            self.extraction_file = os.path.join(db, self.name + '.pk')
        else:
            self.extraction_file = os.path.join(os.getcwd(), self.name, '.pk')
