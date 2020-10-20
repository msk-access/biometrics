import os


class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, alignment_file=None, name=None, patient=None, sex=None, sample_type=None):
        self.alignment_file = None
        self.name = name
        self.sex = sex
        self.patient = patient
        self.sample_type = sample_type

    def find_titlefile_alignment(self, basedir, bam_type=None):

        bam_basedir = os.path.join(basedir, self.patient, self.name, 'latest')
