import glob
import os

from utils import exit_error


class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, name, alignment_file=None, patient=None, sex=None, sample_type=None, db=None):
        self.alignment_file = alignment_file
        self.name = name
        self.sex = sex
        self.patient = patient
        self.sample_type = sample_type
        self.pileup = None

        if db is not None:
            self.extraction_file = os.path.join(db, self.name + '.pk')
        else:
            self.extraction_file = os.path.join(os.getcwd(), self.name, '.pk')

    def find_titlefile_alignment(self, basedir, bam_type=None):

        bam_basedir = os.path.join(basedir, self.patient, self.name, 'latest')

        bam = glob.glob(
            os.path.join(
                bam_basedir,
                '*_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX.bam'))

        if not bam:
            exit_error('Could not find BAM file for {}.'.format(self.name))
        elif len(bam) > 1:
            exit_error(
                'Found more than one BAM file for {}.'.format(self.name))

        self.alignment_file = bam[0]
