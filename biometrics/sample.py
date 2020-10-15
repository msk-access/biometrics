

class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, alignment_file, name, patient=None, sex=None):
        self.alignment_file = alignment_file
        self.name = name
        self.sex = sex
        self.patient = patient
