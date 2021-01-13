import sys


def exit_error(msg):
    print("ERROR: {}".format(msg))
    sys.exit(1)


def standardize_sex_nomenclature(val):

    if val is None:
        return None

    # Potential inputs
    female = ['female', 'f', 'Female', 'F']
    male = ['Male', 'M', 'male', 'm']

    if val in female:
        return 'F'
    elif val in male:
        return 'M'

    return None
