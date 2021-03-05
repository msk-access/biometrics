import logging


def get_logger(debug=False):
    FORMAT = '%(levelname)s - %(asctime)-15s: %(message)s'
    logging.basicConfig(format=FORMAT)
    logger = logging.getLogger("biometrics")

    if debug:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    return logger


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
