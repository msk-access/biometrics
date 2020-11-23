import os
import glob

import pandas as pd

# from biometrics.sample import Sample
# from biometrics.extract import Extract
# from biometrics.genotype import Genotyper
# from biometrics.minor_contamination import MinorContamination
# from biometrics.major_contamination import MajorContamination
# from biometrics.utils import standardize_sex_nomenclature, exit_error

from sample import Sample
from extract import Extract
from genotype import Genotyper
from minor_contamination import MinorContamination
from major_contamination import MajorContamination
from utils import standardize_sex_nomenclature, exit_error


def load_extra_database_samples(args, existing_samples):

    samples = {}

    if args.no_db_comparison:
        return samples

    for pickle_file in glob.glob(os.path.join(args.database, '*pk')):

        sample_name = os.path.basename(pickle_file).replace('.pk', '')

        if sample_name in existing_samples:
            continue

        sample = Sample(db=args.database, is_in_db=False)
        sample.load_from_file(extraction_file=pickle_file)

        samples[sample.name] = sample

    return samples


def run_extract(args, samples):
    extractor = Extract(args=args)
    samples = extractor.extract(samples)
    existing_samples = set([i for i in samples.keys()])

    samples.update(load_extra_database_samples(args, existing_samples))

    return samples


def run_sexmismatch(args, samples):
    pass


def run_minor_contamination(args, samples):
    minor_contamination = MinorContamination(args)
    samples = minor_contamination.estimate(samples)

    return samples


def run_major_contamination(args, samples):
    major_contamination = MajorContamination(args)
    samples = major_contamination.estimate(samples)

    return samples


def run_genotyping(args, samples):
    genotyper = Genotyper(args)
    genotyper.genotype(samples)

    return samples


def get_samples_from_input(args):
    samples = {}

    for fpath in args.input:

        input = pd.read_csv(fpath, sep=',')

        for i in input.index:

            alignment_file = input.at[i, 'alignment_file']

            if not os.path.exists(alignment_file):
                exit_error('Alignment file does not exist: {}.'.format(
                    alignment_file))

            sample = Sample(
                alignment_file=alignment_file,
                group=input.at[i, 'group'],
                name=input.at[i, 'sample_name'],
                sample_type=input.at[i, 'type'],
                sex=standardize_sex_nomenclature(input.at[i, 'sex']),
                db=args.db)

            samples[sample.name] = sample

    return samples


def get_samples_list(args):
    samples = {}

    for i, bam in enumerate(args.sample_bam):

        sex = standardize_sex_nomenclature(
            args.sample_sex[i] if args.sample_sex is not None else None)
        name = args.sample_name[i] if args.sample_name is not None else None
        group = args.sample_group[i] \
            if args.sample_group is not None else None
        sample_type = args.sample_type[i] \
            if args.sample_type is not None else None

        sample = Sample(
            alignment_file=bam, group=group, name=name,
            sample_type=sample_type, sex=sex, db=args.database)

        samples[sample.name] = sample

    return samples


def get_samples(args):

    samples = {}

    if args.input:
        samples.update(get_samples_from_input(args))

    if args.sample_bam:
        samples.update(get_samples_list(args))

    return samples


def run_biometrics(args):

    samples = get_samples(args)
    samples = run_extract(args, samples)

    if args.subparser_name == 'sexmismatch':
        run_sexmismatch(args, samples)
    elif args.subparser_name == 'minor':
        run_minor_contamination(args, samples)
    elif args.subparser_name == 'major':
        run_major_contamination(args, samples)
    elif args.subparser_name == 'genotype':
        run_genotyping(args, samples)
