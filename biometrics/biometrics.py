import os
import glob

import pandas as pd

from biometrics.sample import Sample
from biometrics.extract import Extract
from biometrics.genotype import Genotyper
from biometrics.minor_contamination import MinorContamination
from biometrics.major_contamination import MajorContamination
from biometrics.sex_mismatch import SexMismatch
from biometrics.utils import standardize_sex_nomenclature, exit_error


def write_to_file(args, data, basename):

    outdir = os.path.abspath(args.outdir)

    outpath = os.path.join(outdir, basename + '.csv')
    data.to_csv(outpath, index=False)

    if args.json:
        outpath = os.path.join(outdir, basename + '.json')
        data.to_json(outpath)


def load_database_samples(database, existing_samples):

    samples = {}

    for pickle_file in glob.glob(os.path.join(database, '*pk')):

        sample_name = os.path.basename(pickle_file).replace('.pk', '')

        if sample_name in existing_samples:
            continue

        sample = Sample(db=database, query_group=True)
        sample.load_from_file(extraction_file=pickle_file)

        samples[sample.name] = sample

    return samples


def run_extract(args, samples):
    extractor = Extract(args=args)
    samples = extractor.extract(samples)

    return samples


def run_sexmismatch(args, samples):
    sex_mismatch = SexMismatch(50)

    results = sex_mismatch.detect_mismatch(samples)
    write_to_file(args, results, 'sex_mismatch')


def run_minor_contamination(args, samples):
    minor_contamination = MinorContamination(threshold=args.minor_threshold)
    samples = minor_contamination.estimate(samples)

    data = minor_contamination.to_dataframe(samples)
    write_to_file(args, data, 'minor_contamination')

    if args.plot:
        minor_contamination.plot(data, args.outdir)

    return samples


def run_major_contamination(args, samples):
    major_contamination = MajorContamination(threshold=args.major_threshold)
    samples = major_contamination.estimate(samples)

    data = major_contamination.to_dataframe(samples)
    write_to_file(args, data, 'major_contamination')

    if args.plot:
        major_contamination.plot(data, args.outdir)

    return samples


def run_genotyping(args, samples):
    genotyper = Genotyper(args.no_db_compare, args.discordance_threshold)
    data = genotyper.genotype(samples)

    write_to_file(args, data, 'genotype_comparison')

    if args.plot:
        genotyper.plot(data, args.outdir)

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
                db=args.database)

            samples[sample.name] = sample

    return samples


def get_samples_from_bam(args):
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


def get_samples_from_name(args):
    samples = {}

    for i, name in enumerate(args.sample_name):

        extraction_file = os.path.join(args.database, name + '.pk')

        if not os.path.exists(extraction_file):
            exit_error(
                'Could not find: {}. Please rerun the extraction step.'.format(
                    extraction_file))

        sample = Sample(query_group=True)
        sample.load_from_file(extraction_file)

        samples[sample.name] = sample

    return samples


def get_samples(args, extraction_mode=False):

    samples = {}

    if args.input:
        samples.update(get_samples_from_input(args))

    if extraction_mode:
        if args.sample_bam:
            samples.update(get_samples_from_bam(args))
    else:
        if args.sample_name:
            samples.update(get_samples_from_name(args))

        for sample_name in samples.keys():
            extration_file = os.path.join(args.database, sample_name + '.pk')
            samples[sample_name].load_from_file(extration_file)

        existing_samples = set([i for i in samples.keys()])

        if not args.no_db_compare:
            samples.update(load_database_samples(
                args.database, existing_samples))

    return samples


def create_outdir(outdir):
    os.makedirs(outdir, exist_ok=True)


def run_biometrics(args):

    extraction_mode = args.subparser_name == 'extract'

    samples = get_samples(args, extraction_mode=extraction_mode)

    if extraction_mode:
        create_outdir(args.database)
        run_extract(args, samples)
    elif args.subparser_name == 'sexmismatch':
        create_outdir(args.outdir)
        run_sexmismatch(args, samples)
    elif args.subparser_name == 'minor':
        create_outdir(args.outdir)
        run_minor_contamination(args, samples)
    elif args.subparser_name == 'major':
        create_outdir(args.outdir)
        run_major_contamination(args, samples)
    elif args.subparser_name == 'genotype':
        create_outdir(args.outdir)
        run_genotyping(args, samples)
