import os
import glob

import pandas as pd

from biometrics.sample import Sample
from biometrics.extract import Extract
from biometrics.genotype import Genotyper
from biometrics.cluster import Cluster
from biometrics.minor_contamination import MinorContamination
from biometrics.major_contamination import MajorContamination
from biometrics.sex_mismatch import SexMismatch
from biometrics.utils import standardize_sex_nomenclature, get_logger

logger = get_logger()


def write_to_file(args, data, basename):
    """
    Generic function to save output to a file.
    """

    outdir = os.path.abspath(args.outdir)

    outpath = os.path.join(outdir, basename + '.csv')
    data.to_csv(outpath, index=False)

    if args.json:
        outpath = os.path.join(outdir, basename + '.json')
        data.to_json(outpath)


def run_extract(args, samples):
    """
    Extract the pileup and region information from the samples. Then
    save to the database.
    """

    extractor = Extract(args=args)
    samples = extractor.extract(samples)

    return samples


def run_sexmismatch(args, samples):
    """
    Find and sex mismatches and save the output
    """

    sex_mismatch = SexMismatch(args.coverage_threshold)

    results = sex_mismatch.detect_mismatch(samples)

    basename = 'sex_mismatch'
    if args.prefix:
        basename = args.prefix + '_' + basename

    write_to_file(args, results, basename)


def run_minor_contamination(args, samples):
    """
    Compute minor contamination and save the output and figure
    """

    minor_contamination = MinorContamination(threshold=args.minor_threshold)
    samples = minor_contamination.estimate(samples)

    data = minor_contamination.to_dataframe(samples)

    basename = 'minor_contamination'
    if args.prefix:
        basename = args.prefix + '_' + basename

    write_to_file(args, data, basename)

    if args.plot:
        if len(samples) > 1000:
            logger.warning('Turning off plotting functionality. You are trying to plot more than 1000 samples, which is too cumbersome.')
        else:
            minor_contamination.plot(samples, args.outdir)

    return samples


def run_major_contamination(args, samples):
    """
    Compute major contamination and save the output and figure.
    """

    major_contamination = MajorContamination(threshold=args.major_threshold)
    samples = major_contamination.estimate(samples)

    data = major_contamination.to_dataframe(samples)

    basename = 'major_contamination'
    if args.prefix:
        basename = args.prefix + '_' + basename

    write_to_file(args, data, basename)

    if args.plot:
        if len(samples) > 1000:
            logger.warning('Turning off plotting functionality. You are trying to plot more than 1000 samples, which is too cumbersome.')
        else:
            major_contamination.plot(samples, args.outdir)

    return samples


def run_genotyping(args, samples):
    """
    Run the genotyper and save the output and figure.
    """

    genotyper = Genotyper(
        no_db_compare=args.no_db_compare,
        discordance_threshold=args.discordance_threshold,
        threads=args.threads,
        zmin=args.zmin,
        zmax=args.zmax)
    cluster_handler = Cluster(args.discordance_threshold)
    comparisons = genotyper.compare_samples(samples)

    # save genotyping output

    basename = 'genotype_comparison'
    if args.prefix:
        basename = args.prefix + '_' + basename

    write_to_file(args, comparisons, basename)

    # cluster just the input samples

    samples_input = dict(filter(
        lambda x: not x[1].query_group, samples.items()))

    samples_names = [i.sample_name for i in samples_input.values()]
    comparisons_input = comparisons[
        (comparisons['ReferenceSample'].isin(samples_names)) &
        (comparisons['QuerySample'].isin(samples_names))].copy()

    logger.info('Clustering input samples...')
    clusters = cluster_handler.cluster(comparisons_input)

    if clusters is not None:
        basename = 'genotype_clusters_input'
        if args.prefix:
            basename = args.prefix + '_' + basename

        write_to_file(args, clusters, basename)

    # cluster all the samples

    if not args.no_db_compare:

        are_there_db_samples = len(samples_input) != len(samples)

        if not are_there_db_samples:
            logger.warning(
                'The set of database and input samples are the same. Will only cluster the samples once.')
        else:
            logger.info('Clustering input and database samples...')
            clusters = cluster_handler.cluster(comparisons)

            if clusters is not None:
                basename = 'genotype_clusters_database'
                if args.prefix:
                    basename = args.prefix + '_' + basename

                write_to_file(args, clusters, basename)

    # save plots

    if args.plot:
        if len(samples) > 1000:
            logger.warning('Turning off plotting functionality. You are trying to plot more than 1000 samples, which is too cumbersome.')
        else:
            genotyper.plot(comparisons, args.outdir)

    return samples


def run_cluster(args):

    comparisons = []
    for input in args.input:
        comparisons.append(
            pd.read_csv(input)
        )
    comparisons = pd.concat(comparisons)
    comparisons = comparisons.drop_duplicates(['ReferenceSample', 'QuerySample'])

    cluster_handler = Cluster(args.discordance_threshold)

    logger.info('Clustering input samples...')
    clusters = cluster_handler.cluster(comparisons)

    if clusters is not None:
        clusters.to_csv(args.output, index=False)


def load_input_sample_from_db(sample_name, database):
    """
    Loads any the given (that the user specified via the CLI) from the
    database.
    """

    extraction_file = os.path.join(database, sample_name + '.pk')

    assert os.path.exists(extraction_file), 'Could not find: {}. Please rerun the extraction step.'.format(
        extraction_file)

    sample = Sample(query_group=False)
    sample.load_from_file(extraction_file)

    return sample


def load_database_samples(database, existing_samples):
    """
    Loads any samples that are already present in the database AND
    which were not specified as input via the CLI.
    """

    samples = {}

    for pickle_file in glob.glob(os.path.join(database, '*pk')):

        sample_name = os.path.basename(pickle_file).replace('.pk', '')

        if sample_name in existing_samples:
            continue

        sample = Sample(db=database, query_group=True)
        sample.load_from_file(extraction_file=pickle_file)

        samples[sample.sample_name] = sample

    return samples


def get_samples_from_input(inputs, database, extraction_mode):
    """
    Parse the sample information from the user-supplied CSV file.
    """

    if type(inputs) != list:
        inputs = [inputs]

    samples = {}

    for fpath in inputs:

        input = pd.read_csv(fpath, sep=',')

        # check if some required columns are present

        assert 'sample_bam' in input.columns, 'Input file does not have the \'sample_bam\' column.'
        assert 'sample_name' in input.columns, 'Input does not have \'sample_name\' column.'

        input = input.to_dict(orient='records')

        for row in input:

            if not extraction_mode:
                # if not running extract tool, then just need to get
                # the sample name

                sample_name = row['sample_name']

                sample = load_input_sample_from_db(sample_name, database)
                samples[sample.sample_name] = sample

                continue

            # parse in the input

            sample = Sample(
                sample_name=row['sample_name'],
                sample_bam=row['sample_bam'],
                sample_group=row.get('sample_group'),
                sample_type=row.get('sample_type'),
                sample_sex=standardize_sex_nomenclature(row.get('sample_sex')),
                db=database)

            samples[sample.sample_name] = sample

    return samples


def get_samples_from_bam(args):
    """
    Parse the sample information the user supplied via the CLI.
    """

    samples = {}

    for i, sample_bam in enumerate(args.sample_bam):

        sample_sex = standardize_sex_nomenclature(
            args.sample_sex[i] if args.sample_sex is not None else None)
        sample_name = args.sample_name[i] if args.sample_name is not None else None
        sample_group = args.sample_group[i] \
            if args.sample_group is not None else args.sample_name[i]
        sample_type = args.sample_type[i] \
            if args.sample_type is not None else None

        sample = Sample(
            sample_bam=sample_bam, sample_group=sample_group,
            sample_name=sample_name, sample_type=sample_type,
            sample_sex=sample_sex, db=args.database)

        samples[sample.sample_name] = sample

    return samples


def get_samples_from_name(sample_names, database):
    """
    Parse the sample information the user supplied via the CLI.
    """

    if type(sample_names) != list:
        sample_names = [sample_names]

    samples = {}

    for i, sample_name in enumerate(sample_names):
        sample = load_input_sample_from_db(sample_name, database)
        samples[sample.sample_name] = sample

    return samples


def get_samples(args, extraction_mode=False):
    """
    Parse the sample information the user supplied via the CLI.
    """

    samples = {}

    if extraction_mode:
        if args.input:
            samples.update(get_samples_from_input(
                args.input, args.database, extraction_mode))

        if args.sample_bam:
            samples.update(get_samples_from_bam(args))
    else:

        for input in args.input:
            if input.endswith('.pk'):
                sample = Sample(db=args.database, query_group=False)
                sample.load_from_file(extraction_file=input)
                samples[sample.sample_name] = sample
            elif input.endswith('.csv') or input.endswith('.txt'):
                samples.update(get_samples_from_input(
                    input, args.database, extraction_mode))
            else:
                samples.update(get_samples_from_name(input, args.database))

        existing_samples = set([i for i in samples.keys()])

        if not args.no_db_compare:
            samples.update(load_database_samples(
                args.database, existing_samples))

    return samples


def create_outdir(outdir):
    if outdir is None:
        return

    os.makedirs(outdir, exist_ok=True)


def run_biometrics(args):
    """
    Decide what tool to run based in CLI input.
    """

    if args.subparser_name == 'cluster':
        run_cluster(args)
        return

    extraction_mode = args.subparser_name == 'extract'

    samples = get_samples(args, extraction_mode=extraction_mode)

    # if not extraction_mode and args.plot:

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
