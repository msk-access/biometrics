#!/usr/bin/env python
import sys
import os

import pandas as pd

from cli import get_args
from sample import Sample
from extract import Extract
from genotype import Genotyper
from minor_contamination import MinorContamination
from major_contamination import MajorContamination
from utils import standardize_sex_nomenclature, exit_error


def run_extract(args, samples):
    extractor = Extract(args=args)
    extractor.extract(samples)

    return samples


def run_sexmismatch(args, samples):
    pass


def run_minor_contamination(args, samples):
    minor_contamination = MinorContamination(args)
    samples = minor_contamination.estimate(samples)


def run_major_contamination(args, samples):
    major_contamination = MajorContamination(args)
    samples = major_contamination.estimate(samples)


def run_genotyping(args, samples):
    genotyper = Genotyper(args)
    genotyper.genotype(samples)


def get_samples_from_input(args):
    samples = []

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

            samples.append(sample)

    return samples


def get_samples_list(args):
    samples = []

    for i, bam in enumerate(args.sample_bam):

        sex = standardize_sex_nomenclature(
            args.sample_sex[i] if args.sample_sex is not None else None)
        name = args.sample_name[i] if args.sample_name is not None else None
        group = args.sample_group[i] \
            if args.sample_patient is not None else None
        sample_type = args.sample_type[i] \
            if args.sample_type is not None else None

        sample = Sample(
            alignment_file=bam, group=group, name=name,
            sample_type=sample_type, sex=sex, db=args.database)

        samples.append(sample)

    return samples


def get_samples(args):

    samples = []

    if args.input:
        samples += get_samples_from_input(args)

    if args.sample_bam:
        samples += get_samples_list(args)

    return samples


def main():

    args = get_args()

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


if __name__ == "__main__":
    sys.exit(main())
