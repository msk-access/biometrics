#!/usr/bin/env python
import sys

from cli import get_args
from sample import Sample
import pandas as pd


def run_extract(args):
    pass


def run_sexmismatch(args):
    pass


def run_minor_contamination(args):
    pass


def run_major_contamination(args):
    pass


def run_genotyping(args):
    pass


def get_samples_from_titlefile(args):
    samples = []

    for fpath in args.titlefile:

        titlefile = pd.read_csv(fpath, sep='\t')

        for i in titlefile.index:
            sample = Sample(
                patient=titlefile.at[i, 'Patient_ID'],
                name=titlefile.at[i, 'Patient_ID'],
                sample_type=titlefile.at[i, 'Sample_type'],
                sex=titlefile.at[i, 'Sex'])

            sample.find_titlefile_alignment(args.bam_basedir)
            samples.append(sample)

    return samples


def get_samples_list(args):
    samples = []

    for i, bam in enumerate(args.sample_bams):

        patient = args.sample_patient[i] \
            if args.sample_patient[i] is not None else None
        sex = args.sample_sex[i] \
            if args.sample_sex[i] is not None else None
        sample_type = args.sample_type[i] \
            if args.sample_type[i] is not None else None
        name = args.sample_name[i] \
            if args.sample_name[i] is not None else None

        sample = Sample(
            patient=patient, name=name, sample_type=sample_type(), sex=sex)

        samples.append(sample)

    return samples


def get_samples(args):

    samples = []

    if args.titlefile:
        samples += get_samples_from_titlefile(args)

    if args.sample_bam:
        samples += get_samples_list(args)

    return samples


def main():

    args = get_args()

    samples = get_samples(args)

    if args.subparser_name == 'extract':
        run_extract(args)
    elif args.subparser_name == 'sexmismatch':
        run_sexmismatch(args)
    elif args.subparser_name == 'minor':
        run_minor_contamination(args)
    elif args.subparser_name == 'major':
        run_major_contamination(args)
    elif args.subparser_name == 'genotype':
        run_genotyping(args)


if __name__ == "__main__":
    sys.exit(main())
