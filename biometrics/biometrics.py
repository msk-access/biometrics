#!/usr/bin/env python
import sys

import pandas as pd

from cli import get_args
from sample import Sample
from extract import Extract
from utils import standardize_sex_nomenclature


def run_extract(args, samples):
    extractor = Extract(args=args)
    extractor.extract(samples)

    return samples


def run_sexmismatch(args, samples):
    pass


def run_minor_contamination(args, samples):
    pass


def run_major_contamination(args, samples):
    pass


def run_genotyping(args, samples):
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
                sex=standardize_sex_nomenclature(titlefile.at[i, 'Sex']),
                db=args.db)

            sample.find_titlefile_alignment(args.bam_basedir)
            samples.append(sample)

    return samples


def get_samples_list(args):
    samples = []

    for i, bam in enumerate(args.sample_bam):

        sex = standardize_sex_nomenclature(
            args.sample_sex[i] if args.sample_sex is not None else None)
        name = args.sample_name[i] if args.sample_name is not None else None
        patient = args.sample_patient[i] \
            if args.sample_patient is not None else None
        sample_type = args.sample_type[i] \
            if args.sample_type is not None else None

        sample = Sample(
            alignment_file=bam, patient=patient, name=name,
            sample_type=sample_type, sex=sex, db=args.database)

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
