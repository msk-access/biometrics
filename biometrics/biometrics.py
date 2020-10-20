#!/usr/bin/env python
import sys

import pandas as pd
import vcf

from cli import get_args
from sample import Sample
from extract import Extract


def run_extract(args, samples, sites):
    extractor = Extract(db=args.db, sites=sites)
    extractor.extract(samples)


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
                sex=titlefile.at[i, 'Sex'],
                db=args.db)

            sample.find_titlefile_alignment(args.bam_basedir)
            samples.append(sample)

    return samples


def get_samples_list(args):
    samples = []

    for i, bam in enumerate(args.sample_bams):

        sex = args.sample_sex[i] if args.sample_sex is not None else None
        name = args.sample_name[i] if args.sample_name is not None else None
        patient = args.sample_patient[i] \
            if args.sample_patient is not None else None
        sample_type = args.sample_type[i] \
            if args.sample_type is not None else None

        sample = Sample(
            alignment_file=bam, patient=patient, name=name,
            sample_type=sample_type, sex=sex, db=args.db)

        samples.append(sample)

    return samples


def get_samples(args):

    samples = []

    if args.titlefile:
        samples += get_samples_from_titlefile(args)

    if args.sample_bam:
        samples += get_samples_list(args)

    return samples


def parse_vcf(vcf_file):
    sites = []

    for record in vcf.Reader(open(vcf_file, 'r')):
        sites.append({
            'chrom': record.CHROM,
            'pos': record.POS
        })

    return sites


def main():

    args = get_args()

    samples = get_samples(args)
    sites = parse_vcf(args.vcf)

    run_extract(args, samples, sites)

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
