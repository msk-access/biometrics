#!/usr/bin/env python
"""Console script for biometrics."""
import sys
import argparse

from utils import exit_error
# from biometrics.biometrics import run_biometrics
from biometrics import run_biometrics


def add_common_args(parser):

    parser.add_argument(
        '-i', '--input', action="append", required=False,
        help='''Path to file containing sample information (one per line). For example: sample_name,alignment_file,type,sex,group''')
    parser.add_argument(
        '-sb', '--sample-bam', nargs="+", required=False,
        help='''Space-delimited list of BAM files.''')
    parser.add_argument(
        '-st', '--sample-type', nargs="+", required=False,
        help='''Space-delimited list of sample types: Normal or Tumor.
        Must be in the same order as --sample-bam.''')
    parser.add_argument(
        '-ss', '--sample-sex', nargs="+", required=False,
        help='''Space-delimited list of sample sex (i.e. M or F). Must be
        in the same order as --sample-bam.''')
    parser.add_argument(
        '-sg', '--sample-group', nargs="+", required=False,
        help='''Space-delimited list of sample group information
        (e.g. sample patient ID). Must be in the same order as --sample-bam.''')
    parser.add_argument(
        '-sn', '--sample-name', nargs="+", required=False,
        help='''Space-delimited list of sample names. If not specified,
        sample name is automatically figured out from the BAM file. Must
        be in the same order as --sample-bam.''')
    parser.add_argument(
        '--vcf', required=True,
        help='''VCF file containing the sites to be queried.''')
    parser.add_argument(
        '--bed', required=False,
        help='''BED file containing the intervals to be queried.''')
    parser.add_argument(
        '-db', '--database',
        help='''Directory to store the intermediate files after
        running the extraction step.''')
    parser.add_argument(
        '-ov', '--overwrite', action='store_true',
        help='''Overwrite any existing extraction results.''')
    parser.add_argument(
        '-nc', '--no-db-compare', action='store_true',
        help='''Do not compare the sample(s) you provided to all samples in the database.''')
    parser.add_argument(
        '--fafile', required=True,
        help='''Path to reference fasta file.''')
    parser.add_argument(
        '-q', '--min-mapping-quality', default=1, type=int,
        help='''Minimum mapping quality of reads to be used for pileup.''')
    parser.add_argument(
        '-Q', '--min-base-quality', default=1, type=int,
        help='''Minimum base quality of reads to be used for pileup.''')
    parser.add_argument(
        '-mc', '--min-coverage', default=10, type=int,
        help='''Minimum coverage to count a site.''')

    return parser


def add_common_tool_args(parser):
    parser.add_argument(
        '-o', '--outdir', default='.',
        help='''Output directory for results.''')
    parser.add_argument(
        '-p', '--plot', action='store_true',
        help='''Also output plots of the data.''')
    parser.add_argument(
        '-j', '--json', action='store_true',
        help='''Also output data in JSON format.''')

    return parser


def check_arg_equal_len(vals1, vals2, name):
    if vals2 is not None and len(vals1) != len(vals2):
        exit_error(
            '{} does not have the same number of items as --sample-bam'.format(
                name))


def check_args(args):

    if not args.input and not args.sample_bam:
        exit_error('You must specify either --input or --sample-bam')

    if args.sample_bam:
        check_arg_equal_len(args.sample_bam, args.sample_name, '--sample-name')
        check_arg_equal_len(args.sample_bam, args.sample_type, '--sample-type')
        check_arg_equal_len(
            args.sample_bam, args.sample_group, '--sample-group')
        check_arg_equal_len(args.sample_bam, args.sample_sex, '--sample-sex')


def get_args():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Various tools for fingerprinting samples from BAM files.
        Sample information to each sub command is supplied via input file(s)
        and/or as individual samples.''')
    subparsers = parser.add_subparsers(help='', dest="subparser_name")

    # extract parser

    parser_extract = subparsers.add_parser(
        'extract',
        help='''Intermediate step to extract genotype info from one or more
        samples. The output from this step is required for the rest of the
        fingerprinting tools. However, you do not need to run this step
        manually since it will run automatically if the necessary files
        are missing.''')
    parser_extract = add_common_args(parser_extract)

    # sex mismatch parser

    parser_sexmismatch = subparsers.add_parser(
        'sexmismatch', help='Check for sex mismatches.')
    parser_sexmismatch = add_common_args(parser_sexmismatch)
    parser_sexmismatch = add_common_tool_args(parser_sexmismatch)

    # minor contamination parser

    parser_minor = subparsers.add_parser(
        'minor', help='Check for minor contamination.')
    parser_minor = add_common_args(parser_minor)
    parser_minor = add_common_tool_args(parser_minor)

    # major contamination parser

    parser_major = subparsers.add_parser(
        'major', help='Check for major contamination.')
    parser_major = add_common_args(parser_major)
    parser_major = add_common_tool_args(parser_major)

    # genotyping parser

    parser_genotype = subparsers.add_parser(
        'genotype', help='Genotype a set of samples.')
    parser_genotype = add_common_args(parser_genotype)
    parser_genotype = add_common_tool_args(parser_genotype)

    args = parser.parse_args()

    check_args(args)

    return args


def main():

    args = get_args()

    run_biometrics(args)


if __name__ == "__main__":
    sys.exit(main())
