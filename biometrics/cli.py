"""Console script for biometrics."""
import argparse

from .utils import exit_error


def add_common_args(parser):

    parser.add_argument(
        '-t', '--titlefile', action="append",
        help='''Path to title file. Can specify more than once.''')
    parser.add_argument(
        '-sb', '--sample-bam', help='''Space-delimited list of BAM files.''')
    parser.add_argument(
        '-st', '--sample-type', help='''Space-delimited list of sample
        types: Normal or Tumor. Must be in the same order as --sample-bam.''')
    parser.add_argument(
        '-ss', '--sample-sex', help='''Space-delimited list of sample
        sex (i.e. M or F). Must be in the same order as --sample-bam.''')
    parser.add_argument(
        '-sg', '--sample-group', help='''Space-delimited list of sample
        grouping information (e.g. P-D012F). Must be in the same order
        as --sample-bam.''')
    parser.add_argument(
        '-sn', '--sample-name', help='''Space-delimited list of sample names.
        If not specified, sample name is automatically figured out from the
        BAM file. Must be in the same order as --sample-bam.''')
    parser.add_argument(
        '-db', '--database',
        default=['/juno/work/access/production/data/bams'],
        help='''Directory to store the intermediate files after
        running the extraction step.''')
    parser.add_argument(
        '--bam-basedir', default='/juno/work/access/production/data/bams',
        help='''Base directory where BAM files are located if you specified
        title files.''')

    return parser


def add_outdir(parser):
    parser.add_argument(
        '-o', '--outdir', default='.',
        help='''Output directory for results.''')

    return parser


def check_args(args):

    if not args.titlefile and not args.sample_bam:
        exit_error('You must specify either --titlefile or --sample-bam')


def get_args():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''Various tools for fingerprinting samples from BAM files.
        Sample information to each sub command is supplied via title file(s)
        and/or as individual samples. For title files, the sample annotation
        and BAM location is automatically figure out.''')
    subparsers = parser.add_subparsers(help='', dest="")

    parser_extract = subparsers.add_parser(
        'extract',
        help='''Extract genotype info from one or more samples.
        The output from this step is required for the rest of the
        fingerprinting tools. However, you do not need to run this step
        manually since it will run automatically if the necessary files
        are missing.''')
    parser_extract = add_common_args(parser_extract)

    parser_sexmismatch = subparsers.add_parser(
        'sexmismatch', help='Check for sex mismatches.')
    parser_sexmismatch = add_common_args(parser_sexmismatch)
    parser_sexmismatch = add_outdir(parser_sexmismatch)

    parser_minor = subparsers.add_parser(
        'minor', help='Check for minor contamination.')
    parser_minor = add_common_args(parser_minor)
    parser_minor = add_outdir(parser_minor)

    parser_major = subparsers.add_parser(
        'major', help='Check for major contamination.')
    parser_major = add_common_args(parser_major)
    parser_major = add_outdir(parser_major)

    parser_genotype = subparsers.add_parser(
        'genotype', help='Genotype a set of samples.')
    parser_genotype = add_common_args(parser_genotype)
    parser_genotype = add_outdir(parser_genotype)

    args = parser.parse_args()

    check_args(args)

    return args
