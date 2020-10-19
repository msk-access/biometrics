"""Console script for biometrics."""
import sys
import argparse


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

    return parser


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
        help='Extract genotype info from one or more samples.')

    parser_sexmismatch = subparsers.add_parser(
        'sexmismatch',
        help='Check for sex mismatches.')

    parser_minor = subparsers.add_parser(
        'minor',
        help='Check for minor contamination.')

    parser_major = subparsers.add_parser(
        'major',
        help='Check for major contamination.')

    parser_genotype = subparsers.add_parser(
        'genotype',
        help='Genotype a set of samples.')

    # project_parser = add_common_args(project_parser)

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
