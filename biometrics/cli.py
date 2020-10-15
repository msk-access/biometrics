"""Console script for biometrics."""
import sys
import argparse


def get_args():

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='''
            Various tools for fingerprinting samples.
            ''')
    subparsers = parser.add_subparsers(
        help='Creates BAM/BAI file links.',
        dest="subparser_name")

    project_parser = subparsers.add_parser(
        'extract',
        help='Extract genotype info from one or more samples.')
    project_parser.add_argument(
        '-t', '--titlefile', action="append",
        help='''
            Title file.
            Can specify more than once. Optional.''')
    # project_parser = add_common_args(project_parser)

    args = parser.parse_args()

    return args

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
