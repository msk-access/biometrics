"""Top-level package for biometrics."""

import os

__author__ = """Charlie Murphy"""
__email__ = 'murphyc4@mskcc.org'

# version info

library_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(library_path, "VERSION"), "r") as fh:
    __version__ = fh.read().strip()
