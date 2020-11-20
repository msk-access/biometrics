#!/usr/bin/env python

"""Tests for `biometrics` package."""


import os
import argparse
from unittest import TestCase
from unittest import mock

from biometrics.biometrics import get_samples, run_minor_contamination, run_major_contamination
from biometrics.cli import get_args
from biometrics.extract import Extract

CUR_DIR = os.path.dirname(os.path.abspath(__file__))


class TestBiometrics(TestCase):
    """Tests for `biometrics` package."""

    @mock.patch(
        'argparse.ArgumentParser.parse_args',
        return_value=argparse.Namespace(
            input=None,
            sample_bam=[os.path.join(CUR_DIR, 'test_data/test_sample1_golden.bam')],
            sample_name=['test_sample1'],
            sample_type=['tumor'],
            sample_group=['test_sample1'],
            sample_sex=None,
            database=os.path.join(CUR_DIR, 'test_data/'),
            min_mapping_quality=1,
            min_base_quality=1,
            vcf=os.path.join(CUR_DIR, 'test_data/test.vcf'),
            fafile=os.path.join(CUR_DIR, 'test_data/ref.fasta'),
            overwrite=True,
            min_coverage=10))
    def setUp(self, mock_args):
        """Set up test fixtures, if any."""

        self.args = get_args()

    def test_load_vcf(self):
        """Test loading the VCF file."""

        extractor = Extract(self.args)

        self.assertGreater(
            len(extractor.sites), 0, msg="Could not parse VCF sites.")
        self.assertEqual(
            len(extractor.sites), 15, msg="Did not parse right number of sites.")

    def test_extract_sample(self):

        extractor = Extract(self.args)
        samples = get_samples(self.args)
        samples = extractor.extract(samples)

        self.assertEqual(samples[0].name, 'test_sample1', msg='Sample was not loaded correctly.')
        self.assertIsNotNone(samples[0].pileup, msg='Sample pileup was not loaded correctly.')
        self.assertEqual(samples[0].pileup.shape[0], 15, msg='Did not find pileup for 4 variants. Found: {}.'.format(samples[0].pileup))

    def test_sample_contamination(self):

        extractor = Extract(self.args)
        samples = get_samples(self.args)
        samples = extractor.extract(samples)

        samples = run_minor_contamination(self.args, samples)
        samples = run_major_contamination(self.args, samples)

        self.assertAlmostEqual(
            samples[0].metrics['minor_contamination'], 0.0045, places=4, msg='Minor contamination is wrong.')

        self.assertAlmostEqual(
            samples[0].metrics['major_contamination'], 0.2, places=1, msg='Major contamination is wrong.')
