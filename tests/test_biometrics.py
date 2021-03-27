#!/usr/bin/env python

"""Tests for `biometrics` package."""


import os
import argparse
from unittest import TestCase
from unittest import mock

from biometrics.biometrics import get_samples, run_minor_contamination, run_major_contamination
from biometrics.cli import get_args
from biometrics.extract import Extract
from biometrics.genotype import Genotyper
from biometrics.sex_mismatch import SexMismatch
from biometrics.minor_contamination import MinorContamination
from biometrics.major_contamination import MajorContamination


CUR_DIR = os.path.dirname(os.path.abspath(__file__))


class TestExtract(TestCase):
    """Tests for the extract tool in the `biometrics` package."""

    @mock.patch(
        'argparse.ArgumentParser.parse_args',
        return_value=argparse.Namespace(
            subparser_name='extract',
            input=None,
            sample_bam=[
                os.path.join(CUR_DIR, 'test_data/test_sample1_golden.bam'),
                os.path.join(CUR_DIR, 'test_data/test_sample2_golden.bam')],
            sample_name=['test_sample1', 'test_sample2'],
            sample_type=['tumor', 'tumor'],
            sample_group=['patient1', 'patient1'],
            sample_sex=['M', 'M'],
            database=os.path.join(CUR_DIR, 'test_data/'),
            vcf=os.path.join(CUR_DIR, 'test_data/test.vcf'),
            fafile=os.path.join(CUR_DIR, 'test_data/ref.fasta'),
            bed=os.path.join(CUR_DIR, 'test_data/test.bed'),
            min_mapping_quality=1,
            min_base_quality=1,
            min_coverage=10,
            minor_threshold=0.002,
            major_threshold=0.6,
            discordance_threshold=0.05,
            coverage_threshold=50,
            min_homozygous_thresh=0.1,
            zmin=None,
            zmax=None,
            outdir='.',
            json=None,
            plot=True,
            default_genotype=None,
            overwrite=True,
            no_db_compare=False,
            prefix='test',
            version=False,
            threads=1))
    def setUp(self, mock_args):
        """Set up test fixtures, if any."""

        self.args = get_args()

    def test_load_vcf(self):
        """Test loading the VCF file."""

        extractor = Extract(self.args)

        self.assertGreater(
            len(extractor.sites), 0, msg="Could not parse VCF sites.")
        self.assertEqual(
            len(extractor.sites), 15,
            msg="Did not parse right number of sites.")

    def test_load_bed(self):
        """Test loading the BED file."""

        extractor = Extract(self.args)

        self.assertEqual(
            len(extractor.regions), 1, msg="Expected 1 region in BED file.")

    def test_extract_sample(self):

        extractor = Extract(self.args)
        samples = get_samples(self.args, extraction_mode=True)
        samples = extractor.extract(samples)

        self.assertEqual(len(samples), 2, msg='Did not load 2 samples.')
        self.assertEqual(samples['test_sample1'].sample_name, 'test_sample1', msg='Sample was not loaded correctly.')
        self.assertIsNotNone(samples['test_sample1'].pileup, msg='Sample pileup was not loaded correctly.')
        self.assertEqual(samples['test_sample1'].pileup.shape[0], 15, msg='Did not find pileup for 4 variants. Found: {}.'.format(samples['test_sample1'].pileup))


class TestLoadData(TestCase):
    """Tests load data by sample name in `biometrics` package."""

    @mock.patch(
        'argparse.ArgumentParser.parse_args',
        return_value=argparse.Namespace(
            subparser_name='extract',
            input=['test_sample1', 'test_sample2'],
            sample_bam=None,
            sample_name=None,
            sample_type=None,
            sample_group=None,
            sample_sex=None,
            database=os.path.join(CUR_DIR, 'test_data/'),
            vcf=None,
            fafile=None,
            bed=None,
            min_mapping_quality=None,
            min_base_quality=None,
            min_coverage=None,
            minor_threshold=0.002,
            major_threshold=0.6,
            discordance_threshold=0.05,
            coverage_threshold=50,
            min_homozygous_thresh=0.1,
            zmin=None,
            zmax=None,
            outdir='.',
            json=None,
            plot=True,
            default_genotype=None,
            overwrite=True,
            no_db_compare=False,
            prefix='test',
            version=False,
            threads=1))
    def setUp(self, mock_args):
        """Set up test fixtures, if any."""

        self.args = get_args()

    def test_load_by_sample_name(self):
        samples = get_samples(self.args, extraction_mode=False)

        self.assertEqual(len(samples), 2, msg='Did not load two samples')


class TestLoadDataPickle(TestCase):
    """Tests load data by sample name in `biometrics` package."""

    @mock.patch(
        'argparse.ArgumentParser.parse_args',
        return_value=argparse.Namespace(
            subparser_name='extract',
            input=[
                os.path.join(CUR_DIR, 'test_data/test_sample1.pk'),
                os.path.join(CUR_DIR, 'test_data/test_sample2.pk')],
            sample_bam=None,
            sample_name=None,
            sample_type=None,
            sample_group=None,
            sample_sex=None,
            database=os.path.join(CUR_DIR, 'test_data/'),
            vcf=None,
            fafile=None,
            bed=None,
            min_mapping_quality=None,
            min_base_quality=None,
            min_coverage=None,
            minor_threshold=0.002,
            major_threshold=0.6,
            discordance_threshold=0.05,
            coverage_threshold=50,
            min_homozygous_thresh=0.1,
            zmin=None,
            zmax=None,
            outdir='.',
            json=None,
            plot=True,
            default_genotype=None,
            overwrite=True,
            no_db_compare=False,
            prefix='test',
            version=False,
            threads=1))
    def setUp(self, mock_args):
        """Set up test fixtures, if any."""

        self.args = get_args()

    def test_load_by_pickle_file(self):
        samples = get_samples(self.args, extraction_mode=False)

        self.assertEqual(len(samples), 2, msg='Did not load two samples')


class TestDownstreamTools(TestCase):
    """Tests for downstream tools in the `biometrics` package."""

    @mock.patch(
        'argparse.ArgumentParser.parse_args',
        return_value=argparse.Namespace(
            subparser_name='extract',
            input=['test_sample1', 'test_sample2'],
            sample_bam=None,
            sample_name=None,
            sample_type=None,
            sample_group=None,
            sample_sex=None,
            database=os.path.join(CUR_DIR, 'test_data/'),
            vcf=None,
            fafile=None,
            bed=None,
            min_mapping_quality=None,
            min_base_quality=None,
            min_coverage=None,
            minor_threshold=0.002,
            major_threshold=0.6,
            discordance_threshold=0.05,
            coverage_threshold=50,
            min_homozygous_thresh=0.1,
            zmin=None,
            zmax=None,
            outdir='.',
            json=None,
            plot=True,
            default_genotype=None,
            overwrite=True,
            no_db_compare=False,
            prefix='test',
            version=False,
            threads=1))
    def setUp(self, mock_args):
        """Set up test fixtures, if any."""

        self.args = get_args()

    def test_minor_contamination(self):
        samples = get_samples(self.args, extraction_mode=False)
        samples = run_minor_contamination(self.args, samples)

        self.assertAlmostEqual(
            samples['test_sample1'].metrics['minor_contamination']['val'], 0.0043,
            places=4, msg='Minor contamination is wrong.')

        self.assertEqual(
            samples['test_sample1'].metrics['minor_contamination']['n_contributing_sites'], 1,
            msg='Count of contributing sites for minor contamination is wrong.')

    def test_plot_minor_contamination(self):
        samples = get_samples(self.args, extraction_mode=False)
        minor_contamination = MinorContamination(threshold=self.args.minor_threshold)
        samples = minor_contamination.estimate(samples)
        minor_contamination.plot(samples, self.args.outdir)

    def test_major_contamination(self):
        samples = get_samples(self.args, extraction_mode=False)
        samples = run_major_contamination(self.args, samples)

        self.assertAlmostEqual(
            samples['test_sample1'].metrics['major_contamination']['val'], 0.2,
            places=1, msg='Major contamination is wrong.')

    def test_plot_major_contamination(self):
        samples = get_samples(self.args, extraction_mode=False)
        major_contamination = MajorContamination(threshold=self.args.minor_threshold)
        samples = major_contamination.estimate(samples)
        major_contamination.plot(samples, self.args.outdir)

    def test_genotyper(self):
        samples = get_samples(self.args, extraction_mode=False)

        genotyper = Genotyper(
            no_db_compare=self.args.no_db_compare,
            discordance_threshold=self.args.discordance_threshold,
            threads=self.args.threads,
            zmin=self.args.zmin,
            zmax=self.args.zmax)
        data = genotyper.compare_samples(samples)

        self.assertEqual(len(data), 4, msg='There were not four comparisons done.')
        self.assertEqual(set(data['Status']), set(['Expected Match']), msg='All sample comparisons were expected to match.')

    def test_genotyper_plot(self):
        samples = get_samples(self.args, extraction_mode=False)

        genotyper = Genotyper(
            no_db_compare=self.args.no_db_compare,
            discordance_threshold=self.args.discordance_threshold,
            threads=self.args.threads,
            zmin=self.args.zmin,
            zmax=self.args.zmax)
        data = genotyper.compare_samples(samples)
        genotyper.plot(data, self.args.outdir)

    def test_sexmismatch(self):
        samples = get_samples(self.args, extraction_mode=False)

        sex_mismatch = SexMismatch(self.args.coverage_threshold)
        results = sex_mismatch.detect_mismatch(samples)

        self.assertEqual(set(results['expected_sex']), set(['M']), msg='Expected all samples to have an expected sex of M.')
        self.assertEqual(set(results['predicted_sex']), set(['M']), msg='Expected all samples to not have a sex mismatch.')
