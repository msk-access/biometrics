import os

import pandas as pd
import numpy as np
import vcf
import pysamstats
from pysam import AlignmentFile


HETEROZYGOUS_THRESHOLD = 0.1


class Extract:
    """
    Class for extracting genotype information from alignment file using
    the user supplied VCF file.
    """

    def __init__(self, args):
        self.db = args.database
        self.min_mapping_quality = args.min_mapping_quality
        self.min_base_quality = args.min_base_quality
        self.vcf = args.vcf
        self.bed = args.bed
        self.fafile = args.fafile
        self.overwrite = args.overwrite
        self.min_coverage = args.min_coverage
        self.sites = []
        self.regions = None

        self.parse_vcf()
        self.parse_bed_file()

    def parse_vcf(self):

        if self.vcf is None:
            return

        for record in vcf.Reader(open(self.vcf, 'r')):
            self.sites.append({
                'chrom': record.CHROM,
                'start': record.POS-1,
                'end': record.POS,
                'ref_allele': str(record.REF),
                'alt_allele': str(record.ALT[0])
            })

    def parse_bed_file(self):

        if self.bed is None:
            return

        self.regions = pd.read_csv(self.bed, sep='\t')
        self.regions.columns = range(self.regions.shape[1])

    def _get_minor_allele_freq(self, allele_counts):
        if sum(allele_counts) <= self.min_coverage:
            return np.nan
        else:
            return min(allele_counts) / sum(allele_counts)

    def _get_genotype_class(self, minor_allele_freq):
        if pd.isna(minor_allele_freq):
            return np.nan
        else:
            if minor_allele_freq <= HETEROZYGOUS_THRESHOLD:
                return 'Hom'
            else:
                return 'Het'

    def _get_genotype(self, genotype, allele_counts, alleles):

        if pd.isna(genotype):
            return np.nan
        elif genotype == 'Het':
            return ''.join(alleles)
        else:
            if allele_counts[0] > allele_counts[1]:
                return alleles[0]
            else:
                return alleles[1]

    def _get_genotype_info(self, pileup_site):

        allele_counts = [
            pileup_site[pileup_site['ref_allele'][0]][0],
            pileup_site[pileup_site['alt_allele'][0]][0]]

        pileup_site['minor_allele_freq'] = self._get_minor_allele_freq(
            allele_counts)

        pileup_site['genotype_class'] = self._get_genotype_class(
            pileup_site['minor_allele_freq'][0])

        pileup_site['genotype'] = self._get_genotype(
            pileup_site['genotype_class'][0], allele_counts,
            [pileup_site['ref_allele'][0], pileup_site['alt_allele'][0]])

        return pileup_site

    def _extract_regions(self, sample):

        if self.regions is None:
            return sample

        # get the pileup

        bam = AlignmentFile(sample.alignment_file)
        region_counts = []

        for i in self.regions.index:

            chrom = self.regions.at[i, 0]
            start = int(self.regions.at[i, 1])
            end = int(self.regions.at[i, 2])

            count = bam.count(chrom, start, end)

            region_counts.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'count': count})

        region_counts = pd.DataFrame(region_counts)

        sample.region_counts = region_counts

        return sample

    def _extract_sites(self, sample):

        if not self.sites:
            return sample

        # get the pileup

        bam = AlignmentFile(sample.alignment_file)
        pileup = pd.DataFrame()

        for site in self.sites:

            pileup_site = pysamstats.load_pileup(
                'variation', bam, chrom=site['chrom'], start=site['start'],
                end=site['end'], truncate=True, fafile=self.fafile,
                max_depth=30000, min_baseq=self.min_base_quality,
                min_mapq=self.min_mapping_quality)

            pileup_site = pd.DataFrame(pileup_site)

            pileup_site['ref_allele'] = site['ref_allele']
            pileup_site['alt_allele'] = site['alt_allele']

            pileup_site = self._get_genotype_info(pileup_site)

            pileup = pd.concat([pileup, pileup_site])

        sample.pileup = pileup
        sample.pileup = sample.pileup[[
            'chrom', 'pos', 'ref', 'reads_all', 'matches', 'mismatches', 'A',
            'C', 'T', 'G', 'N', 'minor_allele_freq', 'genotype_class', 'genotype']]

        # because pysamstats works in 0-based coordinates. so ned to convert to
        # 1-based
        sample.pileup['pos'] = sample.pileup['pos'] + 1

        return sample

    def extract(self, samples):

        for i, sample_name in enumerate(samples):

            sample = samples[sample_name]

            # if extraction file exists then load it

            if os.path.exists(sample.extraction_file) and not self.overwrite:

                samples[sample_name].load_from_file()

                continue

            sample = self._extract_sites(sample)
            sample = self._extract_regions(sample)
            sample.save_to_file()

            samples[sample_name] = sample

        return samples
