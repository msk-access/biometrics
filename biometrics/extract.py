import os
from multiprocessing import Pool

import pandas as pd
import numpy as np
import vcf
from pysam import AlignmentFile


class Extract:
    """
    Class for extracting genotype information from alignment file using
    the user supplied VCF file.
    """

    def __init__(self, args):
        self.db = args.database
        self.threads = args.threads
        self.min_mapping_quality = args.min_mapping_quality
        self.min_base_quality = args.min_base_quality
        self.default_genotype = args.default_genotype
        self.vcf = args.vcf
        self.bed = args.bed
        self.fafile = args.fafile
        self.overwrite = args.overwrite
        self.min_coverage = args.min_coverage
        self.min_homozygous_thresh = args.min_homozygous_thresh
        self.sites = []
        self.regions = None

        self._parse_vcf()
        self._parse_bed_file()

    def _parse_vcf(self):

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

    def _parse_bed_file(self):

        if self.bed is None:
            return

        self.regions = pd.read_csv(self.bed, sep='\t', header=None)
        self.regions.columns = range(self.regions.shape[1])

    def _extract_regions(self, sample):
        """
        Code to extract the coverage information for the regions listed
        in the BED file.
        """

        if self.regions is None:
            return sample

        # get the pileup

        bam = AlignmentFile(sample.sample_bam)
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

    def _get_minor_allele_freq(self, allele_counts):

        coverage = sum(allele_counts)

        if coverage < self.min_coverage or coverage == 0:
            return np.nan
        else:
            return min(allele_counts) / coverage

    def _get_genotype_class(self, minor_allele_freq):
        """
        Determine if Het, Hom, or unknown/NA.
        """

        if pd.isna(minor_allele_freq):

            if self.default_genotype is not None:
                return self.default_genotype

            return np.nan
        else:
            if minor_allele_freq <= self.min_homozygous_thresh:
                return 'Hom'
            else:
                return 'Het'

    def _get_genotype(self, genotype, allele_counts, alleles):
        """
        Get the genotype in terms of the allele(s) (e.g. A, T, AT, GC, etc.)
        """

        if pd.isna(genotype):
            return np.nan
        elif genotype == 'Het':
            return ''.join(alleles)
        else:
            if allele_counts[0] > allele_counts[1]:
                return alleles[0]
            else:
                return alleles[1]

    def _get_genotype_info(self, pileup_site, ref_allele, alt_allele):
        """
        Plot minor contamination data.
        """

        allele_counts = [pileup_site[ref_allele], pileup_site[alt_allele]]

        pileup_site['minor_allele_freq'] = self._get_minor_allele_freq(
            allele_counts)

        pileup_site['genotype_class'] = self._get_genotype_class(
            pileup_site['minor_allele_freq'])

        pileup_site['genotype'] = self._get_genotype(
            pileup_site['genotype_class'], allele_counts,
            [ref_allele, alt_allele])

        return pileup_site

    def _add_base(self, site, old_base, old_base_qual, new_base,
                  new_base_qual):
        """
        This function is for dealing with the various scenarios that can
        arise when a read pair overlaps and how to handle when the
        bases mismatch. The 'old_base' refers to the first base observed when
        computing pileup information (usually the forward read). Then the
        'new_base' is from the second read in the overlaping pair.
        """

        if old_base is None:
            return [new_base, new_base_qual]

        if old_base == new_base:
            return [old_base, old_base_qual]

        if old_base != 'N' and new_base != 'N':
            if new_base == site['ref_allele']:
                return [new_base, new_base_qual]
            else:
                return [old_base, old_base_qual]

        if old_base == site['ref_allele']:
            return [old_base, old_base_qual]
        elif new_base == site['ref_allele']:
            return [new_base, new_base_qual]
        elif old_base == site['alt_allele'] and old_base_qual >= '!':
            return [old_base, old_base_qual]
        elif new_base == site['alt_allele'] and new_base_qual >= '!':
            return [new_base, new_base_qual]
        else:
            return ['N', '&']

    def _pileup(self, bam, site):
        """
        Get the per-site pileup information.
        """

        read_data = {}
        allele_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}

        for pileupcolumn in bam.pileup(
                contig=site['chrom'], start=site['start'], end=site['end'],
                truncate=True, max_depth=30000, stepper='nofilter',
                min_base_quality=self.min_base_quality):

            for pileupread in pileupcolumn.pileups:

                if pileupread.query_position is None:
                    continue

                mapq = pileupread.alignment.mapping_quality
                read_name = pileupread.alignment.qname
                base_qual = pileupread.alignment.qual[pileupread.query_position]
                base = pileupread.alignment.query_sequence[pileupread.query_position]

                if (mapq < self.min_mapping_quality) or pileupread.is_refskip or pileupread.is_del:
                    # skip the read if its mapping quality is too low
                    # or if the site is part of an indel
                    continue

                if read_name in read_data and read_data[read_name][0] == 'N':
                    continue
                elif read_name in read_data:
                    vals = self._add_base(
                        site, read_data[read_name][0],
                        read_data[read_name][1], base, base_qual)
                    read_data[read_name] = vals[0:2]
                else:
                    read_data[read_name] = [base, base_qual]

        total = 0
        matches = 0
        mismatches = 0
        for read, base_data in read_data.items():

            allele_counts[base_data[0]] += 1
            total += 1

            if base_data[0] == site['ref_allele']:
                matches += 1
            else:
                mismatches += 1

        return {
            'chrom': site['chrom'],
            'pos': site['end'],
            'ref': site['ref_allele'],
            'alt': site['alt_allele'],
            'reads_all': total,
            'matches': matches,
            'mismatches': mismatches,
            'A': allele_counts['A'],
            'C': allele_counts['C'],
            'T': allele_counts['T'],
            'G': allele_counts['G'],
            'N': allele_counts['N']
        }

    def _extract_sites(self, sample):
        """
        Loop through all positions and get pileup information.
        """

        if not self.sites:
            return sample

        # get the pileup

        bam = AlignmentFile(sample.sample_bam)
        pileup = pd.DataFrame()

        for site in self.sites:

            pileup_site = self._pileup(bam, site)

            pileup_site = self._get_genotype_info(
                pileup_site, site['ref_allele'], site['alt_allele'])

            pileup = pileup.append(pileup_site, ignore_index=True)

        pileup = pileup[[
            'chrom', 'pos', 'ref', 'alt', 'reads_all', 'matches', 'mismatches',
            'A', 'C', 'T', 'G', 'N', 'minor_allele_freq', 'genotype_class',
            'genotype']]

        for col in ['pos', 'A', 'C', 'T', 'G', 'N', 'matches', 'mismatches', 'reads_all']:
            pileup[col] = pileup[col].astype(int)

        sample.pileup = pileup

        return sample

    def _extraction_job(self, sample):
        """
        Function to do the extraction steps for a single sample.
        Supposed to be called by multiprocessing functions to parallelize it.
        """

        sample = self._extract_sites(sample)
        sample = self._extract_regions(sample)
        sample.save_to_file()

        return sample

    def extract(self, samples):
        """
        Function to call to extract the pileup and region information
        for the given samples.
        """

        if type(samples) != dict:
            samples = {samples.sample_name: samples}

        # determine with samples need to be extracted, and put them in a list

        samples_to_extract = []

        for sample_name, sample in samples.items():

            # if extraction file exists then load it

            if os.path.exists(sample.extraction_file) and not self.overwrite:
                sample.load_from_file()
                continue

            samples_to_extract.append(sample)

        # if any samples need to be extracted, then do so
        # (using multiprocessing)

        if len(samples_to_extract) > 0:

            thread_pool = Pool(self.threads)

            samples_processed = thread_pool.map(
                self._extraction_job, samples_to_extract)

            for sample in samples_processed:
                samples[sample.sample_name] = sample

        return samples
