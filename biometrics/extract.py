import os
from multiprocessing import Pool

import pandas as pd
import numpy as np
import vcf
from pysam import AlignmentFile


HETEROZYGOUS_THRESHOLD = 0.1


class Extract:
    """
    Class for extracting genotype information from alignment file using
    the user supplied VCF file.
    """

    def __init__(self, args):
        self.db = args.database
        self.num_threads_samples = args.num_threads_samples
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

        coverage = sum(allele_counts)

        if coverage <= self.min_coverage or coverage == 0:
            return np.nan
        else:
            return min(allele_counts) / coverage

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

    def _get_genotype_info(self, pileup_site, ref_allele, alt_allele):

        # allele_counts = [pileup_site[ref_allele][0], pileup_site[alt_allele][0]]
        allele_counts = [pileup_site[ref_allele], pileup_site[alt_allele]]

        pileup_site['minor_allele_freq'] = self._get_minor_allele_freq(
            allele_counts)

        pileup_site['genotype_class'] = self._get_genotype_class(
            pileup_site['minor_allele_freq'])

        pileup_site['genotype'] = self._get_genotype(
            pileup_site['genotype_class'], allele_counts,
            [ref_allele, alt_allele])

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

    def _pileup_pysam(self, bam, site):

        ignore_overlap = True

        read_quals = {}
        read_bases = {}
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

                if (mapq < self.min_mapping_quality) or \
                        (pileupread.is_del or pileupread.is_refskip):
                    # skip the read if its mapping quality is too low
                    # or if the site is part of an indel
                    continue

                if ignore_overlap and \
                        read_name in read_quals and \
                        read_quals[read_name] > base_qual:
                    # if the read's mate pair has already been counted,
                    # then skip the current read if its base quality is less
                    # than the base quality of its mate pair
                    continue

                read_quals[read_name] = base_qual

                if read_name in read_bases and not ignore_overlap:
                    read_bases[read_name].append(base)
                else:
                    read_bases[read_name] = [base]

        total = 0
        matches = 0
        mismatches = 0
        for read, bases in read_bases.items():
            for base in bases:
                allele_counts[base] += 1
                total += 1

                if base == site['ref_allele']:
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

        if not self.sites:
            return sample

        # get the pileup

        bam = AlignmentFile(sample.alignment_file)
        pileup = pd.DataFrame()

        for site in self.sites:

            pileup_site = self._pileup_pysam(bam, site)

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

    def _extract(self, sample):
        sample = self._extract_sites(sample)
        sample = self._extract_regions(sample)
        sample.save_to_file()

        return sample

    def extract(self, samples):

        # determine with samples need to be extracted

        samples_to_extract = []

        for i, sample_name in enumerate(samples):

            sample = samples[sample_name]

            # if extraction file exists then load it

            if os.path.exists(sample.extraction_file) and not self.overwrite:
                samples[sample_name].load_from_file()
                continue

            samples_to_extract.append(sample)

        # if any samples need to be extracted, then do so

        if len(samples_to_extract) > 0:

            thread_pool = Pool(self.num_threads_samples)

            samples_processed = thread_pool.map(
                self._extract, samples_to_extract)

            for sample in samples_processed:
                samples[sample.name] = sample

        return samples
