import os
import pickle

import pandas as pd
import vcf
import pysamstats
from pysam import AlignmentFile


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
        self.fafile = args.fafile
        self.overwrite = args.overwrite

        self.parse_vcf()

    def parse_vcf(self):
        self.sites = []

        for record in vcf.Reader(open(self.vcf, 'r')):
            self.sites.append({
                'chrom': record.CHROM,
                'start': record.POS,
                'end': record.POS + 1,
                'ref_allele': record.REF,
                'alt_allele': record.ALT[0]
            })

    def _save_to_file(self, sample):
        pickle.dump(
            sample.pileup.to_dict(),
            open(sample.extraction_file, "wb"))

    def _load_from_file(self, sample):
        pileup = pickle.load(open(sample.extraction_file, "rb"))
        sample.pileup = pileup

        return pileup

    def _extract_sample(self, sample):

        # get the pileup

        bam = AlignmentFile(sample.alignment_file)
        pileup = pd.DataFrame()

        for site in self.sites:
            pileup_site = pysamstats.load_pileup(
                'variation', bam, chrom=site['chrom'], start=site['start'],
                end=site['end'], truncate=True, fafile=self.fafile,
                max_depth=30000, min_baseq=self.min_base_quality,
                min_mapq=self.min_mapping_quality)

            pileup = pd.concat([pileup, pd.DataFrame(pileup_site)])

        sample.pileup = pileup
        sample.pileup = sample.pileup[[
            'chrom', 'pos', 'ref', 'reads_all', 'matches', 'mismatches', 'A',
            'C', 'T', 'G', 'N']]

        # compute some metrics

        self._save_to_file(sample)

    def extract(self, samples):

        for i, sample in enumerate(samples):

            # if extraction file exists then load it

            if os.path.exists(sample.extraction_file) and not self.overwrite:

                samples[i] = self._load_from_file(sample)

                continue

            samples[i] = self._extract_sample(sample)

        return samples
