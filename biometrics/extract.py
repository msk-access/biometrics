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
        self.db = args.db
        self.min_mapping_quality = args.min_mapping_quality
        self.min_base_quality = args.min_base_quality
        self.vcf = args.vcf

        self.parse_vcf()

    def parse_vcf(self):
        self.sites = []

        for record in vcf.Reader(open(self.vcf, 'r')):
            self.sites.append({
                'chrom': record.CHROM,
                'start': record.POS,
                'end': record.POS + 1
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

        bam = AlignmentFile(sample.alignment_file)
        pileup = pd.DataFrame()

        for site in self.sites:

            pileup_site = pysamstats.load_pileup(
                'variation', bam, chrom=site['chrom'], start=site['start'],
                end=site['end'], truncate=True,
                max_depth=30000, min_baseq=self.min_base_quality,
                min_mapq=self.min_mapping_quality)

            pileup = pd.concat([pileup, pd.DataFrame(pileup_site)])

        sample.pileup = pileup

        self._save_to_file(sample)

    def extract(self, samples):

        for i, sample in enumerate(samples):

            # if extraction file exists then skip

            if os.path.exists(sample.extraction_file):

                samples[i] = self._load_from_file(sample)

                continue

            samples[i] = self._extract_sample(sample)

        return samples
