import os
import pickle

import pandas as pd
import numpy as np
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

        pileup_data = sample.pileup.to_dict()
        sample_data = {
            'alignment_file': sample.alignment_file,
            'name': sample.name,
            'sex': sample.sex,
            'patient': sample.patient,
            'sample_type': sample.sample_type,
            'pileup_data': pileup_data
        }

        pickle.dump(sample_data, open(sample.extraction_file, "wb"))

    def _load_from_file(self, sample):
        sample_data = pickle.load(open(sample.extraction_file, "rb"))
        sample.pileup = sample_data['pileup']
        sample.alignment_file = sample_data['alignment_file']
        sample.name = sample_data['name']
        sample.sex = sample_data['sex']
        sample.patient = sample_data['patient']
        sample.sample_type = sample_data['sample_type']

        return sample

    def _get_pileup_allele_count(self, pileup_site, allele):
        return pileup_site[allele][0]

    def _get_minor_allele_freq(self, pileup_site, site):
        allele_counts = [
            self._get_pileup_allele_count(pileup_site, site['ref_allele']),
            self._get_pileup_allele_count(pileup_site, site['alt_allele'])
        ]

        if sum(allele_counts) == 0:
            return np.nan
        else:
            return min(allele_counts) / sum(allele_counts)

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

            pileup_site = pd.DataFrame(pileup_site)
            pileup_site['minor_allele_freq'] = self._get_minor_allele_freq(
                pileup_site, site)

            pileup = pd.concat([pileup, pileup_site])

        sample.pileup = pileup
        sample.pileup = sample.pileup[[
            'chrom', 'pos', 'ref', 'reads_all', 'matches', 'mismatches', 'A',
            'C', 'T', 'G', 'N', 'minor_allele_freq']]

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
