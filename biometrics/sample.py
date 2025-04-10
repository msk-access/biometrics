import pickle
import os

import pandas as pd
import pdb


class Sample:
    """
    Class to hold information related to a single sample.
    """

    def __init__(self, sample_name=None, sample_bam=None, sample_group=None,
                 sample_sex=None, sample_type=None, db=None, query_group=False):
        self.sample_bam = sample_bam
        self.sample_name = sample_name
        self.sample_sex = sample_sex
        self.sample_type = sample_type

        if sample_group is None:
            self.sample_group = sample_name
        else:
            self.sample_group = sample_group

        self.pileup = None
        self.region_counts = None
        self.extraction_file = None
        self.query_group = query_group
        self.metrics = {}

        if self.sample_name is not None:
            if db is not None:
                self.extraction_file = os.path.join(db, self.sample_name + '.pickle')
            else:
                self.extraction_file = self.sample_name + '.pickle'

    def save_to_file(self):

        pileup_data = self.pileup.to_dict("records")

        if self.region_counts is not None:
            region_counts = self.region_counts.to_dict('records')
        else:
            region_counts = None

        sample_data = {
            'sample_bam': self.sample_bam,
            'sample_name': self.sample_name,
            'sample_sex': self.sample_sex,
            'sample_group': self.sample_group,
            'sample_type': self.sample_type,
            'pileup_data': pileup_data,
            'region_counts': region_counts
        }
        pickle.dump(sample_data, open(self.extraction_file, "wb"))

        #################################################
        #converting pileup data to FP summary like csv file
        #get the alt count based on the alt on the genotype
        def match_column_letters(row):
            matched_values = []
            for letter in row['genotype']:  # Iterate through each character in 'genotype'
                for col in row.index:
                    if letter in col:  # Check if character is in column name
                        matched_values.append(f"{letter}:{row[col]}")
            return ",".join(matched_values) if matched_values else None

        summary_file = "ALL_FPsummary.txt"
        #make sure to remove older txt file. But this is neccessary to add the samples to existing df
        if os.path.exists(summary_file):
            fp_summary = pd.read_csv(summary_file)
        else:
            fp_summary = pd.DataFrame()

        #convert existing pileup data as dataframe
        pileup_df =pd.DataFrame.from_dict(self.pileup)
        pileup_df = pileup_df[pileup_df['genotype_class'].notna()]
        sample_name= self.sample_name

        new_sample_data = pd.DataFrame()
        new_sample_data['Locus'] = pileup_df['chrom'].astype(str) + ":" + pileup_df['pos'].astype(str)
        new_sample_data[sample_name + '_Counts'] = pileup_df.apply(match_column_letters, axis=1)
        new_sample_data[sample_name + '_Genotypes'] = pileup_df['genotype']
        new_sample_data[sample_name + '_MinorAlleleFreq'] = pileup_df['minor_allele_freq']

        #merge the dataframe with the empty frame to create a single df with all the samples
        if not fp_summary.empty:
            fp_summary = pd.merge(fp_summary, new_sample_data, on='Locus', how='outer')
        else:
            fp_summary = new_sample_data
        fp_summary.to_csv(summary_file, index=False)

    def load_from_file(self, extraction_file=None):

        if extraction_file is not None:
            self.extraction_file = extraction_file

        assert self.extraction_file is not None, 'Extraction file path is None.'
        assert os.path.exists(self.extraction_file), 'Extraction file does not exist: {}'.format(
            self.extraction_file)

        sample_data = pickle.load(open(self.extraction_file, "rb"))

        region_counts = None
        if sample_data.get('region_counts') is not None:
            region_counts = pd.DataFrame(
                sample_data['region_counts'], dtype=object)

        self.pileup = pd.DataFrame(sample_data['pileup_data'])
        self.sample_bam = sample_data['sample_bam']
        self.sample_name = sample_data['sample_name']
        self.sample_sex = sample_data['sample_sex']
        self.sample_group = sample_data['sample_group'] if sample_data['sample_group'] is not None else sample_data['sample_name']
        self.sample_type = sample_data['sample_type']
        self.region_counts = region_counts
