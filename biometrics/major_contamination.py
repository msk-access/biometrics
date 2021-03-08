import os

import pandas as pd
import numpy as np
import plotly.graph_objects as go


class MajorContamination():
    """
    Major contamination.
    """

    def __init__(self, threshold):
        self.threshold = threshold

    def to_dataframe(self, samples):

        data = pd.DataFrame(
            columns=['sample_name', 'sample_group', 'sample_sex', 'sample_type',
                     'total_sites', 'total_heterozygous_sites',
                     'major_contamination'])

        for sample_name, sample in samples.items():

            row = {
                'sample_name': sample.sample_name,
                'sample_group': sample.sample_group,
                'sample_sex': sample.sample_sex,
                'sample_type': sample.sample_type,
                'total_sites': sample.metrics['major_contamination']['total_sites'],
                'total_heterozygous_sites': sample.metrics['major_contamination']['total_heterozygous_sites'],
                'major_contamination': sample.metrics['major_contamination']['val']
            }

            data = data.append(row, ignore_index=True)

        data = data.sort_values('major_contamination', ascending=False)
        return data

    def plot(self, samples, outdir):
        """
        Plot minor contamination data.
        """

        data = self.to_dataframe(samples)
        data['major_contamination'] = data['major_contamination'].map(
            lambda x: round(x, 5))

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=data['sample_name'],
                y=data['major_contamination'],
                customdata=data.to_numpy(),
                hovertemplate='<b>Sample group:</b> %{customdata[1]}' +
                              '<br><b>Sample name:</b> %{customdata[0]}' +
                              '<br><b>Sample sex:</b> %{customdata[2]}' +
                              '<br><b>Sample type:</b> %{customdata[3]}' +
                              '<br><b>Total sites:</b> %{customdata[5]}' +
                              '<br><b>Total heterozygous sites:</b> %{customdata[6]}' +
                              '<br><b>Major contamination:</b> %{y:E}' +
                              '<extra></extra>',
            ))
        fig.update_layout(
            yaxis_title="Major contamination",
            title_text="Major contamination across samples")
        fig.add_hline(y=self.threshold, line_color='red')

        fig.write_html(os.path.join(outdir, 'major_contamination.html'))

    def estimate(self, samples):
        """
        Estimate major contamination.
        """

        for sample_name, sample in samples.items():

            sites = sample.pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            het_sites = sites_notna[sites_notna['genotype_class'] == 'Het']

            sample.metrics['major_contamination'] = {}
            sample.metrics['major_contamination']['total_sites'] = len(sites_notna)
            sample.metrics['major_contamination']['total_heterozygous_sites'] = len(het_sites)

            if sample.metrics['major_contamination']['total_sites'] == 0:
                sample.metrics['major_contamination']['val'] = np.nan
            else:
                sample.metrics['major_contamination']['val'] = \
                    len(het_sites) / len(sites_notna)

        return samples
