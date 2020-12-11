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
            columns=['sample', 'sample_group', 'sample_sex', 'sample_type',
                     'major_contamination', 'total_sites',
                     'total_heterozygous_sites'])

        for sample_name, sample in samples.items():

            row = {
                'sample': sample_name,
                'sample_group': sample.group,
                'sample_sex': sample.sex,
                'sample_type': sample.sample_type,
                'major_contamination': sample.metrics['major_contamination'],
                'total_sites': sample.metrics['total_sites'],
                'total_heterozygous_sites': sample.metrics['total_heterozygous_sites']
            }

            data = data.append(row, ignore_index=True)

        return data

    def plot(self, data, outdir):

        ymax = max(self.threshold, max(data['major_contamination'])) * 1.05

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=data['sample'],
                y=data['major_contamination'],
                customdata=data.to_numpy(),
                hovertemplate='<b>Sample group:</b> %{customdata[1]}' +
                              '<br><b>Sample name:</b> %{customdata[0]}' +
                              '<br><b>Sample sex:</b> %{customdata[2]}' +
                              '<br><b>Sample type:</b> %{customdata[3]}' +
                              '<br><b>Total sites:</b> %{customdata[5]}' +
                              '<br><b>Total heterozygous sites:</b> %{customdata[6]}' +
                              '<br><b>Major contamination:</b> %{y:E}',
            ))
        fig.update_layout(
            yaxis_title="Major contamination",
            title_text="Major contamination across samples",
            yaxis=dict(range=[0, ymax]))
        fig.add_shape(
            type='line',
            x0=-1,
            y0=self.threshold,
            x1=data.shape[0],
            y1=self.threshold,
            line=dict(color='Red'),
            xref='x',
            yref='y')

        fig.write_html(os.path.join(outdir, 'major_contamination.html'))

    def estimate(self, samples):

        for sample_name, sample in samples.items():

            sites = sample.pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            het_sites = sites_notna[sites_notna['genotype_class'] == 'Het']

            sample.metrics = {
                'total_sites': sites_notna.shape[0],
                'total_heterozygous_sites': het_sites.shape[0]
            }

            if sites_notna.shape[0] == 0:
                sample.metrics['major_contamination'] = np.nan
            else:
                sample.metrics['major_contamination'] = \
                    round(het_sites.shape[0] / sites_notna.shape[0], 4)

        return samples
