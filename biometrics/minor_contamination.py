import os

import pandas as pd
import numpy as np
import plotly.graph_objects as go


class MinorContamination():
    """
    Minor contamination.
    """

    def __init__(self, threshold):
        self.threshold = threshold

    def to_dataframe(self, samples):

        data = pd.DataFrame(
            columns=['sample_name', 'sample_group', 'sample_sex', 'sample_type',
                     'total_homozygous_sites', 'minor_contamination'])

        for sample_name, sample in samples.items():

            row = {
                'sample_name': sample.sample_name,
                'sample_group': sample.sample_group,
                'sample_sex': sample.sample_sex,
                'sample_type': sample.sample_type,
                'total_homozygous_sites': sample.metrics['total_homozygous_sites'],
                'minor_contamination': sample.metrics['minor_contamination']
            }

            data = data.append(row, ignore_index=True)

        return data

    def plot(self, data, outdir):
        """
        Plot major contamination data.
        """

        ymax = max(self.threshold, max(data['minor_contamination'])) * 1.05
        data['minor_contamination'] = data['minor_contamination'].map(
            lambda x: round(x, 4))

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=data['sample_name'],
                y=data['minor_contamination'],
                customdata=data.to_numpy(),
                hovertemplate='<b>Sample group:</b> %{customdata[1]}' +
                              '<br><b>Sample name:</b> %{customdata[0]}' +
                              '<br><b>Sample sex:</b> %{customdata[2]}' +
                              '<br><b>Sample type:</b> %{customdata[3]}' +
                              '<br><b>Total homozygous sites:</b> %{customdata[4]}' +
                              '<br><b>Minor contamination:</b> %{y:E}' +
                              '<extra></extra>'))

        fig.update_layout(
            yaxis_title="Minor contamination",
            title_text="Minor contamination across samples",
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

        fig.write_html(os.path.join(outdir, 'minor_contamination.html'))

    def estimate(self, samples):
        """
        Estimate minor contamination.
        """

        for sample_name, sample in samples.items():

            sites = sample.pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            hom_sites = sites_notna[sites_notna['genotype_class'] == 'Hom']

            sample.metrics['total_homozygous_sites'] = len(hom_sites)

            if sample.metrics['total_homozygous_sites'] == 0:
                sample.metrics['minor_contamination'] = np.nan
            else:
                sample.metrics['minor_contamination'] = \
                    hom_sites['minor_allele_freq'].mean()

        return samples
