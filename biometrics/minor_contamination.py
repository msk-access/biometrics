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
            columns=['sample', 'sample_group', 'sample_sex', 'sample_type',
                     'total_homozygous_sites', 'minor_contamination'])

        for sample_name, sample in samples.items():

            row = {
                'sample': sample_name,
                'sample_group': sample.group,
                'sample_sex': sample.sex,
                'sample_type': sample.sample_type,
                'total_homozygous_sites': sample.metrics['total_homozygous_sites'],
                'minor_contamination': sample.metrics['minor_contamination']
            }

            data = data.append(row, ignore_index=True)

        return data

    def plot(self, data, outdir):

        ymax = max(self.threshold * 1.05, max(data['minor_contamination']))

        fig = go.Figure()
        fig = go.Figure(go.Heatmap(
            x=data['n1'],
            y=data['n2'],
            z=data['val'],
            customdata=data.to_numpy(),
            hovertemplate='<b>Sample group:</b> %{customdata[1]}' +
                          '<br><b>Sample name:</b> %{customdata[0]}' +
                          '<br><b>Sample sex:</b> %{customdata[2]}' +
                          '<br><b>Sample type:</b> %{customdata[3]}' +
                          '<br><b>Total homozygous sites:</b> %{customdata[5]}' +
                          '<br><b>Minor contamination:</b> %{y:E}'))
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
            line=dict(color='Red',),
            xref='x',
            yref='y')
        fig.write_html(os.path.join(outdir, 'minor_contamination.html'))

    def estimate(self, samples):

        for sample_name, sample in samples.items():

            sites = samples[sample_name].pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]

            hom_sites = sites_notna[sites_notna['genotype_class'] == 'Hom']

            sample.metrics = {
                'total_homozygous_sites': hom_sites.shape[0]
            }

            if hom_sites.shape[0] == 0:
                sample.metrics['minor_contamination'] = np.nan
            else:
                sample.metrics['minor_contamination'] = \
                    round(hom_sites['minor_allele_freq'].mean(), 4)

        return samples
