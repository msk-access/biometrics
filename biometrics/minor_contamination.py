import os

import pandas as pd
import numpy as np
import plotly.graph_objects as go


class MinorContamination():
    """
    Minor contamination.
    """

    def __init__(self, args):
        pass

    def to_dataframe(self, samples):

        data = pd.DataFrame(
            columns=['sample', 'sample_group', 'sample_sex', 'sample_type',
                     'minor_contamination', 'total_homozygous_sites'])

        for sample_name, sample in samples.items():

            row = {
                'sample': sample_name,
                'sample_group': sample.group,
                'sample_sex': sample.sex,
                'sample_type': sample.sample_type,
                'minor_contamination': sample.metrics['minor_contamination'],
                'total_homozygous_sites': sample.metrics['total_homozygous_sites']
            }

            data = data.append(row, ignore_index=True)

        return data

    def plot(self, data, outdir):

        ymax = max(0.0025, max(data['minor_contamination']))

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                x=data['sample'],
                y=data['minor_contamination'],
                customdata=data.to_numpy(),
                hovertemplate='''<b>Sample group:</b> %{customdata[1]}<br><b>Sample name:</b> %{customdata[0]}<br><b>Sample sex:</b> %{customdata[2]}<br><b>Sample type:</b> %{customdata[3]}<br><b>Total homozygous sites:</b> %{customdata[5]}''',
            ))
        fig.update_layout(
            yaxis_title="Minor contamination",
            title_text="Minor contamination across samples",
            yaxis=dict(range=[0, ymax]))
        fig.add_shape(
            type='line',
            x0=-1,
            y0=0.002,
            x1=8,
            y1=0.002,
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
