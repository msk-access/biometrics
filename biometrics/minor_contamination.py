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
                     'total_homozygous_sites', 'n_contributing_sites', 'minor_contamination'])

        for sample_name, sample in samples.items():

            row = {
                'sample_name': sample.sample_name,
                'sample_group': sample.sample_group,
                'sample_sex': sample.sample_sex,
                'sample_type': sample.sample_type,
                'total_homozygous_sites': sample.metrics['minor_contamination']['n_homozygous_sites'],
                'n_contributing_sites': sample.metrics['minor_contamination']['n_contributing_sites'],
                'minor_contamination': sample.metrics['minor_contamination']['val']
            }

            data = data.append(row, ignore_index=True)

        data = data.sort_values('minor_contamination', ascending=False)
        return data

    def plot(self, samples, outdir):
        """
        Plot major contamination data.
        """

        data = self.to_dataframe(samples)
        data['minor_contamination'] = data['minor_contamination'].map(
            lambda x: round(x, 5))

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
                              '<br><b>Total contributing sites:</b> %{customdata[5]}' +
                              '<br><b>Minor contamination:</b> %{y:E}' +
                              '<extra></extra>'))

        fig.update_layout(
            yaxis_title="Minor contamination",
            title_text="Minor contamination across samples")
        fig.add_hline(y=self.threshold, line_color='red')

        fig.write_html(os.path.join(outdir, 'minor_contamination.html'))

        # plot VAF of contributing sites

        plot_data = []

        for i, sample_name in enumerate(data['sample_name']):
            contributing_sites = samples[sample_name].metrics['minor_contamination']['contributing_sites']

            for site_id, site_data in contributing_sites.items():
                site_data['sample_name'] = sample_name
                site_data['MAF'] = site_data['minor_allele_freq']
                site_data['index'] = i

                del site_data['minor_allele_freq']
                plot_data.append(site_data)

        plot_data = pd.DataFrame(plot_data)
        plot_data = plot_data[[
            'sample_name', 'chrom', 'pos', 'ref', 'alt', 'MAF', 'reads_all', 'A', 'C', 'T', 'G',
            'N', 'index']]
        plot_data['MAF'] = plot_data['MAF'].map(lambda x: round(x, 5))


        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=plot_data['index'],
                y=plot_data['MAF'],
                mode='markers',
                showlegend=False,
                customdata=plot_data.to_numpy(),
                hovertemplate='<b>Sample:</b> %{customdata[0]}' +
                              '<br><b>Chrom:</b> %{customdata[1]}' +
                              '<br><b>Pos:</b> %{customdata[2]}' +
                              '<br><b>Ref allele:</b> %{customdata[3]}' +
                              '<br><b>Alt allele:</b> %{customdata[4]}' +
                              '<br><b>MAF:</b> %{customdata[5]}' +
                              '<br><b>Total reads:</b> %{customdata[6]}' +
                              '<br><b>Count A:</b> %{customdata[7]}' +
                              '<br><b>Count C:</b> %{customdata[8]}' +
                              '<br><b>Count T:</b> %{customdata[9]}' +
                              '<br><b>Count G:</b> %{customdata[10]}' +
                              '<br><b>Count N:</b> %{customdata[11]}' +
                              '<extra></extra>'))

        data = data[data['sample_name'].isin(plot_data['sample_name'])]
        data['index'] = range(len(data))

        for i in data.index:
            fig.add_shape(go.layout.Shape(
                type="line",
                x0=data.at[i, 'index'] - 0.5,
                y0=data.at[i, 'minor_contamination'],
                x1=data.at[i, 'index'] + 0.5,
                y1=data.at[i, 'minor_contamination'],
                line=dict(color='black', width=2)
            ))

        fig.add_trace(go.Scatter(
            x=[1],
            y=[3],
            name='Minor contamination',
            line_color='black'
        ))
        fig.add_trace(go.Scatter(
            x=[1],
            y=[3],
            name='Threshold',
            line_color='red'
        ))
        fig.add_trace(go.Scatter(
            x=[1],
            y=[3],
            name='Minor allele sites (MAF > 0)',
            line_color='#636EFA',
            mode='markers'
        ))


        ticks = plot_data[['sample_name', 'index']].drop_duplicates()

        fig.update_layout(
            yaxis_title="Minor allele frequency",
            title_text="Statistics of sites that contribute to minor contamination",
            yaxis = dict(range=(-0.005, plot_data['MAF'].max()*1.05)),
            xaxis = dict(
                tickmode = 'array',
                tickvals = ticks['index'],
                ticktext = ticks['sample_name']
            ))
        fig.add_hline(y=self.threshold, line_color='red')

        fig.write_html(os.path.join(outdir, 'minor_contamination_sites.html'))

        return data

    def estimate(self, samples):
        """
        Estimate minor contamination.
        """

        for sample_name, sample in samples.items():

            sites = sample.pileup
            sites_notna = sites[~pd.isna(sites['genotype_class'])]
            hom_sites = sites_notna[sites_notna['genotype_class'] == 'Hom']

            sample.metrics['minor_contamination'] = {}
            sample.metrics['minor_contamination']['n_homozygous_sites'] = len(hom_sites)

            if len(hom_sites) == 0:
                sample.metrics['minor_contamination']['val'] = np.nan
                sample.metrics['minor_contamination']['n_contributing_sites'] = np.nan
                sample.metrics['minor_contamination']['contributing_sites'] = np.nan
            else:

                contributing_sites = hom_sites[hom_sites['minor_allele_freq']>0]
                contributing_sites.index = contributing_sites['chrom'].astype(str) + ':' + \
                    contributing_sites['pos'].astype(str)

                sample.metrics['minor_contamination']['val'] = \
                    hom_sites['minor_allele_freq'].mean()
                sample.metrics['minor_contamination']['n_contributing_sites'] = \
                    len(contributing_sites)
                sample.metrics['minor_contamination']['contributing_sites'] = \
                    contributing_sites.to_dict(orient='index')


        return samples
