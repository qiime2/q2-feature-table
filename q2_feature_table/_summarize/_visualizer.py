# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import shutil

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from q2_types.feature_data import DNAIterator
import q2templates
import skbio
import qiime2

_blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                       "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                       "=nt&CMD=Put&QUERY=%s")

TEMPLATES = pkg_resources.resource_filename('q2_feature_table', '_summarize')


def tabulate_seqs(output_dir: str, data: DNAIterator) -> None:
    sequences = []
    with open(os.path.join(output_dir, 'sequences.fasta'), 'w') as fh:
        for sequence in data:
            skbio.io.write(sequence, format='fasta', into=fh)
            str_seq = str(sequence)
            sequences.append({'id': sequence.metadata['id'],
                              'url': _blast_url_template % str_seq,
                              'seq': str_seq})

    index = os.path.join(TEMPLATES, 'tabulate_seqs_assets', 'index.html')
    q2templates.render(index, output_dir, context={'data': sequences})

    js = os.path.join(
        TEMPLATES, 'tabulate_seqs_assets', 'js', 'tsorter.min.js')
    os.mkdir(os.path.join(output_dir, 'js'))
    shutil.copy(js, os.path.join(output_dir, 'js', 'tsorter.min.js'))


def summarize(output_dir: str, table: pd.DataFrame,
              sample_metadata: qiime2.Metadata=None) -> None:
    number_of_samples, number_of_features = table.shape
    counts = table.sum(axis=1)
    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("load(")
        if sample_metadata:
            sample_metadata.to_dataframe().to_json(fh)
        else:
            fh.write('{}')
        fh.write(', ')
        counts.to_json(fh)
        fh.write(');')

    sample_summary, sample_frequencies = _frequency_summary(
        table, axis=1)
    if number_of_samples > 1:
        # Freedman–Diaconis rule
        IQR = sample_summary['3rd quartile'] - sample_summary['1st quartile']
        bin_width = (2 * IQR) / (number_of_samples ** (1/3))

        # Calculate the bin count, with a minimum of 5 bins
        bins = max((sample_summary['Maximum frequency'] -
                    sample_summary['Minimum frequency']) / bin_width, 5)

        sample_frequencies_ax = sns.distplot(sample_frequencies, kde=False,
                                             rug=True, bins=round(bins))
        sample_frequencies_ax.get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        sample_frequencies_ax.set_xlabel('Frequency per sample')
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.pdf'))
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.png'))
        plt.gcf().clear()

    feature_summary, feature_frequencies = _frequency_summary(
        table, axis=0)
    if number_of_features > 1:
        feature_frequencies_ax = sns.distplot(feature_frequencies, kde=False,
                                              rug=True)
        feature_frequencies_ax.set_xlabel('Frequency per feature')
        feature_frequencies_ax.set_xscale('log')
        feature_frequencies_ax.set_yscale('log')
        feature_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-frequencies.pdf'))
        feature_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-frequencies.png'))

    sample_summary_table = _format_html_table(
        sample_summary.apply('{:,}'.format).to_frame('Frequency'))
    feature_summary_table = _format_html_table(
        feature_summary.apply('{:,}'.format).to_frame('Frequency'))

    index = os.path.join(TEMPLATES, 'summarize_assets', 'index.html')
    context = {
        'number_of_samples': number_of_samples,
        'number_of_features': number_of_features,
        'total_frequencies': int(np.sum(sample_frequencies)),
        'sample_summary_table': sample_summary_table,
        'feature_summary_table': feature_summary_table,
    }

    sample_frequencies.sort_values(inplace=True, ascending=False)
    feature_frequencies.sort_values(inplace=True, ascending=False)
    sample_frequencies.to_csv(
        os.path.join(output_dir, 'sample-frequency-detail.csv'))
    feature_frequencies.to_csv(
        os.path.join(output_dir, 'feature-frequency-detail.csv'))

    feature_frequencies_table = _format_html_table(
        feature_frequencies.to_frame('Frequency'))
    overview_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'overview.html')
    sample_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'sample-frequency-detail.html')
    feature_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'feature-frequency-detail.html')

    context.update({'max_count': counts.max(),
                    'feature_frequencies_table': feature_frequencies_table})
    templates = [index, sample_frequency_template,
                 feature_frequency_template, overview_template]
    q2templates.render(templates, output_dir, context=context)

    shutil.copytree(os.path.join(TEMPLATES, 'summarize_assets', 'app'),
                    os.path.join(output_dir, 'app'))


def _format_html_table(df):
    table = df.to_html(classes="table table-striped table-hover")
    return table.replace('border="1"', 'border="0"')


def _frequencies(table, axis):
    return table.sum(axis=axis)


def _frequency_summary(table, axis='sample'):
    frequencies = _frequencies(table, axis=axis)

    summary = pd.Series([frequencies.min(), frequencies.quantile(0.25),
                         frequencies.median(), frequencies.quantile(0.75),
                         frequencies.max(), frequencies.mean()],
                        index=['Minimum frequency', '1st quartile',
                               'Median frequency', '3rd quartile',
                               'Maximum frequency', 'Mean frequency'])
    return summary, frequencies
