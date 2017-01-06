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

import biom
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from q2_types.feature_data import DNAIterator
import q2templates
import skbio

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


def summarize(output_dir: str, table: biom.Table) -> None:
    number_of_samples = table.shape[1]
    number_of_features = table.shape[0]

    sample_summary, sample_frequencies = _frequency_summary(
        table, axis='sample')
    if number_of_samples > 1:
        sample_frequencies_ax = sns.distplot(sample_frequencies, kde=False,
                                             rug=True)
        sample_frequencies_ax.get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        sample_frequencies_ax.set_xlabel('Frequency per sample')
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.pdf'))
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.png'))
        plt.gcf().clear()

    feature_summary, feature_frequencies = _frequency_summary(
        table, axis='observation')
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

    sample_frequencies.sort_values(inplace=True)
    sample_frequencies.to_csv(
        os.path.join(output_dir, 'sample-frequency-detail.csv'))

    sample_frequencies_table = _format_html_table(
        sample_frequencies.to_frame('Frequency'))
    sample_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'sample-frequency-detail.html')

    context.update({'sample_frequencies_table': sample_frequencies_table})
    templates = [index, sample_frequency_template]
    q2templates.render(templates, output_dir, context=context)


def _format_html_table(df):
    table = df.to_html(classes="table table-striped table-hover")
    return table.replace('border="1"', 'border="0"')


def _frequencies(table, axis):
    return pd.Series(data=table.sum(axis=axis), index=table.ids(axis=axis))


def _frequency_summary(table, axis='sample'):
    frequencies = _frequencies(table, axis=axis)

    summary = pd.Series([frequencies.min(), frequencies.quantile(0.25),
                         frequencies.median(), frequencies.quantile(0.75),
                         frequencies.max()],
                        index=['Minimum frequency', '1st quartile',
                               'Median frequency', '3rd quartile',
                               'Maximum frequency'])
    mean = pd.Series([frequencies.mean()], index=['Mean frequency'])
    summary.sort_values(ascending=False, inplace=True)
    summary = mean.append(summary)
    return summary, frequencies
