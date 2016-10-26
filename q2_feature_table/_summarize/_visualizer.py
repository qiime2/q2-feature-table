# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import biom
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from q2_types.feature_data import DNAIterator
import q2templates

_blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                       "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                       "=nt&CMD=Put&QUERY=%s")

TEMPLATES = pkg_resources.resource_filename('q2_feature_table', '_summarize')


def view_seq_data(output_dir: str, data: DNAIterator) -> None:
    sequences = []
    for sequence in data:
        str_seq = str(sequence)
        sequences.append({'id': sequence.metadata['id'],
                          'url': _blast_url_template % str_seq,
                          'seq': str_seq})

    index = os.path.join(TEMPLATES, 'view_seq_data_assets', 'index.html')
    q2templates.render(index, output_dir, context={'data': sequences})


def view_taxa_data(output_dir: str, data: pd.Series) -> None:
    prepped = []
    for _id, taxa in data.iteritems():
        prepped.append({'id': _id, 'taxa': taxa})

    index = os.path.join(TEMPLATES, 'view_taxa_data_assets', 'index.html')
    q2templates.render(index, output_dir, context={'data': prepped})


def summarize(output_dir: str, table: biom.Table) -> None:
    number_of_samples = len(table.ids(axis='sample'))
    number_of_features = len(table.ids(axis='observation'))

    sample_summary, sample_counts = _count_summary(table, axis='sample')
    if number_of_samples > 1:
        sample_counts_ax = sns.distplot(sample_counts, kde=False, rug=True)
        sample_counts_ax.set_xlabel('Frequency of sample')
        sample_counts_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-counts.pdf'))
        sample_counts_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-counts.png'))
        plt.gcf().clear()

    feature_summary, feature_counts = _count_summary(table, axis='observation')
    if number_of_features > 1:
        feature_counts_ax = sns.distplot(feature_counts, kde=False, rug=True)
        feature_counts_ax.set_xlabel('Frequency of feature')
        feature_counts_ax.set_xscale('log')
        feature_counts_ax.set_yscale('log')
        feature_counts_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-counts.pdf'))
        feature_counts_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-counts.png'))

    sample_summary_table = _format_html_table(
        sample_summary.to_frame('Frequency'))
    feature_summary_table = _format_html_table(
        feature_summary.to_frame('Frequency'))

    index = os.path.join(TEMPLATES, 'summarize_assets', 'index.html')
    context = {
        'number_of_samples': number_of_samples,
        'number_of_features': number_of_features,
        'total_counts': int(np.sum(sample_counts)),
        'sample_summary_table': sample_summary_table,
        'feature_summary_table': feature_summary_table,
    }

    sample_counts.sort_values(inplace=True)
    sample_counts.to_csv(os.path.join(output_dir, 'sample-count-detail.csv'))

    sample_counts_table = _format_html_table(sample_counts.to_frame('Counts'))
    sample_count_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'sample-count-detail.html')

    context.update({'sample_counts_table': sample_counts_table})
    templates = [index, sample_count_template]
    q2templates.render(templates, output_dir, context=context)


def _format_html_table(df):
    table = df.to_html(classes="table table-striped table-hover")
    return table.replace('border="1"', 'border="0"')


def _counts(table, axis):
    result = {}
    for count_vector, id_, _ in table.iter(axis=axis):
        result[id_] = float(count_vector.sum())
    return pd.Series(result)


def _count_summary(table, axis='sample'):
    counts = _counts(table, axis=axis)

    summary = pd.Series([counts.min(), counts.quantile(0.25), counts.median(),
                         counts.quantile(0.75), counts.max(), counts.mean()],
                        index=['Minimum frequency', '1st quartile',
                               'Median frequency', '3rd quartile',
                               'Maximum frequency', 'Mean frequency'])
    summary.sort_values()
    return summary, counts
