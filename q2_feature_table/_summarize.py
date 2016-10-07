# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import pkg_resources

import biom
import numpy as np
import pandas as pd
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
from q2_types.feature_data import DNAIterator
from trender import TRender

_blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                       "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                       "=nt&CMD=Put&QUERY=%s")


def view_seq_data(output_dir: str, data: DNAIterator) -> None:
    sequences = []
    for sequence in data:
        str_seq = str(sequence)
        sequences.append({'id': sequence.metadata['id'],
                          'url': _blast_url_template % str_seq,
                          'seq': str_seq})

    TEMPLATES = pkg_resources.resource_filename('q2_feature_table', 'assets')
    index = TRender('view_seq_data.template', path=TEMPLATES)
    rendered_index = index.render({'data': sequences})

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(rendered_index)

    for fn in ['bootstrap.min.css', 'qiime_logo_large.png']:
        shutil.copy(os.path.join(TEMPLATES, fn), os.path.join(output_dir, fn))


def view_taxa_data(output_dir: str, data: pd.Series) -> None:
    prepped = []
    for _id, taxa in data.iteritems():
        prepped.append({'id': _id, 'taxa': taxa})

    TEMPLATES = pkg_resources.resource_filename('q2_feature_table', 'assets')
    index = TRender('view_taxa_data.template', path=TEMPLATES)
    rendered_index = index.render({'data': prepped})

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(rendered_index)

    for fn in ['bootstrap.min.css', 'qiime_logo_large.png']:
        shutil.copy(os.path.join(TEMPLATES, fn), os.path.join(output_dir, fn))


def summarize(output_dir: str, table: biom.Table) -> None:
    sample_summary, sample_counts = _count_summary(table, axis='sample')

    max_count_even_sampling_depth = \
        _get_max_count_even_sampling_depth(sample_counts)
    sample_counts_ax = sns.distplot(sample_counts, kde=False, rug=True)
    sample_counts_ax.set_title('Counts per sample')
    sample_counts_ax.set_xlabel('Counts')
    sample_counts_ax.set_ylabel('Frequency')
    sample_counts_ax.get_figure().savefig(
        os.path.join(output_dir, 'sample-counts.pdf'))
    sample_counts_ax.get_figure().savefig(
        os.path.join(output_dir, 'sample-counts.png'))

    plt.gcf().clear()
    feature_summary, feature_counts = _count_summary(table, axis='observation')
    feature_counts_ax = sns.distplot(feature_counts, kde=False, rug=True)
    feature_counts_ax.set_title('Counts per feature')
    feature_counts_ax.set_xlabel('Counts')
    feature_counts_ax.set_ylabel('Frequency')
    feature_counts_ax.set_xscale('log')
    feature_counts_ax.set_yscale('log')
    feature_counts_ax.get_figure().savefig(
        os.path.join(output_dir, 'feature-counts.pdf'))
    feature_counts_ax.get_figure().savefig(
        os.path.join(output_dir, 'feature-counts.png'))

    sample_summary_table = _format_html_table(sample_summary.to_frame('Count'))
    feature_summary_table = _format_html_table(
        feature_summary.to_frame('Count'))

    TEMPLATES = pkg_resources.resource_filename('q2_feature_table', 'assets')
    index = TRender('summarize.template', path=TEMPLATES)
    rendered_index = index.render({
        'number_of_samples': len(table.ids(axis='sample')),
        'number_of_features': len(table.ids(axis='observation')),
        'total_counts': int(np.sum(sample_counts)),
        'max_count_even_sampling_depth': max_count_even_sampling_depth,
        'sample_summary_table': sample_summary_table,
        'feature_summary_table': feature_summary_table,
    })

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write(rendered_index)

    for fn in ['bootstrap.min.css', 'qiime_logo_large.png']:
        shutil.copy(os.path.join(TEMPLATES, fn), os.path.join(output_dir, fn))

    sample_counts.sort_values(inplace=True)
    sample_counts.to_csv(os.path.join(output_dir, 'sample-count-detail.csv'))

    sample_counts_table = _format_html_table(sample_counts.to_frame('Counts'))
    sample_count_template = TRender('sample-count-detail.template',
                                    path=TEMPLATES)
    rendered_sample_count_template = sample_count_template.render({
        'sample_counts_table': sample_counts_table})

    with open(os.path.join(output_dir, 'sample-count-detail.html'), 'w') as fh:
        fh.write(rendered_sample_count_template)


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
                        index=['Minimum count', '1st quartile', 'Median count',
                               '3rd quartile', 'Maximum count', 'Mean count'])
    summary.sort_values()
    return summary, counts


def _get_max_count_even_sampling_depth(counts):
    return _get_depth_for_max_sequence_count(counts)


def _get_depth_for_max_sequence_count(counts):
    """Find the even sampling depth that retains the most sequences."""
    def f(d):
        return -1 * _summarize_even_sampling_depth(d, counts)[2]

    res = scipy.optimize.minimize_scalar(
            f, bounds=(counts.min(), counts.max()), method='bounded')
    return int(np.floor(res.x))


def _summarize_even_sampling_depth(even_sampling_depth, counts):
    samples_retained = (counts >= even_sampling_depth)
    num_samples_retained = samples_retained.sum()
    num_sequences_retained = num_samples_retained * even_sampling_depth
    return samples_retained, num_samples_retained, num_sequences_retained
