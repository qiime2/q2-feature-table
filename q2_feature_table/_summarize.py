# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

import biom
import numpy as np
import pandas as pd
import scipy.optimize
import seaborn as sns
import matplotlib.pyplot as plt
from q2_types.feature_data import DNAIterator

_blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                       "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                       "=nt&CMD=Put&QUERY=%s")


def view_seq_data(output_dir: str, data: DNAIterator) -> None:
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        # style adapted from http://stackoverflow.com/a/6667942/3424666
        fh.write(' <style>\n'
                 '  table {table-layout:fixed; width:1250px;}\n'
                 '  table td {word-wrap:break-word;}\n'
                 ' </style>\n')
        fh.write('To BLAST a sequence against the NCBI nt database, click the '
                 'BLAST link next to that sequence and then click the '
                 '<i>View report</i> button on the resulting page.<p>\n')
        fh.write(' <table border=1>\n')
        fh.write('  <tr><td width=300px>Feature ID</td>'
                 '<td width=100px>BLAST link</td>'
                 '<td>Sequence</td></tr>\n')
        for sequence in data:
            blast_url = _blast_url_template % (str(sequence))
            fh.write('  <tr><td>%s</td><td><a target="_blank" href="%s">'
                     'BLAST</a></td><td>%s</td></tr>\n' %
                     (sequence.metadata['id'], blast_url, str(sequence)))
        fh.write(' </table>\n')
        fh.write('</body></html>')


def view_taxa_data(output_dir: str, data: pd.Series) -> None:
    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        # style adapted from http://stackoverflow.com/a/6667942/3424666
        fh.write(' <style>\n'
                 '  table {table-layout:fixed; width:1250px;}\n'
                 '  table td {word-wrap:break-word;}\n'
                 ' </style>\n')
        # Note: not using to_html here as that can result in the
        # taxonomy strings being truncated.
        fh.write(' <table border=1>\n')
        fh.write('  <tr><td width=300px>Feature ID</td>'
                 '<td>Taxonomy</td></tr>\n')
        for id_, taxonomy in data.iteritems():
            fh.write('  <tr><td>%s</td><td>%s</td></tr>' % (id_, taxonomy))
        fh.write(' </table>\n')
        fh.write('</body></html>')


def summarize(output_dir: str, table: biom.Table) -> None:
    sample_summary, sample_counts = _count_summary(table, axis='sample')

    sample_md, feature_md = _metadata_summary(table)
    if len(sample_md) == 0:
        sample_md = 'No sample metadata.'
    else:
        sample_md = ', '.join(sample_md)
    if len(feature_md) == 0:
        feature_md = 'No feature metadata.'
    else:
        feature_md = ', '.join(feature_md)

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

    with open(os.path.join(output_dir, 'index.html'), 'w') as fh:
        fh.write('<html><body>\n')
        fh.write('<h1>Table summary</h1>\n')
        fh.write('Number of samples: %d<br>\n' % len(table.ids(axis='sample')))
        fh.write('Number of features: %d<br>\n' %
                 len(table.ids(axis='observation')))
        fh.write('Total counts: %d<br>\n' % np.sum(sample_counts))
        fh.write('An even sampling (i.e., rarefaction) depth of %d will retain'
                 ' the largest number of sequences.<br>\n' %
                 max_count_even_sampling_depth)
        fh.write('Sample metadata categories: %s<br>\n' % sample_md)
        fh.write('Feature metadata categories: %s<br>\n' % feature_md)

        fh.write('<h1>Count per sample summary</h1>\n')
        # This conversion is necessary to call to_html, pending
        # https://github.com/pydata/pandas/issues/5563
        fh.write('%s\n' % sample_summary.to_frame(name='Count').to_html())
        fh.write('<p>\n')
        fh.write('<img src="./sample-counts.png" width="700" height="500"><p>')
        fh.write('\n')
        fh.write('Image source (<a href="./sample-counts.png">png</a> | '
                 '<a href="./sample-counts.pdf">pdf</a>)<p>\n')
        fh.write('Counts per sample detail (<a href="./sample-count-detail.csv'
                 '">csv</a> | <a href="./sample-count-detail.html">html</a>)'
                 '<p>\n')

        fh.write('<h1>Count per feature summary</h1>\n')
        # This conversion is necessary to call to_html, pending
        # https://github.com/pydata/pandas/issues/5563
        fh.write('%s\n' % feature_summary.to_frame(name='Count').to_html())

        fh.write('<h1>Count per feature detail</h1>')
        fh.write('<img src="./feature-counts.png" width="700" height="500">'
                 '<p>\n')
        fh.write('Image source (<a href="./feature-counts.png">png</a> | '
                 '<a href="./feature-counts.pdf">pdf</a>)<p>\n')

        fh.write('</body></html>')

    sample_counts.sort_values(inplace=True)
    sample_counts.to_csv(os.path.join(output_dir, 'sample-count-detail.csv'))
    with open(os.path.join(output_dir, 'sample-count-detail.html'), 'w') as fh:
        fh.write('<html><body>\n')
        # This conversion is necessary to call to_html, pending
        # https://github.com/pydata/pandas/issues/5563
        fh.write(sample_counts.to_frame('Counts').to_html())
        fh.write('\n</body></html>')


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


def _metadata_summary(table):
    if table.metadata() is None:
        sample_md_keys = []
    else:
        sample_md_keys = table.metadata()[0].keys()

    if table.metadata(axis='observation') is None:
        feature_md_keys = []
    else:
        feature_md_keys = table.metadata(axis='observation')[0].keys()

    return list(sample_md_keys), list(feature_md_keys)


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
