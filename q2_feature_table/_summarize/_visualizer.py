# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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
import qiime2
import json
from ._vega_spec import vega_spec

_blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                       "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                       "=nt&CMD=Put&QUERY=%s")

TEMPLATES = pkg_resources.resource_filename('q2_feature_table', '_summarize')


def tabulate_seqs(output_dir: str, data: DNAIterator,
                  taxonomy: pd.DataFrame = None,
                  metadata: qiime2.Metadata = None,
                  merge_method: str = 'strict') -> None:

    display_sequences = set()
    sequences = {}
    seq_lengths = []
    with open(os.path.join(output_dir, 'sequences.fasta'), 'w') as fh:
        for sequence in data:
            skbio.io.write(sequence, format='fasta', into=fh)
            str_seq = str(sequence)
            seq_len = len(str_seq)
            display_sequences.add(sequence.metadata['id'])
            sequences[sequence.metadata['id']]\
                = {'len': seq_len,
                   'url': _blast_url_template % str_seq,
                   'seq': str_seq}
            seq_lengths.append(seq_len)

    if metadata is not None:
        metadata_df = metadata.to_dataframe()
        if merge_method == 'union':
            display_sequences = display_sequences.union(metadata_df.index)
        elif merge_method == 'intersect':
            display_sequences = display_sequences.intersection(
                metadata_df.index)
        elif merge_method == 'strict':
            if set(metadata_df.index) != display_sequences:
                raise ValueError('Merge method is strict and IDs do not match')
    if taxonomy is not None:
        for member in taxonomy.values():
            if merge_method == 'union':
                display_sequences = display_sequences.union(member.index)
            elif merge_method == 'intersect':
                display_sequences = display_sequences.intersection(
                    member.index)
            elif merge_method == 'strict':
                if set(member.index) != display_sequences:
                    raise ValueError(
                                'Merge method is strict and IDs do not match')

    seq_len_stats = _compute_descriptive_stats(seq_lengths)
    _write_tsvs_of_descriptive_stats(seq_len_stats, output_dir)

    index = os.path.join(TEMPLATES, 'tabulate_seqs_assets', 'index.html')
    context = {'data': sequences, 'stats': seq_len_stats}
    if taxonomy is not None:
        context['taxonomy'] = taxonomy
    if metadata is not None:
        context['metadata'] = metadata_df
    context['display_sequences'] = display_sequences
    q2templates.render(index, output_dir, context=context)

    js = os.path.join(
        TEMPLATES, 'tabulate_seqs_assets', 'js', 'tsorter.min.js')
    os.mkdir(os.path.join(output_dir, 'js'))
    shutil.copy(js, os.path.join(output_dir, 'js', 'tsorter.min.js'))


def summarize(output_dir: str, table: biom.Table,
              sample_metadata: qiime2.Metadata = None) -> None:
    # this value is to limit the amount of memory used by seaborn.histplot, for
    # more information see: https://github.com/mwaskom/seaborn/issues/2325
    MAX_BINS = 50
    number_of_features, number_of_samples = table.shape

    sample_summary, sample_frequencies = _frequency_summary(
        table, axis='sample')

    if number_of_samples > 1:

        # Calculate the bin count, with a minimum of 5 bins
        IQR = sample_summary['3rd quartile'] - sample_summary['1st quartile']
        if IQR == 0.0:
            bins = 5
        else:
            # Freedman–Diaconis rule
            bin_width = (2 * IQR) / (number_of_samples ** (1/3))

            bins = max((sample_summary['Maximum frequency'] -
                        sample_summary['Minimum frequency']) / bin_width, 5)
            bins = min(bins, MAX_BINS)

        sample_frequencies_ax = sns.histplot(sample_frequencies, kde=False,
                                             bins=int(round(bins)))
        sample_frequencies_ax.get_xaxis().set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
        sample_frequencies_ax.set_xlabel('Frequency per sample')
        sample_frequencies_ax.set_ylabel('Number of samples')
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.pdf'))
        sample_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'sample-frequencies.png'))
        plt.gcf().clear()

    feature_summary, feature_frequencies = _frequency_summary(
        table, axis='observation')
    if number_of_features > 1:
        IQR = feature_summary['3rd quartile'] - feature_summary['1st quartile']
        if IQR == 0.0:
            bins = 5
        else:
            # Freedman–Diaconis rule
            bin_width = (2 * IQR) / (number_of_features ** (1/3))

            bins = max((feature_summary['Maximum frequency'] -
                        feature_summary['Minimum frequency']) / bin_width, 5)
            bins = min(bins, MAX_BINS)

        feature_frequencies_ax = sns.histplot(feature_frequencies, kde=False,
                                              bins=int(round(bins)))
        feature_frequencies_ax.set_xlabel('Frequency per feature')
        feature_frequencies_ax.set_ylabel('Number of features')
        feature_frequencies_ax.set_xscale('log')
        feature_frequencies_ax.set_yscale('log')
        feature_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-frequencies.pdf'))
        feature_frequencies_ax.get_figure().savefig(
            os.path.join(output_dir, 'feature-frequencies.png'))

    sample_summary_json = pd.DataFrame(
        sample_summary, columns=['Frequency']).to_json()
    feature_summary_json = pd.DataFrame(
        feature_summary, columns=['Frequency']).to_json()

    index = os.path.join(TEMPLATES, 'summarize_assets', 'index.html')
    context = {
        'number_of_samples': number_of_samples,
        'number_of_features': number_of_features,
        'total_frequencies': int(np.sum(sample_frequencies)),
        'sample_summary_table': sample_summary_json,
        'feature_summary_table': feature_summary_json,
    }

    # Create a JSON object containing the Sample Frequencies to build the
    # table in sample-frequency-detail.html
    #
    # Cast to DataFrame to standardize with other tables
    sample_frequencies_json = pd.DataFrame(
        sample_frequencies, columns=['Frequency']).to_json()

    # Create a JSON object containing the Feature Frequencies to build the
    # table in feature-frequency-detail.html
    feature_qualitative_data = _compute_qualitative_summary(table)
    feature_frequencies = feature_frequencies.astype(int).to_frame('Frequency')
    feature_frequencies['# of Samples Observed In'] = \
        pd.Series(feature_qualitative_data).astype(int)
    feature_frequencies_json = feature_frequencies.to_json()

    sample_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'sample-frequency-detail.html')
    feature_frequency_template = os.path.join(
        TEMPLATES, 'summarize_assets', 'feature-frequency-detail.html')

    context.update({'max_count': sample_frequencies.max(),
                    'feature_frequencies_json': feature_frequencies_json,
                    'feature_qualitative_data': feature_qualitative_data,
                    'tabs': [{'url': 'index.html',
                              'title': 'Overview'},
                             {'url': 'sample-frequency-detail.html',
                              'title': 'Interactive Sample Detail'},
                             {'url': 'feature-frequency-detail.html',
                              'title': 'Feature Detail'}]})

    templates = [index, sample_frequency_template, feature_frequency_template]
    context.update({'frequencies_list':
                    json.dumps(sorted(sample_frequencies.values.tolist()))})
    if sample_metadata is not None:
        context.update({'vega_spec':
                        json.dumps(vega_spec(sample_metadata,
                                             sample_frequencies
                                             ))
                        })
    context.update({'sample_frequencies_json': sample_frequencies_json})
    q2templates.util.copy_assets(os.path.join(TEMPLATES,
                                              'summarize_assets',
                                              'vega'),
                                 output_dir)
    q2templates.util.copy_assets(os.path.join(TEMPLATES,
                                              'summarize_assets',
                                              'utils'),
                                 output_dir)
    q2templates.render(templates, output_dir, context=context)

    plt.close('all')


def tabulate_feature_frequencies(table: biom.Table) -> qiime2.Metadata:
    feature_frequencies = _frequencies(table, 'observation')
    feature_frequencies = feature_frequencies.apply(
        '{:,}'.format).to_frame('Frequency')
    feature_qualitative_data = _compute_qualitative_summary(table)
    samples_observed_in =\
        pd.Series(feature_qualitative_data).astype(int).apply('{:,}'.format)
    feature_frequencies["No. of Samples Observed In"] = samples_observed_in
    feature_frequencies.index.name = "Feature ID"
    return qiime2.Metadata(feature_frequencies)


def tabulate_sample_frequencies(table: biom.Table) -> qiime2.Metadata:
    sample_frequencies = _frequencies(table, 'sample')
    sample_frequencies = sample_frequencies.apply(
        '{:,}'.format).to_frame('Frequency')
    sample_qualitative_data = _compute_qualitative_summary_sample(table)
    samples_with_feature =\
        pd.Series(sample_qualitative_data).astype(int).apply('{:,}'.format)
    sample_frequencies["No. of Associated Features"] = samples_with_feature
    sample_frequencies.index.name = "Sample ID"
    return qiime2.Metadata(sample_frequencies)


def summarize_plus(ctx, table, metadata=None):

    try:
        table_dimensions = table.view(pd.DataFrame).shape

    except ValueError:
        raise ValueError('Cannot summarize a table with no features')

    if table_dimensions[0] == 0:
        raise ValueError('Cannot summarize a table with no samples')

    _feature_frequencies = ctx.get_action('feature_table',
                                          'tabulate_feature_frequencies')
    _sample_frequencies = ctx.get_action('feature_table',
                                         'tabulate_sample_frequencies')
    _visualizer = ctx.get_action('feature_table',
                                 'summarize')

    feature_frequencies, = _feature_frequencies(table)
    sample_frequencies, = _sample_frequencies(table)
    summary, = _visualizer(table, metadata)

    return feature_frequencies, sample_frequencies, summary


def _compute_descriptive_stats(lst: list):
    """Basic descriptive statistics and a (parametric) seven-number summary.

    Calculates descriptive statistics for a list of numerical values, including
    count, min, max, mean, and a parametric seven-number-summary. This summary
    includes values for the lower quartile, median, upper quartile, and
    percentiles 2, 9, 91, and 98. If the data is normally distributed, these
    seven percentiles will be equally spaced when plotted.

    Parameters
    ----------
    lst : list of int or float values

    Returns
    -------
    dict
        a dictionary containing the following descriptive statistics:

        count
            int: the number of items in `lst`
        min
            int or float: the smallest number in `lst`
        max
            int or float: the largest number in `lst`
        mean
            float: the mean of `lst`
        range
            int or float: the range of values in `lst`
        std
            float: the standard deviation of values in `lst`
        seven_num_summ_percentiles
            list of floats: the parameter percentiles used to calculate this
            seven-number summary: [2, 9, 25, 50, 75, 91, 98]
        seven_num_summ_values
            list of floats: the calculated percentile values of the summary

    """
    # NOTE: With .describe(), NaN values in passed lst are excluded by default
    if len(lst) == 0:
        raise ValueError('No values provided.')

    seq_lengths = pd.Series(lst)
    seven_num_summ_percentiles = [0.02, 0.09, 0.25, 0.5, 0.75, 0.91, 0.98]
    descriptive_stats = seq_lengths.describe(
        percentiles=seven_num_summ_percentiles)

    return {'count': int(descriptive_stats.loc['count']),
            'min': descriptive_stats.loc['min'],
            'max': descriptive_stats.loc['max'],
            'range': descriptive_stats.loc['max'] -
            descriptive_stats.loc['min'],
            'mean': descriptive_stats.loc['mean'],
            'std': descriptive_stats.loc['std'],
            'seven_num_summ_percentiles': seven_num_summ_percentiles,
            'seven_num_summ_values': descriptive_stats.loc['2%':'98%'].tolist()
            }


def _write_tsvs_of_descriptive_stats(dictionary: dict, output_dir: str):
    descriptive_stats = ['count', 'min', 'max', 'mean', 'range', 'std']
    stat_list = []
    for key in descriptive_stats:
        stat_list.append(dictionary[key])
    descriptive_stats = pd.DataFrame(
        {'Statistic': descriptive_stats, 'Value': stat_list})
    descriptive_stats.to_csv(
        os.path.join(output_dir, 'descriptive_stats.tsv'),
        sep='\t', index=False, float_format='%g')

    seven_number_summary = pd.DataFrame(
        {'Quantile': dictionary['seven_num_summ_percentiles'],
         'Value': dictionary['seven_num_summ_values']})
    seven_number_summary.to_csv(
        os.path.join(output_dir, 'seven_number_summary.tsv'),
        sep='\t', index=False, float_format='%g')


def _compute_qualitative_summary(table):
    table = table.transpose()
    sample_count = {}
    for count_vector, feature_id, _ in table.iter():
        sample_count[feature_id] = (count_vector != 0).sum()
    return sample_count


def _compute_qualitative_summary_sample(table):
    feature_count = {}
    for count_vector, sample_id, _ in table.iter():
        feature_count[sample_id] = (count_vector != 0).sum()
    return feature_count


def _frequencies(table, axis):
    return pd.Series(data=table.sum(axis=axis), index=table.ids(axis=axis))


def _frequency_summary(table: biom.Table, axis='sample'):
    frequencies = _frequencies(table, axis=axis)

    first = frequencies.quantile(0.25)
    third = frequencies.quantile(0.75)
    _fst = round(first, 1)
    _min = round(frequencies.min(), 1)
    _thd = round(third, 1)
    _med = round(frequencies.median(), 1)
    _max = round(frequencies.max(), 1)
    mean = round(frequencies.mean(), 1)

    summary = pd.Series([_min, _fst,
                         _med, _thd,
                         _max, mean],
                        index=['Minimum frequency', '1st quartile',
                               'Median frequency', '3rd quartile',
                               'Maximum frequency', 'Mean frequency'])
    return summary, frequencies
