# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import scipy as sp
import biom
import pandas as pd


def _counts(table, axis):
    result = {}
    for count_vector, id_, _ in table.iter(axis=axis):
        result[id_] = float(count_vector.sum())
    return pd.Series(result)

def count_summary(table, axis='sample'):
    counts = _counts(table, axis=axis)

    summary = pd.Series([counts.min(), counts.quantile(0.25), counts.median(),
                         counts.quantile(0.75), counts.max(), counts.mean()],
                        index=['Minimum count', '1st quartile', 'Median count', '3rd quartile',
                               'Maximum count', 'Mean count'])
    summary.sort()
    return summary, counts

def metadata_summary(table):
    if table.metadata() is None:
        sample_md_keys = []
    else:
        sample_md_keys = table.metadata()[0].keys()

    if table.metadata(axis='observation') is None:
        feature_md_keys = []
    else:
        feature_md_keys = table.metadata(axis='observation')[0].keys()

    return list(sample_md_keys), list(feature_md_keys)

def max_count_even_sampling_depth(counts):
    return _get_depth_for_max_sequence_count(counts)

def _summarize_even_sampling_depth(even_sampling_depth, counts):
    samples_retained = (counts >= even_sampling_depth)
    num_samples_retained = samples_retained.sum()
    num_sequences_retained = num_samples_retained * even_sampling_depth
    return samples_retained, num_samples_retained, num_sequences_retained

def _get_depth_for_max_sequence_count(counts):
    """Find the even sampling depth that retains the most sequences."""
    def f(d):
        return -1 * _summarize_even_sampling_depth(d, counts)[2]

    res = sp.optimize.minimize_scalar(f,
                          bounds=(counts.min(), counts.max()),
                          method='bounded')
    return int(np.floor(res.x))


## old...

import ipywidgets
from ipywidgets import interactive, fixed, IntSlider
from IPython.display import display

def explore_sampling_depth(biom):
    import seaborn as sns
    counts = biom.sum()
    count_summary = counts.describe()
    total_num_samples = len(counts)
    total_num_sequences = counts.sum()
    depth_for_max_sequence_count = _get_depth_for_max_sequence_count(counts)
    sampling_depth_slider = IntSlider(min=count_summary['min'],
                                      max=count_summary['max'],
                                      step=10 ** (math.log(count_summary['max'], 10) - 2),
                                      value=depth_for_max_sequence_count)
    default_samples_retained, default_num_samples_retained, default_num_sequences_retained = \
            _summarize_even_sampling_depth(depth_for_max_sequence_count, counts)

    default_percent_samples_retained = default_num_samples_retained * 100 / total_num_samples
    default_percent_sequences_retained = default_num_sequences_retained * 100 / total_num_sequences

    label_s = "Depth {0}: {1:.2f}% of sequences and {2:.2f}% of samples retained."

    def f(even_sampling_depth):
        samples_retained, num_samples_retained, num_sequences_retained = \
            _summarize_even_sampling_depth(even_sampling_depth, counts)
        percent_samples_retained = num_samples_retained * 100 / total_num_samples
        percent_sequences_retained = num_sequences_retained * 100 / total_num_sequences
        ax = sns.distplot(counts)
        ax.set_xlabel("Number of sequences per sample")
        ax.set_ylabel("Frequency")
        line_label = label_s.format(depth_for_max_sequence_count,
                                    default_percent_sequences_retained,
                                    default_percent_samples_retained)
        ax.plot([depth_for_max_sequence_count, depth_for_max_sequence_count], ax.get_ylim(),
                'k--', label=line_label)

        line_label = label_s.format(even_sampling_depth,
                                    percent_sequences_retained,
                                    percent_samples_retained)
        ax.plot([even_sampling_depth, even_sampling_depth], ax.get_ylim(),
                'k-', label=line_label)
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    def reset_depth(_):
        sampling_depth_slider.value = depth_for_max_sequence_count

    reset = ipywidgets.Button(icon='fa-refresh')
    reset.on_click(reset_depth)

    w = interactive(f, even_sampling_depth=sampling_depth_slider)
    display(ipywidgets.HBox(children=[w, reset]))
