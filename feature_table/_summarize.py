# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import scipy as sp
import scipy.optimize
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
    summary.sort_values()
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
