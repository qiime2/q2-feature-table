# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd


def overlap_methods():
    return ('error_on_overlapping_sample', 'error_on_overlapping_feature')


def _get_overlapping(table1, table2, axis):
    table1_ids = set(table1.ids(axis=axis))
    table2_ids = set(table2.ids(axis=axis))
    return table1_ids & table2_ids


def merge(table1: biom.Table, table2: biom.Table,
          overlap_method: str='error_on_overlapping_sample') -> biom.Table:

    if overlap_method == 'error_on_overlapping_sample':
        overlapping_ids = _get_overlapping(table1, table2, 'sample')
        if len(overlapping_ids) > 0:
            raise ValueError('Some samples are present in both tables: %s' %
                             ', '.join(overlapping_ids))
    elif overlap_method == 'error_on_overlapping_feature':
        overlapping_ids = _get_overlapping(table1, table2, 'observation')
        if len(overlapping_ids) > 0:
            raise ValueError('Some features are present in both tables: %s' %
                             ', '.join(overlapping_ids))
    else:
        raise ValueError('Invalid overlap method: %s. Please provide one of '
                         'the following methods: %s.' %
                         (overlap_method, ', '.join(overlap_methods())))

    return table1.merge(table2)


def _merge_feature_data(data1: pd.Series, data2: pd.Series) \
        -> pd.Series:
    return data1.combine_first(data2)


def merge_seq_data(data1: pd.Series, data2: pd.Series) \
        -> pd.Series:
    return _merge_feature_data(data1, data2)


def merge_taxa_data(data1: pd.Series, data2: pd.Series) \
        -> pd.Series:
    return _merge_feature_data(data1, data2)
