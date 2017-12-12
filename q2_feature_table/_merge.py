# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
import collections


def overlap_methods():
    return ('error_on_overlapping_sample', 'error_on_overlapping_feature',
            'sum')


def _get_overlapping(tables, axis):
    ids = collections.Counter()
    for table in tables:
        ids.update(table.ids(axis=axis))
    return {e for e, c in ids.items() if c > 1}


def merge(tables: biom.Table,
          overlap_method: str='error_on_overlapping_sample') -> biom.Table:
    if overlap_method == 'error_on_overlapping_sample':
        overlapping_ids = _get_overlapping(tables, 'sample')
        if len(overlapping_ids) > 0:
            raise ValueError('Same samples are present in provided tables: %s'
                             % ', '.join(overlapping_ids))
    elif overlap_method == 'error_on_overlapping_feature':
        overlapping_ids = _get_overlapping(tables, 'observation')
        if len(overlapping_ids) > 0:
            raise ValueError('Same features are present in provided tables: %s'
                             % ', '.join(overlapping_ids))
    elif overlap_method == 'sum':
        # This is the default behavior for biom.Table.merge
        pass
    else:
        raise ValueError('Invalid overlap method: %s. Please provide one of '
                         'the following methods: %s.' %
                         (overlap_method, ', '.join(overlap_methods())))
    tables = iter(tables)
    result = next(tables)  # There is always at least 1
    for table in tables:
        result = result.merge(table)

    return result


def _merge_feature_data(data: pd.Series) -> pd.Series:
    data = iter(data)
    result = next(data)  # There is always at least 1
    for d in data:
        result = result.combine_first(d)
    return result


def merge_seqs(data: pd.Series) -> pd.Series:
    return _merge_feature_data(data)


def merge_taxa(data: pd.Series) -> pd.Series:
    return _merge_feature_data(data)
