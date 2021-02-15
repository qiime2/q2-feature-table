# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
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
            'average', 'sum')


def _get_overlapping(tables, axis):
    ids = collections.Counter()
    for table in tables:
        ids.update(table.ids(axis=axis))
    return {e for e, c in ids.items() if c > 1}


def merge(tables: biom.Table,
          overlap_method: str = 'error_on_overlapping_sample') -> biom.Table:
    if len(tables) == 1:
        return tables[0]

    if overlap_method == 'error_on_overlapping_sample':
        try:
            return tables[0].concat(tables[1:], 'sample')
        except biom.exception.DisjointIDError:
            overlapping = _get_overlapping(tables, 'sample')
            raise ValueError('Same samples are present in some of the '
                             'provided tables: %s' % ', '.join(overlapping))
    elif overlap_method == 'error_on_overlapping_feature':
        try:
            return tables[0].concat(tables[1:], 'observation')
        except biom.exception.DisjointIDError:
            overlapping = _get_overlapping(tables, 'observation')
            raise ValueError('Same features are present in some of the '
                             'provided tables: %s' % ', '.join(overlapping))
    elif overlap_method in ('sum', 'average'):
        n = len(tables)
        tables = iter(tables)
        result = next(tables)  # There is always at least 1
        for table in tables:
            result = result.merge(table)
        if overlap_method == 'average':
            result.transform(lambda v, _, __: v / n)
        return result
    else:
        raise ValueError('Invalid overlap method: %s. Please provide one of '
                         'the following methods: %s.' %
                         (overlap_method, ', '.join(overlap_methods())))


def _merge_feature_data(data):
    data = iter(data)
    result = next(data)  # There is always at least 1
    for d in data:
        result = result.combine_first(d)
    return result


def merge_seqs(data: pd.Series) -> pd.Series:
    return _merge_feature_data(data)


def merge_taxa(data: pd.DataFrame) -> pd.DataFrame:
    data = _merge_feature_data(data)
    # merge orders columns alphabetically; Taxon must be first header column
    # as defined here: https://github.com/qiime2/q2-types/blob/
    # 067d83e2aefe98674433e95162336fb5b9d96474/q2_types/feature_data/
    # _format.py#L97
    data = data[data.columns.drop('Taxon').insert(0, 'Taxon')]
    return data
