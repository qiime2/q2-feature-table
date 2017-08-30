# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import qiime2
import numpy as np
import pandas as pd


def _collapse_factory(function):
    def collapse_f(table, axis):
        # axis is always the transpose of the original collapse axis
        return np.array([function(x) for x in table.iter_data(axis=axis)])
    return collapse_f


_mode_lookup = {
    'sum': _collapse_factory(np.sum),
    'median-ceiling': _collapse_factory(lambda x: np.ceil(np.median(x))),
    'mean-ceiling': _collapse_factory(lambda x: np.ceil(np.mean(x)))
}


def _munge_metadata_category(mc, ids, axis):
    # TODO: centralize these ideas in MetadataCategory somehow
    series = mc.to_series()
    table_ids = set(ids)
    series_ids = set(series.index)

    # Check numeric
    if pd.api.types.is_numeric_dtype(pd.to_numeric(series, errors='ignore')):
        raise ValueError("Cannot group by a numeric metadata category.")

    # Check missing
    missing_ids = table_ids - series_ids
    if missing_ids:
        raise ValueError("All %s IDs in the feature table must be present in "
                         "the metadata category, but the following are "
                         "missing: %r" % (axis, missing_ids))

    # While preserving order, get rid of any IDs found only in the metadata
    series = series.drop(series_ids - table_ids)

    # Check for empty values only after filtering down to relevant IDs
    missing_values = series.isnull() | (series == '')
    if missing_values.any():
        missing = set(series[missing_values].index)
        raise ValueError("There are missing metadata category value(s) for %s "
                         "ID(s): %r" % (axis, missing))

    return series


def group(table: biom.Table, axis: str, metadata: qiime2.MetadataCategory,
          mode: str) -> biom.Table:
    if table.is_empty():
        raise ValueError("Cannot group an empty table.")

    if axis == 'feature':
        biom_axis = 'observation'
    else:
        biom_axis = axis

    series = _munge_metadata_category(metadata, table.ids(axis=biom_axis),
                                      axis)

    grouped_table = table.collapse(lambda axis_id, _: series.loc[axis_id],
                                   collapse_f=_mode_lookup[mode],
                                   axis=biom_axis,
                                   norm=False,
                                   include_collapsed_metadata=False)
    # Reorder axis by first unique appearance of each group value in metadata
    # (makes it stable for identity mappings and easier to test)
    return grouped_table.sort_order(series.unique(), axis=biom_axis)
