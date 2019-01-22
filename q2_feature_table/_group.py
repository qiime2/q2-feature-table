# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import qiime2
import numpy as np


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


def _munge_metadata_column(mc, ids, axis):
    mc = mc.filter_ids(ids)

    # Check for empty values only after filtering down to relevant IDs.
    missing = mc.get_ids(where_values_missing=True)
    if missing:
        raise ValueError("There are missing metadata column value(s) for "
                         "these %s ID(s): %s" %
                         (axis, ', '.join(repr(e) for e in sorted(missing))))
    return mc


def group(table: biom.Table, axis: str,
          metadata: qiime2.CategoricalMetadataColumn, mode: str) -> biom.Table:
    if table.is_empty():
        raise ValueError("Cannot group an empty table.")

    if axis == 'feature':
        biom_axis = 'observation'
    else:
        biom_axis = axis

    metadata = _munge_metadata_column(metadata, table.ids(axis=biom_axis),
                                      axis)

    grouped_table = table.collapse(
        lambda axis_id, _: metadata.get_value(axis_id),
        collapse_f=_mode_lookup[mode],
        axis=biom_axis,
        norm=False,
        include_collapsed_metadata=False)
    # Reorder axis by first unique appearance of each group value in metadata
    # (makes it stable for identity mappings and easier to test)
    # TODO use CategoricalMetadataColumn API for retrieving categories/groups,
    # when the API exists.
    series = metadata.to_series()
    return grouped_table.sort_order(series.unique(), axis=biom_axis)
