# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import qiime2


def split(table: biom.Table,
          metadata: qiime2.CategoricalMetadataColumn,
          filter_empty_features: bool = True) -> biom.Table:
    metadata_df = metadata.drop_missing_values().to_dataframe()
    lookup = metadata_df[metadata.name].to_dict()

    def partition_f(i, m):
        return lookup.get(i)

    unique_grps = sorted(set(lookup.values()))
    try:
        qiime2.sdk.util.validate_result_collection_keys(*unique_grps)
    except KeyError as e:
        raise KeyError(
            "One or more invalid metadata column values identified during "
            "feature table splitting. All metadata column values must be "
            "valid ResultCollection keys when used for splitting a feature "
            f"table. The original error message is as follows: {str(e)}")

    result = {}
    for group, tab in table.partition(partition_f):
        if group is None:
            continue

        if filter_empty_features:
            tab.remove_empty(axis='observation', inplace=True)
        result[group] = tab
    return result
