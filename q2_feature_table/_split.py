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

    indices = metadata_df.reset_index(
        ).groupby(metadata.name)[metadata_df.index.name].apply(list).to_dict()

    try:
        qiime2.sdk.util.validate_result_collection_keys(*indices.keys())
    except KeyError as e:
        raise KeyError(
            "One or more invalid metadata column values identified during "
            "feature table splitting. All metadata column values must be "
            "valid ResultCollection keys when used for splitting a feature "
            f"table. The original error message is as follows: {str(e)}")

    result = {}
    for group, sample_ids in indices.items():
        t = table.filter(sample_ids, axis='sample', inplace=False)
        if filter_empty_features:
            t.remove_empty(axis='observation', inplace=True)
        result[group] = t
    return result
