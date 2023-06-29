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

    result = {}
    for group, sample_ids in indices.items():
        t = table.filter(sample_ids, axis='sample', inplace=False)
        if filter_empty_features:
            t.remove_empty(axis='observation', inplace=True)
        result[group] = t
    return result
