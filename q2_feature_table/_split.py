# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import qiime2


_valid_value_characters = set(
    "0123456789"
    "abcdefghijklmnopqrstuvwxyz"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "-_.")


def _check_valid_value(value):
    return len(set(value) - _valid_value_characters) != 0


def _check_valid_values(values):
    invalid_values = [v for v in values if _check_valid_value(v)]
    if len(invalid_values) > 0:
        raise ValueError(
            "Invalid value(s) identified during feature table splitting. "
            "Values in metadata columns used for splitting can only contain "
            "alphanumeric characters, dashes (-), dots (.), and underscores "
            "(_). The problematic value(s) are: "
            f"{', '.join(invalid_values)}")


def split(table: biom.Table,
          metadata: qiime2.CategoricalMetadataColumn,
          filter_empty_features: bool = True) -> biom.Table:
    metadata_df = metadata.drop_missing_values().to_dataframe()

    indices = metadata_df.reset_index(
        ).groupby(metadata.name)[metadata_df.index.name].apply(list).to_dict()

    _check_valid_values(indices.keys())

    result = {}
    for group, sample_ids in indices.items():
        t = table.filter(sample_ids, axis='sample', inplace=False)
        if filter_empty_features:
            t.remove_empty(axis='observation', inplace=True)
        result[group] = t
    return result
