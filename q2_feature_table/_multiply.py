# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import warnings

import biom
import numpy as np
from q2_types.feature_table import (
    FeatureTable, PresenceAbsence, RelativeFrequency
)


def _multiply(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two biom tables."""
    # Subset table1 to only include observations present in table2's samples
    table2_sample_ids = set(table2.ids(axis="sample"))
    table1_obs_to_keep = [
        obs_id
        for obs_id in table1.ids(axis="observation")
        if obs_id in table2_sample_ids
    ]

    if not table1_obs_to_keep:
        raise ValueError(
            "No overlapping features found between table1 observations and "
            "table2 samples."
        )

    if len(table1_obs_to_keep) < len(table1.ids(axis="observation")):
        warnings.warn(
            "Removed "
            f"{len(table1.ids(axis='observation')) - len(table1_obs_to_keep)}"
            " feature(s) from table1 that had no matching samples in table2."
        )

    table1 = table1.filter(table1_obs_to_keep, axis="observation")

    # Reorder table2 samples to match table1 observations
    table2 = table2.sort_order(table1.ids(axis="observation"), axis="sample")

    # Perform sparse matrix multiplication.
    # In biom.Table, matrix_data is stored as observations x samples.
    # After transpose, both table1.matrix_data.T and table2.matrix_data.T are
    # shaped as samples x observations; the dot product is then transposed back
    # when constructing the result biom.Table so that the final table is
    # observations (from table2) x samples (from table1).
    result_matrix = table1.matrix_data.T.dot(table2.matrix_data.T)

    result_table = biom.Table(
        result_matrix.T,
        observation_ids=table2.ids(axis="observation"),
        sample_ids=table1.ids(axis="sample"),
    )

    return result_table


def _multiply_tables(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two feature tables."""
    result = _multiply(table1, table2)
    return result


def _multiply_tables_relative(
        table1: biom.Table, table2: biom.Table
) -> biom.Table:
    """Calculate dot product of two feature tables and convert to
    a relative frequency table."""
    result = _multiply(table1, table2)
    result.norm(axis="sample", inplace=True)
    return result


def _multiply_tables_pa(table1: biom.Table, table2: biom.Table) -> biom.Table:
    """Calculate dot product of two feature tables and convert to
    a presence-absence table."""
    result = _multiply(table1, table2)
    # Convert to presence-absence (1 if non-zero, 0 otherwise)
    result_data = result.matrix_data.copy()
    result_data.data = np.ones_like(result_data.data)
    result = biom.Table(
        result_data,
        observation_ids=result.ids(axis="observation"),
        sample_ids=result.ids(axis="sample"),
    )
    return result


def multiply_tables(ctx, table1, table2):
    """Calculate dot product of two feature tables."""
    if (
        table1.type <= FeatureTable[PresenceAbsence]
        or table2.type <= FeatureTable[PresenceAbsence]
    ):
        multiply = ctx.get_action("feature-table", "_multiply_tables_pa")
    elif (
        table1.type <= FeatureTable[RelativeFrequency]
        or table2.type <= FeatureTable[RelativeFrequency]
    ):
        multiply = ctx.get_action("feature-table", "_multiply_tables_relative")
    else:
        multiply = ctx.get_action("feature-table", "_multiply_tables")
    (result,) = multiply(table1, table2)
    return result
