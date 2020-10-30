# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import biom


def filter_features_conditionally(table: biom.Table,
                                  prevalance: float, abundance: float
                                  ) -> biom.Table:
    """
    A function to perform joint filtering because it makes life better
    """
    num_observations, num_samples = table.shape
    prevalence = prevalence * num_samples

    # Calculates the filteering parameters on the original table
    def _filter_f(v, id_, md):
        return (v > abundance).sum() >= prevalence

    # Normalized the table to get the prevalance
    # Copy is because biom really wants to normalize the original table. By
    # copying and not using inplace, the original table is preserved.
    # Redundant, but better safe that sorry.
    table_norm = table.copy().norm(axis='sample', inplace=False)
    table_norm.filter(_filter_f, axis='observation', inplace=True)
    filter_ids = table_norm.ids(axis='observation')

    new_table = table.filter(filter_ids, axis='observation', inplace=False)

    return new_table
