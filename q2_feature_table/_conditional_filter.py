# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import biom


def filter_features_conditionally(table: biom.Table,
                                  prevelance: float, abundance: float
                                  ) -> biom.Table:
    """
    A function to perform joint filtering because it makes life better
    """
    num_samples, num_observations = table.shape
    prevelance = prevelance * num_samples

    # Calculates the filteering parameters on the original table
    def _filter_f(v, id_, md):
        return (v > abundance).sum() >= prevelance

    # Normalized the table to get the prevelance
    table_n = table.copy().norm(axis='sample', inplace=False)
    table_s = table_n.filter(_filter_f, axis='observation', inplace=True)
    filter_ids = table_s.ids(axis='observation')

    new_table = table.filter(filter_ids, axis='observation', inplace=False)

    return new_table
