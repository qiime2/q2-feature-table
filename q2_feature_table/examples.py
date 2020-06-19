# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from biom import Table

from qiime2.sdk.usage import UsageAction, UsageInputs, UsageOutputNames
from qiime2 import Artifact


def ft1_factory():
    return Artifact.import_data(
        'FeatureTable[Frequency]',
        Table(np.array([[0, 1, 3], [1, 1, 2]]),
              ['O1', 'O2'],
              ['S1', 'S2', 'S3']))


def ft2_factory():
    return Artifact.import_data(
        'FeatureTable[Frequency]',
        Table(np.array([[0, 2, 6], [2, 2, 4]]),
              ['O1', 'O3'],
              ['S4', 'S5', 'S6']))


def ft3_factory():
    return Artifact.import_data(
        'FeatureTable[Frequency]',
        Table(np.array([[0, 4, 9], [4, 4, 8]]),
              ['O1', 'O4'],
              ['S7', 'S8', 'S9']))


def feature_table_merge_example(use):
    feature_table1 = use.init_data('feature_table1', ft1_factory)
    feature_table2 = use.init_data('feature_table2', ft2_factory)
    merged_table = use.init_data_collection('merged_table', list,
                                            feature_table1, feature_table2)

    use.action(
        UsageAction(plugin_id='feature_table',
                    action_id='merge'),
        UsageInputs(tables=merged_table),
        UsageOutputNames(merged_table='merged_table'),
    )


def feature_table_merge_three_tables_example(use):
    feature_table1 = use.init_data('feature_table1', ft1_factory)
    feature_table2 = use.init_data('feature_table2', ft2_factory)
    feature_table3 = use.init_data('feature_table3', ft3_factory)
    merged_table = use.init_data_collection('merged_table', list,
                                            feature_table1,
                                            feature_table2,
                                            feature_table3)

    use.action(
        UsageAction(plugin_id='feature_table',
                    action_id='merge'),
        UsageInputs(tables=merged_table, overlap_method='sum'),
        UsageOutputNames(merged_table='merged_table'),
    )
