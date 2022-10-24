# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources

import numpy as np
from biom import Table

from qiime2 import Artifact


rep_seqs_dada2_url = 'https://docs.qiime2.org/2022.8/data/tutorials/' \
                     'moving-pictures/rep-seqs-dada2.qza'
rep_seqs_deblur_url = 'https://docs.qiime2.org/2022.8/data/tutorials/' \
                      'moving-pictures/rep-seqs-deblur.qza'


def artifact_from_url(url):
    def factory():
        import tempfile
        import requests
        import qiime2

        data = requests.get(url)

        with tempfile.NamedTemporaryFile() as f:
            f.write(data.content)
            f.flush()
            result = qiime2.Artifact.load(f.name)

        return result
    return factory


def _get_data_from_tests(path):
    return pkg_resources.resource_filename('q2_feature_table.tests',
                                           os.path.join('data', path))


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
    feature_table1 = use.init_artifact('feature_table1', ft1_factory)
    feature_table2 = use.init_artifact('feature_table2', ft2_factory)

    merged_table, = use.action(
        use.UsageAction(plugin_id='feature_table',
                        action_id='merge'),
        use.UsageInputs(tables=[feature_table1, feature_table2]),
        use.UsageOutputNames(merged_table='merged_table'),
    )


def feature_table_merge_three_tables_example(use):
    feature_table1 = use.init_artifact('feature_table1', ft1_factory)
    feature_table2 = use.init_artifact('feature_table2', ft2_factory)
    feature_table3 = use.init_artifact('feature_table3', ft3_factory)

    merged_table, = use.action(
        use.UsageAction(plugin_id='feature_table',
                        action_id='merge'),
        use.UsageInputs(
            tables=[feature_table1, feature_table2, feature_table3],
            overlap_method='sum'
        ),
        use.UsageOutputNames(merged_table='merged_table'),
    )


def feature_table_merge_seqs(use):
    dada2_seqs = use.init_artifact_from_url('seqs1', rep_seqs_dada2_url)
    deblur_seqs = use.init_artifact_from_url('seqs2', rep_seqs_deblur_url)

    merged_data, = use.action(
        use.UsageAction('feature_table', 'merge_seqs'),
        use.UsageInputs(
            data=[dada2_seqs, deblur_seqs]
        ),
        use.UsageOutputNames(
            merged_data='merged_data'
        )
    )
