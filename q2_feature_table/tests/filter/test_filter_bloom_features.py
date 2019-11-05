# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2
import numpy as np
import pandas as pd
from biom.table import Table

from q2_feature_table import filter_bloom_features


class FilterBloomFeaturesTests(unittest.TestCase):
    def setUp(self):
        # self.seqs = pd.Series(['ACGT', 'GCTA', 'CCCC', 'TGTT'],
        #                       index=['O1', 'O2', 'O3', 'O4'])
        self.seqs = pd.Series(['ACGT', 'GCTA'],
                              index=['O1', 'O2'])
        self.drop_seqs = pd.Series(['GCTAT', 'TGACGAC'])

    def test_filter_blooms(self):
        table = Table(np.array([[0, 1, 1], [1, 1, 2]]),
                      ['O1', 'O2'],
                      ['S1', 'S2', 'S3'])
        actual = filter_bloom_features(table, self.seqs,
                                       bloom_sequences=self.drop_seqs)
        expected = Table(np.array([[0, 1, 1]]),
                         ['O1'],
                         ['S1', 'S2', 'S3'])
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
