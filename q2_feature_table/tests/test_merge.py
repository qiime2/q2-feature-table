# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
from biom.table import Table

from q2_feature_table import merge_tables


class MergeTests(unittest.TestCase):

    def test_valid_overlapping_feature_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S4', 'S5', 'S6'])
        obs = merge_tables(t1, t2)
        exp = Table(np.array([[0, 1, 3, 0, 2, 6], [1, 1, 2, 0, 0, 0],
                              [0, 0, 0, 2, 2, 4]]),
                    ['O1', 'O2', 'O3'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_valid_non_overlapping_feature_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O3', 'O4'],
                   ['S4', 'S5', 'S6'])
        obs = merge_tables(t1, t2)
        exp = Table(np.array([[0, 1, 3, 0, 0, 0], [1, 1, 2, 0, 0, 0],
                              [0, 0, 0, 0, 2, 6], [0, 0, 0, 2, 2, 4]]),
                    ['O1', 'O2', 'O3', 'O4'],
                    ['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
        self.assertEqual(obs, exp)

    def test_overlapping_sample_ids(self):
        t1 = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                   ['O1', 'O2'],
                   ['S1', 'S2', 'S3'])
        t2 = Table(np.array([[0, 2, 6], [2, 2, 4]]),
                   ['O1', 'O3'],
                   ['S1', 'S5', 'S6'])
        with self.assertRaises(ValueError):
            merge_tables(t1, t2)

if __name__ == "__main__":
    unittest.main()
