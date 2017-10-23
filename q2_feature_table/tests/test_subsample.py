# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
from biom.table import Table

from q2_feature_table import subsample


class SubsampleTests(TestCase):

    def test_subsample_samples(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        a = subsample(t, 2, 'sample')
        self.assertEqual(a.shape, (2, 2))

        sample_ids = frozenset(a.ids(axis='sample'))
        self.assertIn(sample_ids, set([frozenset(['S1', 'S2']),
                                       frozenset(['S1', 'S3']),
                                       frozenset(['S2', 'S3'])]))
        self.assertEqual(set(a.ids(axis='observation')), set(['O1', 'O2']))

        for i in a.ids(axis='sample'):
            npt.assert_equal(t.data(i, axis='sample'),
                             a.data(i, axis='sample'))

    def test_subsample_features(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]).T,
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2'])
        a = subsample(t, 2, 'feature')
        self.assertEqual(a.shape, (2, 2))

        sample_ids = frozenset(a.ids(axis='observation'))
        self.assertIn(sample_ids, set([frozenset(['O1', 'O2']),
                                       frozenset(['O1', 'O3']),
                                       frozenset(['O2', 'O3'])]))
        self.assertEqual(set(a.ids(axis='sample')), set(['S1', 'S2']))

        for i in a.ids(axis='observation'):
            npt.assert_equal(t.data(i, axis='observation'),
                             a.data(i, axis='observation'))

    def test_subsample_samples_oversample(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]).T,
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2'])
        with self.assertRaisesRegex(ValueError, "depth exceeds"):
            subsample(t, 10, 'sample')

    def test_subsample_features_oversample(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]).T,
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2'])
        with self.assertRaisesRegex(ValueError, "depth exceeds"):
            subsample(t, 10, 'feature')

    def test_subsample_samples_empty(self):
        t = Table(np.array([[0, 0, 0], [0, 0, 0]]).T,
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2'])
        with self.assertRaisesRegex(ValueError, "contains no"):
            subsample(t, 2, 'sample')

    def test_subsample_features_empty(self):
        t = Table(np.array([[0, 0, 0], [0, 0, 0]]).T,
                  ['O1', 'O2', 'O3'],
                  ['S1', 'S2'])
        with self.assertRaisesRegex(ValueError, "contains no"):
            subsample(t, 2, 'feature')


if __name__ == "__main__":
    main()
