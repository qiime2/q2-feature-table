# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
from biom.table import Table

from q2_feature_table import relative_frequency, presence_absence


class RelativeFrequencyTests(TestCase):

    def test_relative_frequency(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        a = relative_frequency(t)
        self.assertEqual(a.shape, (2, 3))
        self.assertEqual(set(a.ids(axis='sample')), set(['S1', 'S2', 'S3']))
        self.assertEqual(set(a.ids(axis='observation')), set(['O1', 'O2']))
        npt.assert_array_equal(a.sum(axis='sample'), np.array([1., 1., 1.]))
        npt.assert_array_equal(a.matrix_data.toarray(),
                               np.array([[0, 0.5, 3/5], [1.0, 0.5, 2/5]]))


class PresenceAbsenceTests(TestCase):

    def test_presence_absence(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        a = presence_absence(t)
        self.assertEqual(a.shape, (2, 3))
        self.assertEqual(set(a.ids(axis='sample')), set(['S1', 'S2', 'S3']))
        self.assertEqual(set(a.ids(axis='observation')), set(['O1', 'O2']))
        npt.assert_array_equal(a.matrix_data.toarray(),
                               np.array([[0, 1, 1], [1, 1, 1]]))


if __name__ == "__main__":
    main()
