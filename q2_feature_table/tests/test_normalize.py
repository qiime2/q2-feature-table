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

from q2_feature_table import rarefy


class RarefyTests(TestCase):

    def test_rarefy(self):
        t = Table(np.array([[0, 1, 3], [1, 1, 2]]),
                  ['O1', 'O2'],
                  ['S1', 'S2', 'S3'])
        a = rarefy(t, 2)
        self.assertEqual(a.shape, (2, 2))
        self.assertEqual(set(a.ids(axis='sample')), set(['S2', 'S3']))
        self.assertEqual(set(a.ids(axis='observation')), set(['O1', 'O2']))
        npt.assert_array_equal(a.sum(axis='sample'), np.array([2., 2.]))


if __name__ == "__main__":
    main()
