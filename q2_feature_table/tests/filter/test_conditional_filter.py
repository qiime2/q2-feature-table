# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import biom
import numpy as np
import numpy.testing as npt

from q2_feature_table import filter_features_conditionally


class TestConditional(unittest.TestCase):
    def test_conditional(self):
        table = biom.Table(
            data=np.array([[0,   0,  10,   0,   0],
                           [250, 250, 140,  90, 150],
                           [250,  25, 100, 200, 100],
                           [0, 225, 250, 210, 250]]),
            sample_ids=['A', 'B', 'C', 'D', 'E'],
            observation_ids=['bat', 'cat', 'rat', 'a-tat-tat']
            )
        known = biom.Table(
            data=np.array([[0, 225, 250, 210, 250]]),
            sample_ids=['A', 'B', 'C', 'D', 'E'],
            observation_ids=['a-tat-tat']
            )
        test_ = filter_features_conditionally(table,
                                              0.8, 0.4)
        npt.assert_array_equal(known.matrix_data.toarray(),
                               test_.matrix_data.toarray())
        npt.assert_array_equal(known.ids(axis='sample'),
                               test_.ids(axis='sample'))
        npt.assert_array_equal(known.ids(axis='observation'),
                               test_.ids(axis='observation'))


if __name__ == "__main__":
    unittest.main()
