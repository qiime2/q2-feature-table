# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import qiime2
import pandas as pd
from pandas.util.testing import assert_series_equal

from q2_feature_table import filter_seqs


class FilterSeqsTests(unittest.TestCase):
    def setUp(self):
        self.seqs = pd.Series(['ACGT', 'GCTA', 'CCCC', 'TGTT'],
                              index=['O1', 'O2', 'O3', 'O4'])
        self.df_lite = pd.DataFrame(
            [['A'], ['C'], ['G'], ['T']], index=['O1', 'O2', 'O3', 'O4'],
            columns=['seq'])
        md_full = pd.DataFrame(
            [['foo', '1'], ['bar', '2'], ['baz', '3'], ['foo', '4']],
            index=['O1', 'O2', 'O3', 'O4'], columns=['stuff', 'some_numbers'])
        md_full.index.name = 'FeatureID'
        self.md_full = qiime2.Metadata(md_full)

    def filter_and_assertEqual(self, exp, md=None, exclude_ids=False,
                               where=None):
        if md is None:
            md = self.md_full
        obs = filter_seqs(self.seqs, md, exclude_ids=exclude_ids, where=where)
        assert_series_equal(exp, obs)

    def test_id_based_filtering(self):
        # filter none
        self.filter_and_assertEqual(self.seqs,
                                    md=qiime2.Metadata(self.df_lite))

        # filter one
        md = qiime2.Metadata(self.df_lite.drop(['O1']))
        exp = pd.Series(['GCTA', 'CCCC', 'TGTT'], index=['O2', 'O3', 'O4'])
        self.filter_and_assertEqual(exp, md=md)

        # filter all
        md = qiime2.Metadata(pd.DataFrame({}, index=['foo']))
        with self.assertRaisesRegex(ValueError, 'All.*filtered'):
            filter_seqs(self.seqs, md)

        # exclude none
        md = qiime2.Metadata(pd.DataFrame({}, index=['foo']))
        self.filter_and_assertEqual(self.seqs, md=md, exclude_ids=True)

        # exclude one
        md = qiime2.Metadata(self.df_lite.drop(['O1', 'O2', 'O3']))
        exp = pd.Series(['ACGT', 'GCTA', 'CCCC'], index=['O1', 'O2', 'O3'])
        self.filter_and_assertEqual(exp, md=md, exclude_ids=True)

        # exclude all
        md = qiime2.Metadata(self.df_lite)
        with self.assertRaisesRegex(ValueError, 'All.*filtered'):
            filter_seqs(self.seqs, md, exclude_ids=True)

    def test_id_based_filtering_with_extra_ids(self):
        md = qiime2.Metadata(pd.DataFrame([], index=['O1', 'O3', 'foo']))
        exp = pd.Series(['ACGT', 'CCCC'], index=['O1', 'O3'])
        self.filter_and_assertEqual(exp, md=md)

    def test_where_param(self):
        # filter none
        where = "stuff='foo' OR stuff='bar' OR stuff='baz'"
        self.filter_and_assertEqual(self.seqs, where=where)

        # filter one
        where = "stuff='foo' OR stuff='bar'"
        exp = pd.Series(['ACGT', 'GCTA', 'TGTT'], index=['O1', 'O2', 'O4'])
        self.filter_and_assertEqual(exp, where=where)

        # filter all
        where = "stuff='boo'"
        with self.assertRaisesRegex(ValueError, 'All.*filtered'):
            filter_seqs(self.seqs, self.md_full, where=where)

        # exclude none
        where = 'CAST(some_numbers AS INTEGER) < 0'
        self.filter_and_assertEqual(self.seqs, exclude_ids=True, where=where)

        # exclude one
        where = 'CAST(some_numbers AS INTEGER) > 3'
        exp = pd.Series(['ACGT', 'GCTA', 'CCCC'], index=['O1', 'O2', 'O3'])
        self.filter_and_assertEqual(exp, exclude_ids=True, where=where)

        # exclude all
        where = 'CAST(some_numbers AS INTEGER) BETWEEN 0 AND 5'
        with self.assertRaisesRegex(ValueError, 'All.*filtered'):
            filter_seqs(self.seqs, self.md_full, where=where, exclude_ids=True)


if __name__ == "__main__":
    unittest.main()
