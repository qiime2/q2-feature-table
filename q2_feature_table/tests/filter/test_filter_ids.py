# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import pandas as pd
import qiime2
from biom.table import Table
from rachis.plugin.testing import TestPluginBase

from q2_feature_table import filter_ids


class FilterIDsTests(TestPluginBase):
    package = "q2_feature_table.tests"

    def setUp(self):
        super().setUp()
        self.table = Table(np.array([[1, 0, 0], [0, 2, 0]]),
                           ["O1", "O2"], ["S1", "S2", "S3"])

    def test_filter_samples_by_where(self):
        metadata = qiime2.Metadata(pd.DataFrame(
            {"Group": ["keep", "drop"]},
            index=pd.Index(["S1", "S2"], name="id")))

        actual = filter_ids(self.table.copy(), axis="sample",
                            metadata=metadata, where="\"Group\"='keep'")
        expected = Table(np.array([[1], [0]]), ["O1", "O2"], ["S1"])

        self.assertEqual(actual, expected)

    def test_filter_samples_by_metadata(self):
        metadata = qiime2.Metadata(pd.DataFrame(
            {"Group": ["keep"]}, index=pd.Index(["S1"], name="id")))

        actual = filter_ids(self.table.copy(), axis="sample",
                            metadata=metadata)
        expected = Table(np.array([[1], [0]]), ["O1", "O2"], ["S1"])

        self.assertEqual(actual, expected)

    def test_filter_features_by_metadata(self):
        metadata = qiime2.Metadata(pd.DataFrame(
            {"Group": ["keep"]}, index=pd.Index(["O2"], name="id")))

        actual = filter_ids(self.table.copy(), axis="feature",
                            metadata=metadata)
        expected = Table(np.array([[0, 2, 0]]), ["O2"],
                         ["S1", "S2", "S3"])

        self.assertEqual(actual, expected)

    def test_filter_features_by_where(self):
        metadata = qiime2.Metadata(pd.DataFrame(
            {"Group": ["drop", "keep"]},
            index=pd.Index(["O1", "O2"], name="id")))

        actual = filter_ids(self.table.copy(), axis="feature",
                            metadata=metadata,
                            where="\"Group\"='keep'")
        expected = Table(np.array([[0, 2, 0]]), ["O2"],
                         ["S1", "S2", "S3"])

        self.assertEqual(actual, expected)

    def test_exclude_ids(self):
        actual = filter_ids(self.table.copy(), axis="sample", ids=["S1"],
                            exclude_ids=True)
        expected = Table(np.array([[0, 0], [2, 0]]),
                         ["O1", "O2"], ["S2", "S3"])

        self.assertEqual(actual, expected)

    def test_filter_ids_retains_signed_vectors(self):
        table = Table(np.array([[2, -2, 0], [0, 0, -3], [0, 0, 0]]),
                      ["O1", "O2", "O3"], ["S1", "S2", "S3"])

        actual = filter_ids(table, axis="sample", ids=["S2"])
        expected = Table(np.array([[-2], [0], [0]]),
                         ["O1", "O2", "O3"], ["S2"])

        self.assertEqual(actual, expected)

    def test_empty_result(self):
        metadata = qiime2.Metadata(pd.DataFrame(
            {"Group": ["keep"]},
            index=pd.Index(["not-in-table"], name="id")))

        with self.assertRaisesRegex(ValueError, "table is empty"):
            filter_ids(self.table.copy(), axis="sample", metadata=metadata,
                       where="\"Group\"='keep'")

        actual = filter_ids(self.table.copy(), axis="sample",
                            metadata=metadata, where="\"Group\"='keep'",
                            allow_empty_table=True)
        self.assertTrue(actual.is_empty())

    def test_filter_ids_parameter(self):
        actual = filter_ids(self.table.copy(), axis="sample", ids=["S2"])
        expected = Table(np.array([[0], [2]]), ["O1", "O2"], ["S2"])

        self.assertEqual(actual, expected)

    def test_filter_features_by_ids(self):
        actual = filter_ids(self.table.copy(), axis="feature", ids=["O1"])
        expected = Table(np.array([[1, 0, 0]]), ["O1"],
                         ["S1", "S2", "S3"])

        self.assertEqual(actual, expected)

    def test_ids_and_metadata_are_mutually_exclusive(self):
        metadata = qiime2.Metadata(pd.DataFrame(
            {"Group": ["keep"]}, index=pd.Index(["S1"], name="id")))

        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            filter_ids(self.table.copy(), axis="sample", ids=["S1"],
                       metadata=metadata)

    def test_ids_must_be_in_table(self):
        with self.assertRaisesRegex(ValueError, "not in the table"):
            filter_ids(self.table.copy(), axis="sample",
                       ids=["S1", "not-in-table"])

    def test_where_requires_metadata(self):
        with self.assertRaisesRegex(ValueError, "Metadata must be provided"):
            filter_ids(self.table.copy(), axis="sample", where="id='S1'")

    def test_no_filtering_requested(self):
        with self.assertRaisesRegex(ValueError, "No filtering was requested"):
            filter_ids(self.table.copy(), axis="sample")


if __name__ == "__main__":
    unittest.main()
