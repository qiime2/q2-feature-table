# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import biom
import numpy as np
import pandas as pd
import qiime2 as q2
from q2_types.feature_table import (
    FeatureTable, Frequency, PresenceAbsence, RelativeFrequency
)
from rachis.plugin.testing import TestPluginBase

from q2_feature_table import (
    _multiply_tables_pa, _multiply_tables_relative, _multiply_tables
)


class MultiplyTablesTests(TestPluginBase):
    package = "q2_feature_table.tests"

    @staticmethod
    def _df_to_biom(df):
        return biom.Table(
            df.values.T, sample_ids=df.index, observation_ids=df.columns
        )

    def setUp(self):
        super().setUp()
        table1 = pd.DataFrame(
            {"m1": [1, 4], "m2": [2, 5], "m3": [3, 6]}, index=["s1", "s2"]
        )
        self.table1 = self._df_to_biom(table1)
        table1_pa = pd.DataFrame(
            {"m1": [0, 1], "m2": [0, 1], "m3": [0, 0]}, index=["s1", "s2"]
        )
        self.table1_pa = self._df_to_biom(table1_pa)
        table1_rel = pd.DataFrame(
            {"m1": [0.1667, 0.2667], "m2": [0.3333, 0.3333], "m3": [0.5, 0.4]},
            index=["s1", "s2"],
        )
        self.table1_rel = self._df_to_biom(table1_rel)

        table2 = pd.DataFrame(
            {"a1": [7, 9, 11], "a2": [8, 10, 12]}, index=["m1", "m2", "m3"]
        )
        self.table2 = self._df_to_biom(table2)
        table2_pa = pd.DataFrame(
            {"a1": [0, 1, 0], "a2": [1, 0, 1]}, index=["m1", "m2", "m3"]
        )
        self.table2_pa = self._df_to_biom(table2_pa)
        table2_rel = pd.DataFrame(
            {"a1": [0.4667, 0.4737, 0.4783], "a2": [0.5333, 0.5263, 0.5217]},
            index=["m1", "m2", "m3"],
        )
        self.table2_rel = self._df_to_biom(table2_rel)
        self.multiply = self.plugin.pipelines["multiply_tables"]

    def test_multiply_tables(self):
        obs = _multiply_tables(self.table1, self.table2)
        exp = pd.DataFrame(
            {"a1": [58, 139], "a2": [64, 154]},
            dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_pa(self):
        obs = _multiply_tables_pa(self.table1_pa, self.table2)
        exp = pd.DataFrame(
            {"a1": [0, 1], "a2": [0, 1]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_pa_both(self):
        obs = _multiply_tables_pa(self.table1_pa, self.table2_pa)
        exp = pd.DataFrame(
            {"a1": [0, 1], "a2": [0, 1]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_relative(self):
        obs = _multiply_tables_relative(self.table1_rel, self.table2)
        exp = pd.DataFrame(
            {"a1": [0.4754, 0.4744], "a2": [0.5246, 0.5256]},
            index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(
            obs.to_dataframe(dense=True).T, exp, atol=1e-4, check_exact=False
        )

    def test_multiply_tables_relative_both(self):
        obs = _multiply_tables_relative(self.table1_rel, self.table2_rel)
        exp = pd.DataFrame(
            {"a1": [0.4748, 0.4737], "a2": [0.5252, 0.5263]},
            index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(
            obs.to_dataframe(dense=True).T, exp, atol=1e-4, check_exact=False
        )

    def test_multiply_tables_partial_overlap(self):
        # Test when table1 has some observations not in table2 samples
        table1_partial = pd.DataFrame(
            {"m1": [1, 4], "m2": [2, 5], "m4": [7, 8]}, index=["s1", "s2"]
        )
        table1_partial = self._df_to_biom(table1_partial)

        with self.assertWarnsRegex(UserWarning, r"Removed 1 feature\(s\)"):
            obs = _multiply_tables(table1_partial, self.table2)

        # verify the result only includes overlapping features
        exp = pd.DataFrame(
            {"a1": [25, 73], "a2": [28, 82]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_pipeline_freq_freq(self):
        table1 = q2.Artifact.import_data(
            "FeatureTable[Frequency]", self.table1
        )
        table2 = q2.Artifact.import_data(
            "FeatureTable[Frequency]", self.table2
        )
        (obs,) = self.multiply(table1, table2)
        self.assertEqual(obs.type, FeatureTable[Frequency])

    def test_multiply_tables_pipeline_pa_freq(self):
        table1 = q2.Artifact.import_data(
            "FeatureTable[PresenceAbsence]", self.table1_pa
        )
        table2 = q2.Artifact.import_data(
            "FeatureTable[Frequency]", self.table2
        )
        (obs,) = self.multiply(table1, table2)
        self.assertEqual(obs.type, FeatureTable[PresenceAbsence])

    def test_multiply_tables_pipeline_pa_rel(self):
        table1 = q2.Artifact.import_data(
            "FeatureTable[PresenceAbsence]", self.table1_pa
        )
        table2 = q2.Artifact.import_data(
            "FeatureTable[RelativeFrequency]", self.table2_rel
        )
        (obs,) = self.multiply(table1, table2)
        self.assertEqual(obs.type, FeatureTable[PresenceAbsence])

    def test_multiply_tables_pipeline_pa_pa(self):
        table1 = q2.Artifact.import_data(
            "FeatureTable[PresenceAbsence]", self.table1_pa
        )
        table2 = q2.Artifact.import_data(
            "FeatureTable[PresenceAbsence]", self.table2_pa
        )
        (obs,) = self.multiply(table1, table2)
        self.assertEqual(obs.type, FeatureTable[PresenceAbsence])

    def test_multiply_tables_pipeline_freq_rel(self):
        table1 = q2.Artifact.import_data(
            "FeatureTable[Frequency]", self.table1
        )
        table2 = q2.Artifact.import_data(
            "FeatureTable[RelativeFrequency]", self.table2_rel
        )
        (obs,) = self.multiply(table1, table2)
        self.assertEqual(obs.type, FeatureTable[RelativeFrequency])

    def test_multiply_tables_pipeline_rel_rel(self):
        table1 = q2.Artifact.import_data(
            "FeatureTable[RelativeFrequency]", self.table1_rel
        )
        table2 = q2.Artifact.import_data(
            "FeatureTable[RelativeFrequency]", self.table2_rel
        )
        (obs,) = self.multiply(table1, table2)
        self.assertEqual(obs.type, FeatureTable[RelativeFrequency])

    def test_multiply_tables_no_overlap(self):
        # table1 with observations that don't overlap with table2 samples
        table1_no_overlap = pd.DataFrame(
            {"m4": [1, 4], "m5": [2, 5], "m6": [3, 6]}, index=["s1", "s2"]
        )
        table1_no_overlap = self._df_to_biom(table1_no_overlap)

        with self.assertRaisesRegex(
                ValueError, "No overlapping features found"
        ):
            _multiply_tables(table1_no_overlap, self.table2)

    def test_multiply_tables_empty_table1(self):
        table1_empty = biom.Table(
            np.array([]).reshape(0, 0), observation_ids=[], sample_ids=[]
        )

        with self.assertRaisesRegex(
                ValueError, "No overlapping features found"
        ):
            _multiply_tables(table1_empty, self.table2)

    def test_multiply_tables_empty_table2(self):
        # empty table2 (no observations)
        table2_empty = biom.Table(
            np.array([]).reshape(0, 0), observation_ids=[], sample_ids=[]
        )

        with self.assertRaisesRegex(
                ValueError, "No overlapping features found"
        ):
            _multiply_tables(self.table1, table2_empty)

    def test_multiply_tables_single_sample(self):
        table1_single = pd.DataFrame(
            {"m1": [1], "m2": [2], "m3": [3]}, index=["s1"]
        )
        table1_single = self._df_to_biom(table1_single)

        obs = _multiply_tables(table1_single, self.table2)
        exp = pd.DataFrame(
            {"a1": [58], "a2": [64]}, dtype="float", index=["s1"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_single_observation_table1(self):
        table1_single = pd.DataFrame({"m1": [1, 4]}, index=["s1", "s2"])
        table1_single = self._df_to_biom(table1_single)

        obs = _multiply_tables(table1_single, self.table2)
        exp = pd.DataFrame(
            {"a1": [7, 28], "a2": [8, 32]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_single_observation_table2(self):
        table2_single = pd.DataFrame({"a1": [7], "a2": [8]}, index=["m1"])
        table2_single = self._df_to_biom(table2_single)

        obs = _multiply_tables(self.table1, table2_single)
        exp = pd.DataFrame(
            {"a1": [7, 28], "a2": [8, 32]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_all_zeros_table1(self):
        table1_zeros = pd.DataFrame(
            {"m1": [0, 0], "m2": [0, 0], "m3": [0, 0]}, index=["s1", "s2"]
        )
        table1_zeros = self._df_to_biom(table1_zeros)

        obs = _multiply_tables(table1_zeros, self.table2)
        exp = pd.DataFrame(
            {"a1": [0, 0], "a2": [0, 0]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_all_zeros_table2(self):
        table2_zeros = pd.DataFrame(
            {"a1": [0, 0, 0], "a2": [0, 0, 0]}, index=["m1", "m2", "m3"]
        )
        table2_zeros = self._df_to_biom(table2_zeros)

        obs = _multiply_tables(self.table1, table2_zeros)
        exp = pd.DataFrame(
            {"a1": [0, 0], "a2": [0, 0]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_relative_zero_sum(self):
        # Test that normalization handles zero-sum samples gracefully
        table1_zeros = pd.DataFrame(
            {"m1": [0, 0], "m2": [0, 0], "m3": [0, 0]}, index=["s1", "s2"]
        )
        table1_zeros = self._df_to_biom(table1_zeros)

        obs = _multiply_tables_relative(table1_zeros, self.table2)
        # after norm, zero-sum samples should remain zero or become NaN
        self.assertEqual(obs.shape, (2, 2))  # 2 observations, 2 samples
        obs_df = obs.to_dataframe(dense=True).T
        obs_values = obs_df.to_numpy()
        self.assertFalse(np.isinf(obs_values).any())
        finite_mask = np.isfinite(obs_values)
        if finite_mask.any():
            self.assertTrue((obs_values[finite_mask] == 0).all())

    def test_multiply_tables_table2_extra_samples(self):
        # Test when table2 has extra samples not in table1 observations
        table2_extra = pd.DataFrame(
            {"a1": [7, 9, 11, 13], "a2": [8, 10, 12, 14]},
            index=["m1", "m2", "m3", "m4"],
        )
        table2_extra = self._df_to_biom(table2_extra)

        obs = _multiply_tables(self.table1, table2_extra)
        exp = pd.DataFrame(
            {"a1": [58, 139], "a2": [64, 154]},
            dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_table2_fewer_samples(self):
        # Test when table2 has fewer samples than table1 observations
        table2_fewer = pd.DataFrame(
            {"a1": [7, 9], "a2": [8, 10]}, index=["m1", "m2"]
        )
        table2_fewer = self._df_to_biom(table2_fewer)

        with self.assertWarnsRegex(UserWarning, r"Removed 1 feature\(s\)"):
            obs = _multiply_tables(self.table1, table2_fewer)

        exp = pd.DataFrame(
            {"a1": [25, 73], "a2": [28, 82]}, dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_pa_preserves_zeros(self):
        # Ensure PA conversion keeps zeros as zeros,
        # not just converts non-zeros to 1
        table1_mixed = pd.DataFrame(
            {"m1": [1, 0], "m2": [0, 5], "m3": [3, 0]}, index=["s1", "s2"]
        )
        table1_mixed = self._df_to_biom(table1_mixed)

        table2_mixed = pd.DataFrame(
            {"a1": [7, 0, 0], "a2": [0, 10, 0]}, index=["m1", "m2", "m3"]
        )
        table2_mixed = self._df_to_biom(table2_mixed)

        obs = _multiply_tables_pa(table1_mixed, table2_mixed)

        # s1: m1=1, m2=0, m3=3
        # s1*a1: 1*7 + 0*0 + 3*0 = 7 -> 1
        # s1*a2: 1*0 + 0*10 + 3*0 = 0 -> 0 (CRITICAL: zero stays zero)
        # s2: m1=0, m2=5, m3=0
        # s2*a1: 0*7 + 5*0 + 0*0 = 0 -> 0 (CRITICAL: zero stays zero)
        # s2*a2: 0*0 + 5*10 + 0*0 = 50 -> 1
        exp = pd.DataFrame(
            {"a1": [1.0, 0.0], "a2": [0.0, 1.0]},
            dtype="float", index=["s1", "s2"]
        )
        pd.testing.assert_frame_equal(obs.to_dataframe(dense=True).T, exp)

    def test_multiply_tables_very_sparse(self):
        # create 100x100 table with only 10 non-zero values (<1% density)
        size = 100
        sparse_data = np.zeros((size, size))
        # add a few non-zero values
        for i in range(10):
            sparse_data[i, i] = i + 1

        table1_sparse = biom.Table(
            sparse_data[:50, :].T,  # First 50 rows as observations
            observation_ids=[f"m{i}" for i in range(size)],
            sample_ids=[f"s{i}" for i in range(50)],
        )

        table2_sparse = biom.Table(
            sparse_data[:, :20].T,  # First 20 columns as observations
            observation_ids=[f"a{i}" for i in range(20)],
            sample_ids=[f"m{i}" for i in range(size)],
        )

        obs = _multiply_tables(table1_sparse, table2_sparse)

        self.assertEqual(obs.shape, (20, 50))  # 20 observations, 50 samples

        density = obs.matrix_data.nnz / (obs.shape[0] * obs.shape[1])
        self.assertLessEqual(density, 0.01)  # Less than 1% non-zero
