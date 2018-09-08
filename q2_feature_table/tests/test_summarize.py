# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from unittest import TestCase, main
import tempfile

import skbio
import biom
import pandas as pd
import numpy as np
import qiime2
from q2_types.feature_data import DNAIterator

from q2_feature_table import tabulate_seqs, summarize
from q2_feature_table._summarize._visualizer import _compute_descriptive_stats


class TabulateSeqsTests(TestCase):

    def test_basic(self):
        seqs = DNAIterator(
            (s for s in (skbio.DNA('ACGT', metadata={'id': 'seq1'}),
                         skbio.DNA('AAAA', metadata={'id': 'seq2'}))))

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs)

            expected_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(expected_fp))
            self.assertTrue('ACGT</a>' in open(expected_fp).read())
            self.assertTrue('<td>4</td>' in open(expected_fp).read())
            self.assertTrue('<td>seq2</td>' in open(expected_fp).read())

    def test_descriptive_stats(self):
        seq_lengths = [2, 2, 5, 6, 10]
        exp_stats = {
            'mean': 5.0, 'min': 2,
            'seven_num_summ': [2.0, 2.0, 2.0, 5.0, 6.0, 8.56, 9.68],
            'max': 10, 'count': 5}
        rendered_stats = _compute_descriptive_stats(seq_lengths)
        self.assertEqual(exp_stats['count'], rendered_stats['count'])
        self.assertEqual(exp_stats['min'], rendered_stats['min'])
        self.assertEqual(exp_stats['max'], rendered_stats['max'])
        self.assertAlmostEqual(exp_stats['mean'], rendered_stats['mean'])
        for expected, rendered in zip(exp_stats['seven_num_summ'],
                                      rendered_stats['seven_num_summ']):
            self.assertAlmostEqual(expected, rendered)

    def test_lengths_identical(self):
        seq_lengths = [5, 5, 5, 5, 5]
        exp_stats = {
            'mean': 5.0, 'min': 5,
            'seven_num_summ': [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            'max': 5, 'count': 5}
        rendered_stats = _compute_descriptive_stats(seq_lengths)
        self.assertEqual(exp_stats['count'], rendered_stats['count'])
        self.assertEqual(exp_stats['min'], rendered_stats['min'])
        self.assertEqual(exp_stats['max'], rendered_stats['max'])
        self.assertAlmostEqual(exp_stats['mean'], rendered_stats['mean'])
        for expected, rendered in zip(exp_stats['seven_num_summ'],
                                      rendered_stats['seven_num_summ']):
            self.assertAlmostEqual(expected, rendered)

    def test_no_sequences(self):
        seq_lengths = []
        with self.assertRaisesRegex(ValueError, 'No sequences provided.'):
            _compute_descriptive_stats(seq_lengths)

    def test_descriptive_stats_integration(self):
        seqs = DNAIterator(
            (s for s in (skbio.DNA('A', metadata={'id': 'seq01'}),
                         skbio.DNA('AA', metadata={'id': 'seq02'}),
                         skbio.DNA('AAA', metadata={'id': 'seq03'}),
                         skbio.DNA('AAAA', metadata={'id': 'seq04'}),
                         skbio.DNA('AAAA', metadata={'id': 'seq05'}),
                         skbio.DNA('AAA', metadata={'id': 'seq06'}),
                         skbio.DNA('AA', metadata={'id': 'seq07'}),
                         skbio.DNA('AAAAAAAAAA', metadata={'id': 'seq08'}))))

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs)

            expected_fp = os.path.join(output_dir, 'index.html')

# all expected summary values are unique.
# if they all render in index.html, our function likely worked as expected
            self.assertTrue('<td>8</td>' in open(expected_fp).read())
            self.assertTrue('<td>1</td>' in open(expected_fp).read())
            self.assertTrue('<td>10</td>' in open(expected_fp).read())
            self.assertTrue('<td>3.62</td>' in open(expected_fp).read())
            self.assertTrue('<td>1.14</td>' in open(expected_fp).read())
            self.assertTrue('<td>1.63</td>' in open(expected_fp).read())
            self.assertTrue('<td>2</td>' in open(expected_fp).read())
            self.assertTrue('<td>3</td>' in open(expected_fp).read())
            self.assertTrue('<td>4</td>' in open(expected_fp).read())
            self.assertTrue('<td>6.22</td>' in open(expected_fp).read())
            self.assertTrue('<td>9.16</td>' in open(expected_fp).read())


class SummarizeTests(TestCase):

    def test_basic(self):
        table = biom.Table(np.array([[0, 1, 3],
                                     [1, 1, 2],
                                     [400, 450, 500],
                                     [1000, 10000, 100000],
                                     [52, 42, 99]]),
                           ['O1', 'O2', '03', '04', 'O5'],
                           ['S1', 'S2', 'S3'])

        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))

            feature_freq_fp = os.path.join(output_dir,
                                           'feature-frequency-detail.csv')
            self.assertTrue(os.path.exists(feature_freq_fp))
            self.assertTrue('O1,4' in open(feature_freq_fp).read())

            sample_freq_fp = os.path.join(output_dir,
                                          'sample-frequency-detail.csv')
            self.assertTrue(os.path.exists(sample_freq_fp))
            self.assertTrue('S1,1453' in open(sample_freq_fp).read())

            interactive_sample_detail_fp = \
                os.path.join(output_dir, 'data.jsonp')
            self.assertTrue(os.path.exists(interactive_sample_detail_fp))

    def test_frequency_ranges_are_zero(self):
        table = biom.Table(np.array([[25, 25, 25], [25, 25, 25]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])

        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))

            feature_freq_fp = os.path.join(output_dir,
                                           'feature-frequency-detail.csv')
            self.assertTrue(os.path.exists(feature_freq_fp))
            self.assertTrue('O1,75' in open(feature_freq_fp).read())

            sample_freq_fp = os.path.join(output_dir,
                                          'sample-frequency-detail.csv')
            self.assertTrue(os.path.exists(sample_freq_fp))
            self.assertTrue('S1,50' in open(sample_freq_fp).read())

            interactive_sample_detail_fp = \
                os.path.join(output_dir, 'data.jsonp')
            self.assertTrue(os.path.exists(interactive_sample_detail_fp))

    def test_one_sample(self):
        sample_frequencies_pdf_fn = 'sample-frequencies.pdf'
        # sample-frequencies.pdf should not be written when there is only
        # one sample...
        table = biom.Table(np.array([[0], [1]]),
                           ['O1', 'O2'],
                           ['S1'])
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)
            sample_frequencies_pdf_fp = \
                os.path.join(output_dir, sample_frequencies_pdf_fn)
            self.assertFalse(os.path.exists(sample_frequencies_pdf_fp))

        # but it should be written when there is more than one sample
        table = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)
            sample_frequencies_pdf_fp = \
                os.path.join(output_dir, sample_frequencies_pdf_fn)
            self.assertTrue(os.path.exists(sample_frequencies_pdf_fp))

    def test_one_feature(self):
        feature_frequencies_pdf_fn = 'feature-frequencies.pdf'
        # feature-frequencies.pdf should not be written when there is only
        # one feature...
        table = biom.Table(np.array([[0, 4]]),
                           ['O1'],
                           ['S1', 'S2'])
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)
            feature_frequencies_pdf_fp = \
                os.path.join(output_dir, feature_frequencies_pdf_fn)
            self.assertFalse(os.path.exists(feature_frequencies_pdf_fp))

        # but it should be written when there is more than one feature
        table = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)
            feature_frequencies_pdf_fp = \
                os.path.join(output_dir, feature_frequencies_pdf_fn)
            self.assertTrue(os.path.exists(feature_frequencies_pdf_fp))

    def test_w_sample_metadata(self):
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)
        table = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])

        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table, metadata)

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))

            feature_freq_fp = os.path.join(output_dir,
                                           'feature-frequency-detail.csv')
            self.assertTrue(os.path.exists(feature_freq_fp))
            self.assertTrue('O1,4' in open(feature_freq_fp).read())

            sample_freq_fp = os.path.join(output_dir,
                                          'sample-frequency-detail.csv')
            self.assertTrue(os.path.exists(sample_freq_fp))
            self.assertTrue('S1,1' in open(sample_freq_fp).read())

            interactive_sample_detail_fp = \
                os.path.join(output_dir, 'data.jsonp')
            self.assertTrue(os.path.exists(interactive_sample_detail_fp))
            self.assertTrue('SampleType' in
                            open(interactive_sample_detail_fp).read())


if __name__ == "__main__":
    main()
