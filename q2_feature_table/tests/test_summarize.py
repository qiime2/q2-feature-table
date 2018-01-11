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
            self.assertTrue('<td>seq2</td>' in open(expected_fp).read())


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
                          index=['S1', 'S2', 'S3'])
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
