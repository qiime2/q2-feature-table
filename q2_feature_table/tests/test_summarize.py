# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
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
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact, Metadata
import csv

from q2_feature_table import (
        tabulate_seqs, summarize,
        tabulate_feature_frequencies, tabulate_sample_frequencies)
from q2_feature_table._summarize._visualizer import _compute_descriptive_stats
from q2_feature_table._summarize._visualizer import _frequencies
from q2_feature_table._summarize._vega_spec import vega_spec


class TabulateSeqsTests(TestCase):

    def test_basic(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b) for a, b in (
            ('ACGT', {'id': 'seq1'}),
            ('AAAA', {'id': 'seq2'})))

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs)

            expected_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(expected_fp))
            with open(expected_fp) as fh:
                file_text = fh.read()
                self.assertTrue('ACGT</a>' in file_text)
                self.assertTrue('<td>4</td>' in file_text)
                self.assertTrue('<td>seq2</td>' in file_text)

    def test_descriptive_stats(self):
        seq_lengths = [2, 2, 5, 6, 10]
        exp_stats = {
            'mean': 5.0, 'min': 2,
            'seven_num_summ_values': [2.0, 2.0, 2.0, 5.0, 6.0, 8.56, 9.68],
            'max': 10, 'count': 5, 'range': 8}
        rendered_stats = _compute_descriptive_stats(seq_lengths)
        self.assertEqual(exp_stats['count'], rendered_stats['count'])
        self.assertEqual(exp_stats['min'], rendered_stats['min'])
        self.assertEqual(exp_stats['max'], rendered_stats['max'])
        self.assertEqual(exp_stats['range'], rendered_stats['range'])
        self.assertAlmostEqual(exp_stats['mean'], rendered_stats['mean'])
        for expected, rendered in zip(exp_stats['seven_num_summ_values'],
                                      rendered_stats['seven_num_summ_values']):
            self.assertAlmostEqual(expected, rendered)

    def test_lengths_identical(self):
        seq_lengths = [5, 5, 5, 5, 5]
        exp_stats = {
            'mean': 5.0, 'min': 5,
            'seven_num_summ_values': [5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0],
            'max': 5, 'count': 5, 'range': 0}
        rendered_stats = _compute_descriptive_stats(seq_lengths)
        self.assertEqual(exp_stats['count'], rendered_stats['count'])
        self.assertEqual(exp_stats['min'], rendered_stats['min'])
        self.assertEqual(exp_stats['max'], rendered_stats['max'])
        self.assertEqual(exp_stats['range'], rendered_stats['range'])
        self.assertAlmostEqual(exp_stats['mean'], rendered_stats['mean'])
        for expected, rendered in zip(exp_stats['seven_num_summ_values'],
                                      rendered_stats['seven_num_summ_values']):
            self.assertAlmostEqual(expected, rendered)

    def test_no_sequences(self):
        seq_lengths = []
        with self.assertRaisesRegex(ValueError, 'No values provided.'):
            _compute_descriptive_stats(seq_lengths)

    def test_descriptive_stats_integration(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b)for a, b in (
            ('A', {'id': 'seq01'}),
            ('AA', {'id': 'seq02'}),
            ('AAA', {'id': 'seq03'}),
            ('AAAA', {'id': 'seq04'}),
            ('AAAA', {'id': 'seq05'}),
            ('AAA', {'id': 'seq06'}),
            ('AA', {'id': 'seq07'}),
            ('AAAAAAAAAA', {'id': 'seq08'})))

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs)

            expected_fp = os.path.join(output_dir, 'index.html')

        # all expected values are unique. If they all render in index.html, our
        # function likely worked as expected.
            with open(expected_fp) as fh:
                file_text = fh.read()
                self.assertTrue('<td>8</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>10</td>' in file_text)
                self.assertTrue('<td>3.62</td>' in file_text)
                self.assertTrue('<td>9</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>2</td>' in file_text)
                self.assertTrue('<td>3</td>' in file_text)
                self.assertTrue('<td>4</td>' in file_text)
                self.assertTrue('<td>6</td>' in file_text)
                self.assertTrue('<td>9</td>' in file_text)

    def test_tsv_builder(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b)for a, b in (
            ('A', {'id': 'seq01'}),
            ('AA', {'id': 'seq02'}),
            ('AAA', {'id': 'seq03'}),
            ('AAAA', {'id': 'seq04'}),
            ('AAAA', {'id': 'seq05'}),
            ('AAA', {'id': 'seq06'}),
            ('AA', {'id': 'seq07'}),
            ('AAAAAAAAAA', {'id': 'seq08'})))

        # Do the files exist?
        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs)

            expected_stats_fp = os.path.join(
                output_dir, 'descriptive_stats.tsv')
            expected_summary_fp = os.path.join(
                output_dir, 'seven_number_summary.tsv')
            self.assertTrue(os.path.exists(expected_stats_fp))
            self.assertTrue(os.path.exists(expected_summary_fp))

            # Was data written to the files?
            with open(expected_stats_fp) as stats_tsv:
                tsv_reader = csv.reader(stats_tsv, dialect="excel-tab")
                tsv_text = []
                for row in tsv_reader:
                    tsv_text.append(row)
            self.assertEqual(['Statistic', 'Value'], tsv_text[0])
            self.assertEqual(['count', '8'], tsv_text[1])

            with open(expected_summary_fp) as summ_tsv:
                tsv_reader = csv.reader(summ_tsv, dialect="excel-tab")
                tsv_text = []
                for row in tsv_reader:
                    tsv_text.append(row)
            self.assertEqual(['Quantile', 'Value'], tsv_text[0])
            self.assertEqual(['0.02', '1.14'], tsv_text[1])

            # Does link html generate correctly?
            expected_index_fp = os.path.join(output_dir, 'index.html')
            with open(expected_index_fp) as fh:
                self.assertTrue('href="descriptive_stats.tsv"' in fh.read())

            with open(expected_index_fp) as fh:
                self.assertTrue(
                    'href="seven_number_summary.tsv"' in fh.read())

    def test_optional_inputs(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b)for a, b in (
            ('A', {'id': 'seq01'}),
            ('AA', {'id': 'seq02'}),
            ('AAA', {'id': 'seq03'}),
            ('AAAA', {'id': 'seq04'}),
            ('AAAA', {'id': 'seq05'}),
            ('AAA', {'id': 'seq06'}),
            ('AA', {'id': 'seq07'}),
            ('AAAAAAAAAA', {'id': 'seq08'})))

        metadata = pd.DataFrame(index=['seq01', 'seq02',
                                       'seq03', 'seq04',
                                       'seq05', 'seq06',
                                       'seq07', 'seq08'],
                                columns=['att1', 'att2'],
                                data=[['00', '01'], ['10', '11'],
                                      ['03', '04'], ['12', '13'],
                                      ['05', '06'], ['14', '15'],
                                      ['07', '08'], ['16', '17']])
        metadata.index.name = 'feature id'

        taxonomy = pd.DataFrame([('a;b;c;d', '1.0'), ('a;b;c;f', '0.7'),
                                 ('a;b;h;d', '0.3'), ('a;b;d;f', '0.7'),
                                 ('a;b;e;d', '0.4'), ('a;b;c;f', '0.6'),
                                 ('a;b;t;d', '1.0'), ('a;b;d;f', '0.5')],
                                index=['seq01', 'seq02', 'seq03', 'seq04',
                                       'seq05', 'seq06', 'seq07', 'seq08'],
                                columns=['Taxon', 'Confidence'])

        metadata = qiime2.Metadata(metadata)
        taxonomy = {"Taxon Name": taxonomy}

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs, metadata=metadata,
                          taxonomy=taxonomy)

            expected_fp = os.path.join(output_dir, 'index.html')
            with open(expected_fp) as fh:
                file_text = fh.read()
                self.assertTrue('<td>8</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>10</td>' in file_text)
                self.assertTrue('<td>3.62</td>' in file_text)
                self.assertTrue('<td>9</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>2</td>' in file_text)
                self.assertTrue('<td>3</td>' in file_text)
                self.assertTrue('<td>4</td>' in file_text)
                self.assertTrue('<td>6</td>' in file_text)
                self.assertTrue('<td>a;b;e;d</td>' in file_text)
                self.assertTrue('<td>a;b;c;d</td>' in file_text)
                self.assertTrue('<td>08</td>' in file_text)
                self.assertTrue('<td>17</td>' in file_text)

    def test_optional_input_union(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b)for a, b in (
            ('A', {'id': 'seq01'}),
            ('AA', {'id': 'seq02'}),
            ('AAA', {'id': 'seq03'}),
            ('AAAA', {'id': 'seq04'}),
            ('AAAA', {'id': 'seq05'}),
            ('AAA', {'id': 'seq06'}),
            ('AA', {'id': 'seq07'}),
            ('AAAAAAAAAA', {'id': 'seq08'})))

        metadata = pd.DataFrame(index=['seq01', 'seq02',
                                       'seq15', 'seq04',
                                       'seq05', 'seq90',
                                       'seq48', 'seq08'],
                                columns=['att1', 'att2'],
                                data=[['00', '01'], ['10', '11'],
                                      ['03', '04'], ['12', '13'],
                                      ['05', '06'], ['14', '15'],
                                      ['07', '08'], ['16', '17']])
        metadata.index.name = 'feature id'

        taxonomy = pd.DataFrame([('a;b;c;d', '1.0'), ('a;b;c;f', '0.7'),
                                 ('a;b;h;d', '0.3'), ('a;b;d;f', '0.7'),
                                 ('a;b;e;d', '0.4'), ('a;b;c;f', '0.6'),
                                 ('a;b;t;d', '1.0'), ('a;b;d;f', '0.5')],
                                index=['seq17', 'seq02', 'seq03', 'seq48',
                                       'seq05', 'seq19', 'seq07', 'seq08'],
                                columns=['Taxon', 'Confidence'])
        taxonomy2 = pd.DataFrame([('a;b;c;d', '1.0'), ('a;b;c;f', '0.7'),
                                 ('a;b;h;d', '0.3'), ('a;b;d;f', '0.7'),
                                 ('a;b;e;d', '0.4'), ('a;b;c;f', '0.6'),
                                 ('a;b;t;d', '1.0'), ('a;b;d;f', '0.5')],
                                 index=['seq99', 'seq02', 'seq44', 'seq48',
                                        'seq05', 'seq67', 'seq07', 'seq08'],
                                 columns=['Taxon', 'Confidence'])

        metadata = qiime2.Metadata(metadata)
        taxonomy = {"foo": taxonomy,
                    "bar": taxonomy2}

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs, metadata=metadata,
                          taxonomy=taxonomy, merge_method="union")
            expected_fp = os.path.join(output_dir, 'index.html')
            with open(expected_fp) as fh:
                file_text = fh.read()
                self.assertTrue('<td>8</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('<td>10</td>' in file_text)
                self.assertTrue('<td>3.62</td>' in file_text)
                self.assertTrue('<td>9</td>' in file_text)
                self.assertTrue('<td>1</td>' in file_text)
                self.assertTrue('seq17' in file_text)
                self.assertTrue('<td>2</td>' in file_text)
                self.assertTrue('<td>3</td>' in file_text)
                self.assertTrue('<td>4</td>' in file_text)
                self.assertTrue('<td>6</td>' in file_text)
                self.assertTrue('<td>a;b;e;d</td>' in file_text)
                self.assertTrue('<td>a;b;c;d</td>' in file_text)
                self.assertTrue('<td>08</td>' in file_text)
                self.assertTrue('foo' in file_text)
                self.assertTrue('bar' in file_text)

    def test_optional_inputs_intersect(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b)for a, b in (
            ('A', {'id': 'seq01'}),
            ('AA', {'id': 'seq02'}),
            ('AAA', {'id': 'seq03'}),
            ('AAAA', {'id': 'seq04'}),
            ('AAAA', {'id': 'seq05'}),
            ('AAA', {'id': 'seq06'}),
            ('AA', {'id': 'seq07'}),
            ('AAAAAAAAAA', {'id': 'seq08'})))

        metadata = pd.DataFrame(index=['seq01', 'seq02',
                                       'seq15', 'seq04',
                                       'seq05', 'seq90',
                                       'seq48', 'seq08'],
                                columns=['att1', 'att2'],
                                data=[['00', '01'], ['10', '11'],
                                      ['03', '04'], ['12', '13'],
                                      ['05', '06'], ['14', '15'],
                                      ['07', '08'], ['16', '17']])
        metadata.index.name = 'feature id'

        taxonomy = pd.DataFrame([('a;b;c;d', '1.0'), ('a;b;c;f', '0.7'),
                                 ('a;b;h;d', '0.3'), ('a;b;d;f', '0.7'),
                                 ('a;b;e;d', '0.4'), ('a;b;c;f', '0.6'),
                                 ('a;b;t;d', '1.0'), ('a;b;d;f', '0.5')],
                                index=['seq17', 'seq02', 'seq03', 'seq48',
                                       'seq05', 'seq19', 'seq07', 'seq08'],
                                columns=['Taxon', 'Confidence'])

        metadata = qiime2.Metadata(metadata)
        taxonomy = {"foo": taxonomy}

        with tempfile.TemporaryDirectory() as output_dir:
            tabulate_seqs(output_dir, seqs, metadata=metadata,
                          taxonomy=taxonomy, merge_method="intersect")
            expected_fp = os.path.join(output_dir, 'index.html')
            with open(expected_fp) as fh:
                file_text = fh.read()
                self.assertFalse('<td>a;b;c;d</td>' in file_text)
                self.assertFalse('<td>14</td>' in file_text)
                self.assertTrue('<td>10</td>' in file_text)

    def test_optional_input_strict_fail(self):
        seqs = DNAIterator(skbio.DNA(a, metadata=b)for a, b in (
            ('A', {'id': 'seq01'}),
            ('AA', {'id': 'seq02'}),
            ('AAA', {'id': 'seq03'}),
            ('AAAA', {'id': 'seq04'}),
            ('AAAA', {'id': 'seq05'}),
            ('AAA', {'id': 'seq06'}),
            ('AA', {'id': 'seq07'}),
            ('AAAAAAAAAA', {'id': 'seq08'})))

        metadata = pd.DataFrame(index=['seq01', 'seq02',
                                       'seq15', 'seq04',
                                       'seq05', 'seq90',
                                       'seq48', 'seq08'],
                                columns=['att1', 'att2'],
                                data=[['00', '01'], ['10', '11'],
                                      ['03', '04'], ['12', '13'],
                                      ['05', '06'], ['14', '15'],
                                      ['07', '08'], ['16', '17']])
        metadata.index.name = 'feature id'

        taxonomy = pd.DataFrame([('a;b;c;d', '1.0'), ('a;b;c;f', '0.7'),
                                 ('a;b;h;d', '0.3'), ('a;b;d;f', '0.7'),
                                 ('a;b;e;d', '0.4'), ('a;b;c;f', '0.6'),
                                 ('a;b;t;d', '1.0'), ('a;b;d;f', '0.5')],
                                index=['seq17', 'seq02', 'seq03', 'seq48',
                                       'seq05', 'seq19', 'seq07', 'seq08'],
                                columns=['Taxon', 'Confidence'])
        metadata = qiime2.Metadata(metadata)
        taxonomy = {"Taxon Name": taxonomy}

        with self.assertRaisesRegex(Exception, "Merge method is strict"):
            with tempfile.TemporaryDirectory() as output_dir:
                tabulate_seqs(output_dir, seqs, metadata=metadata,
                              taxonomy=taxonomy)
                # Did not error out, this is a problem


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

    def test_frequency_ranges_are_zero(self):
        table = biom.Table(np.array([[25, 25, 25], [25, 25, 25]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])

        with tempfile.TemporaryDirectory() as output_dir:
            summarize(output_dir, table)

            index_fp = os.path.join(output_dir, 'index.html')
            self.assertTrue(os.path.exists(index_fp))

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

    def test_vega_spec_data(self):
        # test if metadata is converted correctly to vega compatible JSON
        df = pd.DataFrame({'Subject': ['subject-1', 'subject-1', 'subject-2'],
                           'SampleType': ['gut', 'tongue', 'gut']},
                          index=pd.Index(['S1', 'S2', 'S3'], name='id'))
        metadata = qiime2.Metadata(df)
        table = biom.Table(np.array([[0, 1, 3], [1, 1, 2]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        sample_frequencies = _frequencies(table, axis='sample')
        spec = vega_spec(metadata, sample_frequencies)

        self.assertTrue([{'id': 'S1', 'metadata': {'Subject': 'subject-1',
                          'SampleType': 'gut'}, 'frequency': 1.0},
                         {'id': 'S2', 'metadata': {'Subject': 'subject-1',
                          'SampleType': 'tongue'}, 'frequency': 2.0},
                         {'id': 'S3', 'metadata': {'Subject': 'subject-2',
                          'SampleType': 'gut'}, 'frequency': 5.0}],
                        spec['data'][0]['values'])

    def test_vega_spec_nandling(self):
        df = pd.DataFrame({'a': [0.5, float('nan')]})
        df.index = df.index.map(str)
        df.index.name = 'id'
        md = qiime2.Metadata(df)
        sample_freqs = pd.Series([10, 50])
        sample_freqs.index = sample_freqs.index.map(str)

        spec = vega_spec(md, sample_freqs)
        exp = [{'frequency': 10, 'id': '0', 'metadata': {'a': 0.5}},
               {'frequency': 50, 'id': '1', 'metadata': {'a': None}}]

        self.assertEqual(spec['data'][0]['values'], exp)


class TabulateSampleFrequencyTests(TestCase):

    def test_basic_case(self):
        table = biom.Table(np.array([[0, 25, 25], [25, 25, 25]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        obs = tabulate_sample_frequencies(table).to_dataframe()

        exp = pd.DataFrame({'Frequency': ['25.0', '50.0', '50.0'],
                            'No. of Associated Features':
                            ['1', '2', '2']},
                           index=['S1', 'S2', 'S3'])
        exp.index.name = 'Sample ID'
        pd.testing.assert_frame_equal(exp, obs)


class TabulateFeatureFrequencyTests(TestCase):

    def test_basic_case(self):
        table = biom.Table(np.array([[25, 25, 0], [25, 25, 25]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        obs = tabulate_feature_frequencies(table).to_dataframe()

        exp = pd.DataFrame({'Frequency': ['50.0', '75.0'],
                            'No. of Samples Observed In':
                            ['2', '3']},
                           index=['O1', 'O2'])
        exp.index.name = 'Feature ID'
        pd.testing.assert_frame_equal(exp, obs)


class SummarizePlusTests(TestPluginBase):

    package = 'q2_feature_table'

    def setUp(self):
        super().setUp()
        self.summarize_plus = self.plugin.pipelines['summarize_plus']

    def test_basic(self):
        table = biom.Table(np.array([[25, 0, 25], [25, 25, 25]]),
                           ['O1', 'O2'],
                           ['S1', 'S2', 'S3'])
        table = Artifact.import_data('FeatureTable[Frequency]', table)
        results = self.summarize_plus(table)

        self.assertEqual(len(results), 3)
        self.assertEqual(repr(results.feature_frequencies.type),
                         'ImmutableMetadata')
        self.assertEqual(repr(results.sample_frequencies.type),
                         'ImmutableMetadata')
        self.assertEqual(repr(results.summary.type),
                         'Visualization')

        exp_feature = pd.DataFrame({'Frequency': ['50.0', '75.0'],
                                   'No. of Samples Observed In':
                                    ['2', '3']},
                                   index=['O1', 'O2'])
        exp_feature.index.name = "Feature ID"
        obs_feature = results[0].view(Metadata).to_dataframe()
        pd.testing.assert_frame_equal(exp_feature, obs_feature)

        exp_sample = pd.DataFrame({'Frequency': ['50.0', '25.0', '50.0'],
                                  'No. of Associated Features':
                                   ['2', '1', '2']},
                                  index=['S1', 'S2', 'S3'])
        exp_sample.index.name = "Sample ID"
        obs_sample = results[1].view(Metadata).to_dataframe()
        pd.testing.assert_frame_equal(exp_sample, obs_sample)


if __name__ == "__main__":
    main()
