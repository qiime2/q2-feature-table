# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import unittest
import tempfile

import skbio
from q2_types import DNAIterator

from q2_feature_table import tabulate_seqs


class TabulateSeqsTests(unittest.TestCase):

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
