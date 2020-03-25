# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
import warnings

import biom
import qiime2
import pandas as pd
import numpy as np
import numpy.testing as npt

from q2_feature_table import _rename


class TestRename(unittest.TestCase):
    def setUp(self):
        self.old_ids = ['S1', 'S2', 'S3']
        self.name_map = pd.Series({'S1': 'S1_new', 
                                   'S2': 'S2_new', 
                                   'S4': 'S4_name'})
        self.known = {'S1': 'S1_new', 'S2': 'S2_new', 'S3': 'S3'}

    def test_generate_new_names_non_unique(self):
        name_map = pd.Series({'S1': 'S2_new', 'S2': 'S2_new'})
        with self.assertRaises(ValueError) as cm:
            _rename._generate_new_names(self.old_ids, 
                                        name_map, 
                                        strict=True, 
                                        verbose=False)
            self.assertEqual(
                str(cm.exception),
                ('All new ids must be unique.\n'
                 'Try qiime feature-table group if you want '
                 'to combine multipel samples in the same table.')
                )

    def test_generate_new_names_old_disjoint_strict(self):
        with self.assertRaises(ValueError) as cm:
            _rename._generate_new_names(self.old_ids, self.name_map, 
                                        strict=True, 
                                        verbose=False)
            self.assertEqual(
                str(cm.exception),
                ("There are ids in the table which do not have new names.\n"
                 "Either turn off strict mode or provide a remapping for all ids.\n"
                 "The following ids are not mapped:\n    S3")
                )

    def test_generate_new_names_verbose_warnings(self):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            new_names = \
                _rename._generate_new_names(self.old_ids, 
                                            self.name_map, 
                                            strict=False, 
                                            verbose=True)
        self.assertEqual(len(w), 2)
        self.assertTrue(isinstance(w[0].message, UserWarning))
        self.assertEqual(str(w[0].message), 
                         'There are ids in the original table which do not '
                         'have new names.\nThe following ids will not be '
                         'mapped:\n   S3')
        self.assertTrue(isinstance(w[1].message, UserWarning))
        self.assertEqual(str(w[1].message),
                         'There are ids supplied for renaming that are not in'
                         ' the table.\nThe following ids will not be mapped:'
                         '\n   S4'
                         )
        self.assertEqual(new_names.keys(), self.known.keys())
        for k, v in new_names.items():
            self.assertEqual(v, self.known[k])

    def test_generate_new_names_no_verbse(self):
        with warnings.catch_warnings(record=True) as w:
            # Cause all warnings to always be triggered.
            warnings.simplefilter("always")
            new_names = \
                _rename._generate_new_names(self.old_ids, 
                                            self.name_map, 
                                            strict=False, 
                                            verbose=False)
        self.assertEqual(len(w), 0)
        self.assertEqual(new_names.keys(), self.known.keys())
        for k, v in new_names.items():
            self.assertEqual(v, self.known[k])

    def test_rename_samples(self):
        table = biom.Table(np.array([[0, 1, 2], [3, 4, 5]]),
                            observation_ids=['01', '02'],
                            sample_ids=['S1', 'S2', 'S3'])
        meta = qiime2.Metadata(pd.DataFrame(
            data=np.array([['cat'], ['rat'], ['dog']]),
            index=pd.Index(['S1', 'S2', 'S3'], name='sample-id'),
            columns=['animal']
            ))
        updated = _rename.rename_samples(table, meta.get_column('animal'))
        
        npt.assert_array_equal(np.array(updated.ids(axis='sample')), 
                               np.array(['cat', 'rat', 'dog']))



if __name__ == "__main__":
    unittest.main()