# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def merge_tables(table1: biom.Table, table2: biom.Table) -> biom.Table:
    table1_sids = set(table1.ids(axis='sample'))
    table2_sids = set(table2.ids(axis='sample'))
    if len(table1_sids & table2_sids) > 0:
        raise ValueError('Some samples are present in both tables: %s' %
                         ', '.join(table1_sids & table2_sids))
    return table1.merge(table2)
