# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def rarefy(table, even_sampling_depth):
    return table.subsample(even_sampling_depth, axis='sample', by_id=False)
