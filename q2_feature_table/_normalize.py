# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def rarefy(table: biom.Table, counts_per_sample: int) -> biom.Table:
    return table.subsample(counts_per_sample, axis='sample', by_id=False)
