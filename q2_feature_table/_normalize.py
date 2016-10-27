# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def rarefy(table: biom.Table, sampling_depth: int) -> biom.Table:
    return table.subsample(sampling_depth, axis='sample', by_id=False)
