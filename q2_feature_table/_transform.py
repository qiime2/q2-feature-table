# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def presence_absence(table: biom.Table) -> biom.Table:
    """ Convert feature table to presence/absence data
    """
    return table.pa(inplace=False)


def relative_frequency(table: biom.Table, axis: str='sample') -> biom.Table:
    """ Convert feature table from frequencies to relative frequencies
    """
    return table.norm(axis=axis, inplace=False)
