# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def presence_absence(table: biom.Table) -> biom.Table:
    """ Convert feature table in-place to presence/absence data
    """
    table.pa(inplace=True)
    return table


def relative_frequency(table: biom.Table) -> biom.Table:
    """ Convert feature table in-place from frequencies to relative frequencies
    """
    table.norm(axis='sample', inplace=True)
    return table
