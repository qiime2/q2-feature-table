# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def presence_absence(table):
    """ Convert feature table to presence/absence data
    """
    return table.pa(inplace=False)


def relative_frequency(table, axis='sample'):
    """ Convert feature table from frequencies to relative frequencies
    """
    return table.norm(axis=axis, inplace=False)
