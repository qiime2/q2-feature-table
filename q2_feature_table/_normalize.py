# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def rarefy(table: biom.Table, sampling_depth: int,
           with_replacement: bool=False) -> biom.Table:
    table = table.subsample(sampling_depth, axis='sample', by_id=False,
                            with_replacement=with_replacement)

    if table.is_empty():
        raise ValueError('The rarefied table contains no samples or features. '
                         'Verify your table is valid and that you provided a '
                         'shallow enough sampling depth.')

    return table
