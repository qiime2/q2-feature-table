# ----------------------------------------------------------------------------
# Copyright (c) 2016-2019, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom


def subsample(table: biom.Table, subsampling_depth: int,
              axis: str) -> biom.Table:
    if axis == 'feature':
        # we are transposing the table due to biocore/biom-format#759
        table = table.transpose()

    if len(table.ids()) < subsampling_depth:
        raise ValueError('The subsampling depth exceeds the number of '
                         'elements on the desired axis. The maximum depth '
                         'is: %d.' % len(table.ids()))

    # the axis is always 'sample' due to the above transpose
    table = table.subsample(subsampling_depth, axis='sample', by_id=True)

    # the inverted axis is always observation due to the above transpose
    invaxis = 'observation'
    table.filter(lambda v, i, m: v.sum() > 0, axis=invaxis)

    if axis == 'feature':
        # reverse the transpose necessary due to biocore/biom-format#759
        table = table.transpose()

    if table.is_empty():
        raise ValueError('The subsampled table contains no samples or features'
                         ' (samples/features that sum to zero after filtering'
                         ' are automatically removed). It may be a good idea'
                         ' to double check that your table is valid/nonempty.')

    return table
