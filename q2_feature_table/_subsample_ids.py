# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import Union

import biom
import numpy as np


def subsample_ids(table: biom.Table, subsampling_depth: int,
                  axis: str, seed: Union[int, str] = 1) -> biom.Table:
    if axis == 'feature':
        # we are transposing the table due to biocore/biom-format#759
        table = table.transpose()

    if len(table.ids()) < subsampling_depth:
        raise ValueError('The subsampling depth exceeds the number of '
                         'elements on the desired axis. The maximum depth '
                         'is: %d.' % len(table.ids()))

    # Generate a random seed if seed is None
    if seed == "random":
        rng = np.random.default_rng()
        seed = rng.integers(0, 2 ** 32 - 1)

    # the axis is always 'sample' due to the above transpose
    table = table.subsample(subsampling_depth, axis='sample',
                            by_id=True, seed=seed)

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
