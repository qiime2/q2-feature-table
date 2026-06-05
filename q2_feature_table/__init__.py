# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._normalize import rarefy, normalize
from ._subsample_ids import subsample_ids
from ._transform import (presence_absence, relative_frequency, transpose)
from ._summarize import (_summarize, tabulate_seqs,
                         tabulate_sample_frequencies,
                         tabulate_feature_frequencies, summarize)
from ._merge import (merge, merge_seqs, merge_taxa, overlap_methods)
from ._filter import (filter_samples, filter_features, filter_seqs,
                      filter_features_conditionally)
from ._split import split
from ._core_features import core_features
from ._group import group
from ._rename import rename_ids
from ._heatmap import (heatmap, heatmap_choices)
from ._multiply import (multiply_tables, _multiply_tables,
                        _multiply_tables_pa, _multiply_tables_relative)

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = '0.0.0+notfound'

__all__ = ['rarefy', 'presence_absence', 'relative_frequency', 'transpose',
           '_summarize', 'merge', 'merge_seqs', 'filter_samples',
           'filter_features', 'merge_taxa', 'tabulate_seqs', 'overlap_methods',
           'core_features', 'group', 'heatmap', 'heatmap_choices',
           'filter_seqs', 'subsample_ids', 'rename_ids',
           'filter_features_conditionally', 'split',
           'tabulate_feature_frequencies', 'tabulate_sample_frequencies',
           'summarize', 'normalize', 'multiply_tables', '_multiply_tables',
           '_multiply_tables_pa', '_multiply_tables_relative']
