# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._normalize import rarefy
from ._subsample import subsample
from ._transform import (presence_absence, relative_frequency, transpose)
from ._summarize import (summarize, tabulate_seqs)
from ._merge import (merge, merge_seqs, merge_taxa, overlap_methods)
from ._filter import (filter_samples, filter_features, filter_seqs)
from ._core_features import core_features
from ._group import group
from ._heatmap import (heatmap, heatmap_choices)
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

__all__ = ['rarefy', 'presence_absence', 'relative_frequency', 'transpose',
           'summarize', 'merge', 'merge_seqs', 'filter_samples',
           'filter_features', 'merge_taxa', 'tabulate_seqs', 'overlap_methods',
           'core_features', 'group', 'heatmap', 'heatmap_choices',
           'filter_seqs', 'subsample']
