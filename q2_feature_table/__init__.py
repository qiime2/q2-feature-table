# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._normalize import rarefy
from ._transform import (presence_absence, relative_frequency)
from ._summarize import (summarize, tabulate_seqs)
from ._merge import (merge, merge_seq_data, merge_taxa_data, overlap_methods)
from ._filter import (filter_samples, filter_features)
from ._core_features import core_features
from ._version import get_versions

__version__ = get_versions()['version']
del get_versions

__all__ = ['rarefy', 'presence_absence', 'relative_frequency', 'summarize',
           'merge', 'merge_seq_data', 'filter_samples', 'filter_features',
           'merge_taxa_data', 'tabulate_seqs', 'overlap_methods',
           'core_features']
