# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from ._normalize import rarefy
from ._transform import (presence_absence, relative_frequency)


__version__ = '0.0.0-dev'

__all__ = ['rarefy', 'presence_absence', 'relative_frequency']
