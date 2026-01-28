# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

# Numpy recommends using at least 128 bits of entropy as a seed, and we are
# indirectly seeding numpy in this plugin.
RNG_MAX_SIZE = 2**128
