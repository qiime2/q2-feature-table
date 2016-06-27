# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import Plugin, Int

import q2_feature_table
from q2_types import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)

plugin = Plugin(
    name='feature-table',
    version=q2_feature_table.__version__,
    website='https://github.com/qiime2-plugins/q2-feature-table',
    package='q2_feature_table'
)

# TODO create decorator for promoting functions to methods. This info would
# be moved to the decorator calls.
plugin.methods.register_function(
    function=q2_feature_table.rarefy,
    # TODO use more restrictive primitive type for `depth`
    inputs={'table': FeatureTable[Frequency]},
    parameters={'depth': Int},
    outputs=[('rarefied_table', FeatureTable[Frequency])],
    name='Rarefaction',
    description="Let's rarefy!"
)

plugin.methods.register_function(
    function=q2_feature_table.presence_absence,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency]},
    parameters={},
    outputs=[('presence_absence_table', FeatureTable[PresenceAbsence])],
    name='Convert to presence/absence',
    description="Let's convert to presence/absence!"
)

plugin.methods.register_function(
    function=q2_feature_table.relative_frequency,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[('relative_frequency_table', FeatureTable[RelativeFrequency])],
    name='Convert to relative frequencies',
    description="Let's convert to relative frequencies!"
)

plugin.visualizers.register_function(
    function=q2_feature_table.summarize,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={},
    name='Summarize feature table',
    description='Generate visual and tabular summaries of a feature table.'
)
