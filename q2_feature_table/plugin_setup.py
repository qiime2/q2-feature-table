# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import Plugin, Int, Properties

import q2_feature_table
from q2_types import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence, Phylogeny,
    Rooted, Unrooted)

plugin = Plugin(
    name='feature-table',
    version=q2_feature_table.__version__,
    website='https://github.com/qiime2-plugins/q2-feature-table',
    package='q2_feature_table'
)

plugin.methods.register_function(
    function=q2_feature_table.rarefy,
    # TODO use more restrictive primitive type for `depth`
    inputs={'table': FeatureTable[Frequency]},
    parameters={'counts_per_sample': Int},
    outputs=[('rarefied_table',
              FeatureTable[Frequency] % Properties('uniform-sampling'))],
    name='Rarefy table',
    description="Subsample counts from all samples without replacement so "
                "that the sum of counts in each sample is equal. Samples "
                "where the sum of counts is less than the requested counts "
                "per sample will be not be included in the resulting table."
)

plugin.methods.register_function(
    function=q2_feature_table.presence_absence,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency]},
    parameters={},
    outputs=[('presence_absence_table', FeatureTable[PresenceAbsence])],
    name="Convert to presence/absence",
    description="Convert frequencies to binary values indicating presence or "
                "absence of a feature in a sample."
)

plugin.methods.register_function(
    function=q2_feature_table.relative_frequency,
    inputs={'table': FeatureTable[Frequency]},
    parameters={},
    outputs=[
        ('relative_frequency_table',
         FeatureTable[RelativeFrequency] % Properties('uniform-sampling'))],
    name="Convert to relative frequencies",
    description="Convert frequencies to relative frequencies by dividing each "
                "count in a sample by the sum of counts in that sample."
)

plugin.visualizers.register_function(
    function=q2_feature_table.summarize,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={},
    name="Summarize table",
    description="Generate visual and tabular summaries of a feature table."
)

plugin.methods.register_function(
    function=q2_feature_table.filter,
    inputs={'table': FeatureTable[Frequency],
            'tree': Phylogeny[Rooted | Unrooted]},
    parameters={},
    outputs=[('filtered_table', FeatureTable[Frequency])],
    name="Remove features from table if they're not present in tree.",
    description=("This method is a placeholder and will be generalized to "
                 "support different types of filtering, including obtaining "
                 "ids from different data types and filtering on both axes. "
                 "See https://github.com/qiime2/q2-feature-table/issues/14 to "
                 "track progress on this.")
)
