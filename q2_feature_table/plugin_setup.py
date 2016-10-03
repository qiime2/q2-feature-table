# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime.plugin import Plugin, Int, Properties

import q2_feature_table
from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)
from q2_types.tree import Phylogeny, Rooted, Unrooted
from q2_types.feature_data import FeatureData, Sequence, Taxonomy

plugin = Plugin(
    name='feature-table',
    version=q2_feature_table.__version__,
    website='https://github.com/qiime2/q2-feature-table',
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

plugin.methods.register_function(
    function=q2_feature_table.merge,
    inputs={'table1': FeatureTable[Frequency],
            'table2': FeatureTable[Frequency]},
    parameters={},
    outputs=[
        ('merged_table', FeatureTable[Frequency])],
    name="Combine two tables",
    description="Combines a pair of feature tables which contain different "
                "samples, and which may or may not contain the same features."
)


plugin.methods.register_function(
    function=q2_feature_table.merge_seq_data,
    inputs={'data1': FeatureData[Sequence],
            'data2': FeatureData[Sequence]},
    parameters={},
    outputs=[
        ('merged_data', FeatureData[Sequence])],
    name="Combine two collections of feature sequences",
    description="Combines a pair of feature data objects which may or may not "
                "contain data for the same features. If different feature "
                "data is present for the same feature id in the two inputs, "
                "the data from the first (data1) will be propagated to the "
                "result."
)


plugin.methods.register_function(
    function=q2_feature_table.merge_taxa_data,
    inputs={'data1': FeatureData[Taxonomy],
            'data2': FeatureData[Taxonomy]},
    parameters={},
    outputs=[
        ('merged_data', FeatureData[Taxonomy])],
    name="Combine two collections of feature taxonomies",
    description="Combines a pair of feature data objects which may or may not "
                "contain data for the same features. If different feature "
                "data is present for the same feature id in the two inputs, "
                "the data from the first (data1) will be propagated to the "
                "result."
)

plugin.visualizers.register_function(
    function=q2_feature_table.summarize,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={},
    name="Summarize table",
    description="Generate visual and tabular summaries of a feature table."
)

plugin.visualizers.register_function(
    function=q2_feature_table.view_seq_data,
    inputs={'data': FeatureData[Sequence]},
    parameters={},
    name='View sequence associated with each feature',
    description="Generate tabular view of feature identifier to sequence "
                "mapping, including links to BLAST each sequence against "
                "the NCBI nt database."
)

plugin.visualizers.register_function(
    function=q2_feature_table.view_taxa_data,
    inputs={'data': FeatureData[Taxonomy]},
    parameters={},
    name='View taxonomy associated with each feature',
    description="Generate tabular view of feature identifier to taxonomic "
                "assignment mapping."
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
