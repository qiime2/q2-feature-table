# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime
from qiime.plugin import Plugin, Int, Properties, Metadata, Str

import q2_feature_table
from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)
from q2_types.feature_data import FeatureData, Sequence, Taxonomy

plugin = Plugin(
    name='feature-table',
    version=q2_feature_table.__version__,
    website='https://github.com/qiime2/q2-feature-table',
    package='q2_feature_table',
    citation_text=('The Biological Observation Matrix (BIOM) format or: how '
                   'I learned to stop worrying and love the ome-ome. '
                   'Daniel McDonald, Jose C Clemente, Justin Kuczynski, '
                   'Jai Ram Rideout, Jesse Stombaugh, Doug Wendel, Andreas '
                   'Wilke, Susan Huse, John Hufnagle, Folker Meyer, Rob '
                   'Knight and J Gregory Caporaso. GigaScience 1:7 (2012).'
                   'doi:10.1186/2047-217X-1-7')
)

plugin.methods.register_function(
    function=q2_feature_table.rarefy,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int},
    outputs=[('rarefied_table',
              FeatureTable[Frequency] % Properties('uniform-sampling'))],
    name='Rarefy table',
    description="Subsample frequencies from all samples without replacement "
                "so that the sum of frequencies in each sample is equal to "
                "sampling-depth. Samples where the sum of frequencies is less "
                "than sampling-depth will be not be included in the resulting "
                "table."
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
                "frequency in a sample by the sum of frequencies in that "
                "sample."
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

_where_description = ("The where parameter takes a SQLite WHERE clause. "
                      "See the table filtering tutorial for additional "
                      "detail: https://docs.qiime2.org/%s/tutorials/"
                      "table-filtering.html" % qiime.__version__)

plugin.methods.register_function(
    function=q2_feature_table.filter_samples,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'min_frequency': Int,
                'max_frequency': Int,
                'min_features': Int,
                'max_features': Int,
                'sample_metadata': Metadata,
                'where': Str},
    outputs={'filtered_table': FeatureTable[Frequency]},
    name="Filter samples from table.",
    description="Filter samples from table based on frequency and/or "
                "metadata. Any features with a frequency of zero after sample "
                "filtering will also be removed. If no value(s) are provided "
                "for max_frequency or max_features, they will default to "
                "infinity (i.e., no maximum frequency and/or feature filter "
                "will be applied).\n\n%s" % _where_description
)

plugin.methods.register_function(
    function=q2_feature_table.filter_features,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'min_frequency': Int,
                'max_frequency': Int,
                'min_samples': Int,
                'max_samples': Int,
                'feature_metadata': Metadata,
                'where': Str},
    outputs={'filtered_table': FeatureTable[Frequency]},
    name="Filter features from table.",
    description="Filter features from table based on frequency and/or "
                "metadata. Any samples with a frequency of zero after feature "
                "filtering will also be removed. If no value(s) are provided "
                "for max_frequency and/or max_samples, they will default to "
                "infinity (i.e., no maximum frequency and/or sample filter "
                "will be applied).\n\n%s" % _where_description
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
    function=q2_feature_table.tabulate_seqs,
    inputs={'data': FeatureData[Sequence]},
    parameters={},
    name='View sequence associated with each feature',
    description="Generate tabular view of feature identifier to sequence "
                "mapping, including links to BLAST each sequence against "
                "the NCBI nt database."
)
