# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import Plugin, Int, Properties, Metadata, Str

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
                   'doi:10.1186/2047-217X-1-7'),
    short_description=('Plugin for working with sample x feature tables.'),
    description=('This is a QIIME 2 plugin supporting operations on sample '
                 'x feature tables, such as filtering, merging, and '
                 'transforming tables.')
)

plugin.methods.register_function(
    function=q2_feature_table.rarefy,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int},
    outputs=[('rarefied_table',
              FeatureTable[Frequency] % Properties('uniform-sampling'))],
    input_descriptions={'table': 'The feature table to be rarefied.'},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'rarefied to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table.')
    },
    output_descriptions={
        'rarefied_table': 'The resulting rarefied feature table.'
    },
    name='Rarefy table',
    description=("Subsample frequencies from all samples without replacement "
                 "so that the sum of frequencies in each sample is equal to "
                 "sampling-depth.")
)

plugin.methods.register_function(
    function=q2_feature_table.presence_absence,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency]},
    parameters={},
    outputs=[('presence_absence_table', FeatureTable[PresenceAbsence])],
    input_descriptions={
        'table': ('The feature table to be converted into presence/absence '
                  'abundances.')
    },
    parameter_descriptions={},
    output_descriptions={
        'presence_absence_table': ('The resulting presence/absence feature '
                                   'table.')
    },
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
    input_descriptions={
        'table': 'The feature table to be converted into relative frequencies.'
    },
    parameter_descriptions={},
    output_descriptions={
        'relative_frequency_table': ('The resulting relative frequency '
                                     'feature table.')
    },
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
    input_descriptions={
        'table1': 'The first feature table to be merged.',
        'table2': 'The second feature table to be merged.',
    },
    parameter_descriptions={},
    output_descriptions={
        'merged_table': ('The resulting merged feature table containing all '
                         'features from table1 and table2.')
    },
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
    input_descriptions={
        'data1': 'The first collection of feature sequences to be merged.',
        'data2': 'The second collection of feature sequences to be merged.',
    },
    parameter_descriptions={},
    output_descriptions={
        'merged_data': ('The resulting collection of feature sequences '
                        'containing all feature sequences from data1 and '
                        'data2.')
    },
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
    input_descriptions={
        'data1': 'The first collection of feature taxonomies to be merged.',
        'data2': 'The second collection of feature taxonomies to be merged.',
    },
    parameter_descriptions={},
    output_descriptions={
        'merged_data': ('The resulting collection of feature taxonomies '
                        'containing all feature taxonomies from data1 and '
                        'data2.')
    },
    name="Combine two collections of feature taxonomies",
    description="Combines a pair of feature data objects which may or may not "
                "contain data for the same features. If different feature "
                "data is present for the same feature id in the two inputs, "
                "the data from the first (data1) will be propagated to the "
                "result."
)

plugin.methods.register_function(
    function=q2_feature_table.filter_samples,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'min_frequency': Int,
                'max_frequency': Int,
                'min_features': Int,
                'max_features': Int,
                'sample_metadata': Metadata,
                'where': Str},
    outputs=[('filtered_table', FeatureTable[Frequency])],
    input_descriptions={
        'table': 'The feature table from which samples should be filtered.'
    },
    parameter_descriptions={
        'min_frequency': ('The minimum total frequency that a sample must '
                          'have to be retained.'),
        'max_frequency': ('The maximum total frequency that a sample can '
                          'have to be retained. If no value is provided '
                          'this will default to infinity (i.e., no maximum '
                          'frequency filter will be applied).'),
        'min_features': ('The minimum number of features that a sample must '
                         'have to be retained.'),
        'max_features': ('The maximum number of features that a sample can '
                         'have to be retained. If no value is provided '
                         'this will default to infinity (i.e., no maximum '
                         'feature filter will be applied).'),
        'sample_metadata': 'Sample metadata used in conjuction with `where` '
                           'parameter to select samples to retain.',
        'where': 'SQLite WHERE clause specifying sample metadata criteria '
                 'that must be met to be included in the filtered feature '
                 'table. If not provided, all samples in `sample_metadata` '
                 'that are also in the feature table will be retained.'
    },
    output_descriptions={
        'filtered_table': 'The resulting feature table filtered by sample.'
    },
    name="Filter samples from table.",
    description="Filter samples from table based on frequency and/or "
                "metadata. Any features with a frequency of zero after sample "
                "filtering will also be removed. See the filtering tutorial "
                "on https://docs.qiime2.org for additional details."
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
    outputs=[('filtered_table', FeatureTable[Frequency])],
    input_descriptions={
        'table': 'The feature table from which features should be filtered.'
    },
    parameter_descriptions={
        'min_frequency': ('The minimum total frequency that a feature must '
                          'have to be retained.'),
        'max_frequency': ('The maximum total frequency that a feature can '
                          'have to be retained. If no value is provided '
                          'this will default to infinity (i.e., no maximum '
                          'frequency filter will be applied).'),
        'min_samples': ('The minimum number of samples that a feature must '
                        'be observed in to be retained.'),
        'max_samples': ('The maximum number of samples that a feature can '
                        'be observed in to be retained. If no value is '
                        'provided this will default to infinity (i.e., no '
                        'maximum sample filter will be applied).'),
        'feature_metadata': 'Feature metadata used in conjuction with `where` '
                            'parameter to select features to retain.',
        'where': 'SQLite WHERE clause specifying feature metadata criteria '
                 'that must be met to be included in the filtered feature '
                 'table. If not provided, all features in `feature_metadata` '
                 'that are also in the feature table will be retained.'
    },
    output_descriptions={
        'filtered_table': 'The resulting feature table filtered by feature.'
    },
    name="Filter features from table.",
    description="Filter features from table based on frequency and/or "
                "metadata. Any samples with a frequency of zero after feature "
                "filtering will also be removed. See the filtering tutorial "
                "on https://docs.qiime2.org for additional details."
)

plugin.visualizers.register_function(
    function=q2_feature_table.summarize,
    inputs={'table': FeatureTable[Frequency | RelativeFrequency |
                                  PresenceAbsence]},
    parameters={'sample_metadata': Metadata},
    input_descriptions={'table': 'The feature table to be summarized.'},
    parameter_descriptions={'sample_metadata': 'The sample metadata.'},
    name="Summarize table",
    description="Generate visual and tabular summaries of a feature table."
)

plugin.visualizers.register_function(
    function=q2_feature_table.tabulate_seqs,
    inputs={'data': FeatureData[Sequence]},
    parameters={},
    input_descriptions={'data': 'The feature sequences to be tabulated.'},
    parameter_descriptions={},
    name='View sequence associated with each feature',
    description="Generate tabular view of feature identifier to sequence "
                "mapping, including links to BLAST each sequence against "
                "the NCBI nt database."
)
