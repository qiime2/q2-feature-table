# ----------------------------------------------------------------------------
# Copyright (c) 2016-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import (Plugin, Int, Float, Range, Metadata, Str, Bool,
                           Choices, MetadataColumn, Categorical, List,
                           Citations)

import q2_feature_table
from q2_types.feature_table import (
    FeatureTable, Frequency, RelativeFrequency, PresenceAbsence)
from q2_types.feature_data import FeatureData, Sequence, Taxonomy

citations = Citations.load('citations.bib', package='q2_feature_table')
plugin = Plugin(
    name='feature-table',
    version=q2_feature_table.__version__,
    website='https://github.com/qiime2/q2-feature-table',
    package='q2_feature_table',
    short_description=('Plugin for working with sample by feature tables.'),
    description=('This is a QIIME 2 plugin supporting operations on sample '
                 'by feature tables, such as filtering, merging, and '
                 'transforming tables.')
)

plugin.methods.register_function(
    function=q2_feature_table.rarefy,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'sampling_depth': Int % Range(1, None),
                'with_replacement': Bool},
    outputs=[('rarefied_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be rarefied.'},
    parameter_descriptions={
        'sampling_depth': ('The total frequency that each sample should be '
                           'rarefied to. Samples where the sum of frequencies '
                           'is less than the sampling depth will be not be '
                           'included in the resulting table unless '
                           'subsampling is performed with replacement.'),
        'with_replacement': ('Rarefy with replacement by sampling from the '
                             'multinomial distribution instead of rarefying '
                             'without replacement.')
    },
    output_descriptions={
        'rarefied_table': 'The resulting rarefied feature table.'
    },
    name='Rarefy table',
    description=("Subsample frequencies from all samples so that the sum of "
                 "frequencies in each sample is equal to sampling-depth."),
    citations=[citations['Weiss2017']]
)

plugin.methods.register_function(
    function=q2_feature_table.subsample,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'subsampling_depth': Int % Range(1, None),
                'axis': Str % Choices(['sample', 'feature'])},
    outputs=[('sampled_table', FeatureTable[Frequency])],
    input_descriptions={'table': 'The feature table to be sampled.'},
    parameter_descriptions={
        'subsampling_depth': ('The total number of samples or features to be '
                              'randomly sampled. Samples or features that are '
                              'reduced to a zero sum will not be included in '
                              'the resulting table.'),
        'axis': ('The axis to sample over. If "sample" then samples will be '
                 'randomly selected to be retained. If "feature" then '
                 'a random set of features will be selected to be retained.')
    },
    output_descriptions={
        'sampled_table': 'The resulting subsampled feature table.'
    },
    name='Subsample table',
    description=("Randomly pick samples or features, without replacement, "
                 "from the table.")
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
         FeatureTable[RelativeFrequency])],
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
    function=q2_feature_table.group,
    inputs={'table': FeatureTable[Frequency]},
    parameters={
        'mode': Str % Choices({'sum', 'median-ceiling', 'mean-ceiling'}),
        'metadata': MetadataColumn[Categorical],
        'axis': Str % Choices({'sample', 'feature'})
    },
    outputs=[
        ('grouped_table', FeatureTable[Frequency])
    ],
    input_descriptions={
        'table': 'The table to group samples or features on.'
    },
    parameter_descriptions={
        'mode': 'How to combine samples or features within a group. `sum` '
                'will sum the frequencies across all samples or features '
                'within a group; `mean-ceiling` will take the ceiling of the '
                'mean of these frequencies; `median-ceiling` will take the '
                'ceiling of the median of these frequencies.',
        'metadata': 'A column defining the groups. Each unique value will '
                    'become a new ID for the table on the given `axis`.',
        'axis': 'Along which axis to group. Each ID in the given axis must '
                'exist in `metadata`.'
    },
    output_descriptions={
        'grouped_table': 'A table that has been grouped along the given '
                         '`axis`. IDs on that axis are replaced by values in '
                         'the `metadata` column.'
    },
    name="Group samples or features by a metadata column",
    description="Group samples or features in a feature table using metadata "
                "to define the mapping of IDs to a group."
)

plugin.methods.register_function(
    function=q2_feature_table.merge,
    inputs={'tables': List[FeatureTable[Frequency]]},
    parameters={
        'overlap_method': Str % Choices(q2_feature_table.overlap_methods()),
    },
    outputs=[
        ('merged_table', FeatureTable[Frequency])],
    input_descriptions={
        'tables': 'The collection of feature tables to be merged.',
    },
    parameter_descriptions={
        'overlap_method': 'Method for handling overlapping ids.',
    },
    output_descriptions={
        'merged_table': ('The resulting merged feature table.'),
    },
    name="Combine multiple tables",
    description="Combines feature tables using the `overlap_method` provided."
)


plugin.methods.register_function(
    function=q2_feature_table.merge_seqs,
    inputs={'data': List[FeatureData[Sequence]]},
    parameters={},
    outputs=[
        ('merged_data', FeatureData[Sequence])],
    input_descriptions={
        'data': 'The collection of feature sequences to be merged.',
    },
    parameter_descriptions={},
    output_descriptions={
        'merged_data': ('The resulting collection of feature sequences '
                        'containing all feature sequences provided.')
    },
    name="Combine collections of feature sequences",
    description="Combines feature data objects which may or may not "
                "contain data for the same features. If different feature "
                "data is present for the same feature id in the inputs, "
                "the data from the first will be propagated to the result."
)


plugin.methods.register_function(
    function=q2_feature_table.merge_taxa,
    inputs={'data': List[FeatureData[Taxonomy]]},
    parameters={},
    outputs=[
        ('merged_data', FeatureData[Taxonomy])],
    input_descriptions={
        'data': 'The collection of feature taxonomies to be merged.',
    },
    parameter_descriptions={},
    output_descriptions={
        'merged_data': ('The resulting collection of feature taxonomies '
                        'containing all feature taxonomies provided.')
    },
    name="Combine collections of feature taxonomies",
    description="Combines a pair of feature data objects which may or may not "
                "contain data for the same features. If different feature "
                "data is present for the same feature id in the inputs, "
                "the data from the first will be propagated to the result."
)

plugin.methods.register_function(
    function=q2_feature_table.filter_samples,
    inputs={'table': FeatureTable[Frequency]},
    parameters={'min_frequency': Int,
                'max_frequency': Int,
                'min_features': Int,
                'max_features': Int,
                'metadata': Metadata,
                'where': Str,
                'exclude_ids': Bool},
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
        'metadata': 'Sample metadata used with `where` parameter when '
                    'selecting samples to retain, or with `exclude_ids` '
                    'when selecting samples to discard.',
        'where': 'SQLite WHERE clause specifying sample metadata criteria '
                 'that must be met to be included in the filtered feature '
                 'table. If not provided, all samples in `metadata` that are '
                 'also in the feature table will be retained.',
        'exclude_ids': 'If true, the samples selected by `metadata` or '
                       '`where` parameters will be excluded from the filtered '
                       'table instead of being retained.'
    },
    output_descriptions={
        'filtered_table': 'The resulting feature table filtered by sample.'
    },
    name="Filter samples from table",
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
                'metadata': Metadata,
                'where': Str,
                'exclude_ids': Bool},
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
        'metadata': 'Feature metadata used with `where` parameter when '
                    'selecting features to retain, or with `exclude_ids` '
                    'when selecting features to discard.',
        'where': 'SQLite WHERE clause specifying feature metadata criteria '
                 'that must be met to be included in the filtered feature '
                 'table. If not provided, all features in `metadata` that are '
                 'also in the feature table will be retained.',
        'exclude_ids': 'If true, the features selected by `metadata` or '
                       '`where` parameters will be excluded from the filtered '
                       'table instead of being retained.'
    },
    output_descriptions={
        'filtered_table': 'The resulting feature table filtered by feature.'
    },
    name="Filter features from table",
    description="Filter features from table based on frequency and/or "
                "metadata. Any samples with a frequency of zero after feature "
                "filtering will also be removed. See the filtering tutorial "
                "on https://docs.qiime2.org for additional details."
)

plugin.methods.register_function(
    function=q2_feature_table.filter_seqs,
    inputs={
        'data': FeatureData[Sequence],
        'table': FeatureTable[Frequency],
    },
    parameters={
        'metadata': Metadata,
        'where': Str,
        'exclude_ids': Bool
    },
    outputs=[('filtered_data', FeatureData[Sequence])],
    input_descriptions={
        'data': 'The sequences from which features should be filtered.',
        'table': 'Table containing feature ids used for id-based filtering.'
    },
    parameter_descriptions={
        'metadata': 'Feature metadata used for id-based filtering, with '
                    '`where` parameter when selecting features to retain, or '
                    'with `exclude_ids` when selecting features to discard.',
        'where': 'SQLite WHERE clause specifying feature metadata criteria '
                 'that must be met to be included in the filtered feature '
                 'table. If not provided, all features in `metadata` that are '
                 'also in the sequences will be retained.',
        'exclude_ids': 'If true, the features selected by the `metadata` '
                       '(with or without the `where` parameter) or `table` '
                       'parameter will be excluded from the filtered '
                       'sequences instead of being retained.'
    },
    output_descriptions={
        'filtered_data': 'The resulting filtered sequences.'
    },
    name="Filter features from sequences",
    description="Filter features from sequences based on a feature table or "
                "metadata. See the filtering tutorial on "
                "https://docs.qiime2.org for additional details. This method "
                "can filter based on ids in a table or a metadata file, but "
                "not both (i.e., the table and metadata options are mutually "
                "exclusive)."
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
                "the NCBI nt database.",
    citations=[citations['NCBI'], citations['NCBI-BLAST']]
)

plugin.visualizers.register_function(
    function=q2_feature_table.core_features,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'min_fraction': Float % Range(0.0, 1.0, inclusive_start=False),
        'max_fraction': Float % Range(0.0, 1.0, inclusive_end=True),
        'steps': Int % Range(2, None)
    },
    name='Identify core features in table',
    description=('Identify "core" features, which are features observed in a '
                 'user-defined fraction of the samples. Since the core '
                 'features are a function of the fraction of samples that the '
                 'feature must be observed in to be considered core, this is '
                 'computed over a range of fractions defined by the '
                 '`min_fraction`, `max_fraction`, and `steps` parameters.'),
    input_descriptions={
        'table': 'The feature table to use in core features calculations.'
    },
    parameter_descriptions={
        'min_fraction': 'The minimum fraction of samples that a feature must '
                        'be observed in for that feature to be considered a '
                        'core feature.',
        'max_fraction': 'The maximum fraction of samples that a feature must '
                        'be observed in for that feature to be considered a '
                        'core feature.',
        'steps': 'The number of steps to take between `min_fraction` and '
                 '`max_fraction` for core features calculations. This '
                 'parameter has no effect if `min_fraction` and '
                 '`max_fraction` are the same value.'
    }
)


plugin.visualizers.register_function(
    function=q2_feature_table.heatmap,
    inputs={
        'table': FeatureTable[Frequency]
    },
    parameters={
        'metadata': MetadataColumn[Categorical],
        'normalize': Bool,
        'title': Str,
        'metric': Str % Choices(q2_feature_table.heatmap_choices['metric']),
        'method': Str % Choices(q2_feature_table.heatmap_choices['method']),
        'cluster': Str % Choices(q2_feature_table.heatmap_choices['cluster']),
        'color_scheme': Str % Choices(
            q2_feature_table.heatmap_choices['color_scheme']),
    },
    name='Generate a heatmap representation of a feature table',
    description='Generate a heatmap representation of a feature table with '
                'optional clustering on both the sample and feature axes.\n\n'
                'Tip: To generate a heatmap containing taxonomic annotations, '
                'use `qiime taxa collapse` to collapse the feature table at '
                'the desired taxonomic level.',
    input_descriptions={
        'table': 'The feature table to visualize.'
    },
    parameter_descriptions={
        'metadata': 'Annotate the sample IDs with these metadata values. '
                    'When metadata is present and `cluster`=\'feature\', '
                    'samples will be sorted by the metadata values.',
        'normalize': 'Normalize the feature table by adding a psuedocount '
                     'of 1 and then taking the log10 of the table.',
        'title': 'Optional custom plot title.',
        'metric': 'Metrics exposed by seaborn (see http://seaborn.pydata.org/'
                  'generated/seaborn.clustermap.html#seaborn.clustermap for '
                  'more detail).',
        'method': 'Clustering methods exposed by seaborn (see http://seaborn.'
                  'pydata.org/generated/seaborn.clustermap.html#seaborn.clust'
                  'ermap for more detail).',
        'cluster': 'Specify which axes to cluster.',
        'color_scheme': 'The matplotlib colorscheme to generate the heatmap '
                        'with.',
    },
    citations=[citations['Hunter2007Matplotlib']]
)
