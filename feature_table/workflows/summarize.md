---
name: summarize
inputs:
    table: feature_table.artifact_types.feature_table:FeatureTable
    sample_metadata: null
outputs:
    null
---
## FeatureTable summarization

```python
>>> # this is temporarily in place for testing purposes
... import biom
>>> biom_path = "/Users/caporaso/code/q2d2/example-data/keyboard/q191/otu-table.tsv"
>>> table = biom.load_table(biom_path)
```

```python
>>> # config - we should recommend a standard for this
... import feature_table._summarize as ft
>>> import seaborn as sns
...
>>> %matplotlib inline
```

This workflow allows you to explore the data in your FeatureTable.

### Sample and feature counts

```python
>>> num_features, num_samples = table.shape
>>> print("Number of samples: %d" % num_samples)
>>> print("Number of features: %d" % num_features)
Number of samples: 104
Number of features: 658
```

### Metadata summarization

```python
>>> sample_md_categories, feature_md_categories = ft.metadata_summary(table)
>>> print("Sample metadata categories: %r" % sample_md_categories)
>>> print("Feature metadata categories: %r" % feature_md_categories)
Sample metadata categories: []
Feature metadata categories: ['taxonomy']
```

### Sample count summary

```python
>>> sample_count_summary, sample_counts = ft.count_summary(table, axis='sample')
/Users/caporaso/Dropbox/code/feature-table/feature_table/_summarize.py:28: FutureWarning: sort is deprecated, use sort_values(inplace=True) for for INPLACE sorting
  summary.sort()
```

```python
>>> sample_count_summary
Minimum count     77.000000
1st quartile     277.500000
Median count     314.000000
Mean count       320.634615
3rd quartile     377.250000
Maximum count    634.000000
dtype: float64
```

```python
>>> ax = sns.distplot(sample_counts)
>>> ax.set_xlabel("Number of sequences per sample")
>>> ax.set_ylabel("Frequency")
>>> ax
```

### Feature count summary

```python
>>> feature_count_summary, feature_counts = ft.count_summary(table, axis='observation')
/Users/caporaso/Dropbox/code/feature-table/feature_table/_summarize.py:28: FutureWarning: sort is deprecated, use sort_values(inplace=True) for for INPLACE sorting
  summary.sort()
```

```python
>>> feature_count_summary
Minimum count        2.000000
1st quartile         2.000000
Median count         4.000000
3rd quartile        11.000000
Mean count          50.677812
Maximum count    16248.000000
dtype: float64
```

```python
>>> ax = sns.distplot(feature_counts)
>>> ax.set_xlabel("Number of sequences per feature")
>>> ax.set_ylabel("Frequency")
>>> ax.set(xscale="log")
>>> ax
```

### Explore possible even sampling depths

```python
>>> max_count_sampling_depth = ft.max_count_even_sampling_depth(sample_counts)
>>> print("Rarefaction depth of %d will retain the largest number of sequences." % max_count_sampling_depth)
Rarefaction depth of 265 will retain the largest number of sequences.
```

```python

```
