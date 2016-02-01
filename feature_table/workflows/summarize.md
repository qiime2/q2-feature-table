---
name: summarize
type-imports:
    - feature_table.artifact_types:FeatureTable
inputs:
    table: FeatureTable
outputs: []    
---
## FeatureTable summarization

This workflow allows you to explore the data in your FeatureTable.

### Sample and feature counts

```python
>>> num_features, num_samples = table.shape
>>> print("Number of samples: %d" % num_samples)
>>> print("Number of features: %d" % num_features)
```

### Metadata summarization

```python
>>> from feature_table import metadata_summary
>>> sample_md_categories, feature_md_categories = metadata_summary(table)
>>> print("Sample metadata categories: %r" % sample_md_categories)
>>> print("Feature metadata categories: %r" % feature_md_categories)
```

### Sample count summary

```python
>>> from feature_table import count_summary
>>> sample_count_summary, sample_counts = count_summary(table, axis='sample')
```

```python
>>> sample_count_summary
```

```python
>>> import seaborn as sns
>>> %matplotlib inline
>>> ax = sns.distplot(sample_counts)
>>> ax.set_xlabel("Number of sequences per sample")
>>> ax.set_ylabel("Frequency")
>>> ax
```

### Feature count summary

```python
>>> feature_count_summary, feature_counts = count_summary(table, axis='observation')
```

```python
>>> feature_count_summary
```

```python
>>> ax = sns.distplot(feature_counts, hist=False)
>>> ax.set_xlabel("Number of sequences per feature")
>>> ax.set_ylabel("Frequency")
>>> ax.set_xscale('log')
>>> ax
```

### Explore possible even sampling depths

```python
>>> from feature_table import max_count_even_sampling_depth
>>> max_count_sampling_depth = max_count_even_sampling_depth(sample_counts)
>>> print("Rarefaction depth of %d will retain the largest number of sequences." % max_count_sampling_depth)
```
