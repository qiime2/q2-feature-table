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
>>> #table = biom.util.biom_open(open(biom_path))
... table = biom.load_table(biom_path)
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
```

### Metadata summarization

```python
>>> sample_md_categories, feature_md_categories = ft.metadata_summary(table)
>>> print("Sample metadata categories: %r" % sample_md_categories)
>>> print("Feature metadata categories: %r" % feature_md_categories)
```

### Sample count summary

```python
>>> sample_count_summary, sample_counts = ft.count_summary(table, axis='sample')
```

```python
>>> sample_count_summary
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
>>> max_count_sampling_depth = ft.max_count_even_sampling_depth(sample_counts)
>>> print("Rarefaction depth of %d will retain the largest number of sequences." % max_count_sampling_depth)
```
