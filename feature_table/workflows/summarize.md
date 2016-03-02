---
name: summarize
inputs:
    table: feature_table.artifact_types.feature_table:FeatureTable
    sample_metadata:
outputs:
    null
---
## FeatureTable summarization

This workflow allows you to explore the data in your FeatureTable.

### Generate table summary data

```python
>>> import feature_table
>>> feature_table.summarize(table)
```

### Explore possible even sampling depths

```python
>>> %matplotlib inline
>>> feature_table.explore_even_sampling_depth(table, sample_metadata)
```
