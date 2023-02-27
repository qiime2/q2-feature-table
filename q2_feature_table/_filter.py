# ----------------------------------------------------------------------------
# Copyright (c) 2016-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import qiime2
import numpy as np
import pandas as pd


def _get_biom_filter_function(ids_to_keep, min_frequency, max_frequency,
                              min_nonzero, max_nonzero):
    ids_to_keep = set(ids_to_keep)
    if max_frequency is None:
        max_frequency = np.inf
    if max_nonzero is None:
        max_nonzero = np.inf

    def f(data_vector, id_, metadata):
        return (id_ in ids_to_keep) and \
               (min_frequency <= data_vector.sum() <= max_frequency) and \
               (min_nonzero <= (data_vector > 0).sum() <= max_nonzero)
    return f


_other_axis_map = {'sample': 'observation', 'observation': 'sample'}


def _filter_table(table, min_frequency, max_frequency, min_nonzero,
                  max_nonzero, metadata, where, axis, exclude_ids=False,
                  filter_opposite_axis=True):
    if min_frequency == 0 and max_frequency is None and min_nonzero == 0 and\
       max_nonzero is None and metadata is None and where is None and\
       exclude_ids is False:
        raise ValueError("No filtering was requested.")
    if metadata is None and where is not None:
        raise ValueError("Metadata must be provided if 'where' is "
                         "specified.")
    if metadata is None and exclude_ids is True:
        raise ValueError("Metadata must be provided if 'exclude_ids' "
                         "is True.")
    if metadata is not None:
        ids_to_keep = metadata.get_ids(where=where)
    else:
        ids_to_keep = table.ids(axis=axis)
    if exclude_ids is True:
        ids_to_keep = set(table.ids(axis=axis)) - set(ids_to_keep)

    filter_fn1 = _get_biom_filter_function(
        ids_to_keep, min_frequency, max_frequency, min_nonzero, max_nonzero)
    table.filter(filter_fn1, axis=axis, inplace=True)

    # filter on the opposite axis to remove any entities that now have a
    # frequency of zero
    if filter_opposite_axis:
        filter_fn2 = _get_biom_filter_function(
            ids_to_keep=table.ids(axis=_other_axis_map[axis]), min_frequency=0,
            max_frequency=None, min_nonzero=1, max_nonzero=None)
        table.filter(filter_fn2, axis=_other_axis_map[axis], inplace=True)


def filter_samples(table: biom.Table, min_frequency: int = 0,
                   max_frequency: int = None, min_features: int = 0,
                   max_features: int = None,
                   metadata: qiime2.Metadata = None, where: str = None,
                   exclude_ids: bool = False,
                   filter_empty_features: bool = True)\
                  -> biom.Table:
    _filter_table(table=table, min_frequency=min_frequency,
                  max_frequency=max_frequency, min_nonzero=min_features,
                  max_nonzero=max_features, metadata=metadata,
                  where=where, axis='sample', exclude_ids=exclude_ids,
                  filter_opposite_axis=filter_empty_features)

    return table


def filter_features(table: biom.Table, min_frequency: int = 0,
                    max_frequency: int = None, min_samples: int = 0,
                    max_samples: int = None,
                    metadata: qiime2.Metadata = None, where: str = None,
                    exclude_ids: bool = False,
                    filter_empty_samples: bool = True)\
                   -> biom.Table:
    _filter_table(table=table, min_frequency=min_frequency,
                  max_frequency=max_frequency, min_nonzero=min_samples,
                  max_nonzero=max_samples, metadata=metadata,
                  where=where, axis='observation', exclude_ids=exclude_ids,
                  filter_opposite_axis=filter_empty_samples)

    return table


def filter_seqs(data: pd.Series, table: biom.Table = None,
                metadata: qiime2.Metadata = None, where: str = None,
                exclude_ids: bool = False) -> pd.Series:
    if table is not None and metadata is not None:
        raise ValueError('Filtering with metadata and filtering with a table '
                         'are mutually exclusive.')
    elif table is None and metadata is None:
        raise ValueError('No filtering requested. Must provide either table '
                         'or metadata.')
    elif table is not None:
        ids_to_keep = table.ids(axis='observation')
    else:
        # Note, no need to check for missing feature IDs in the metadata,
        # because that is basically the point of this method.
        ids_to_keep = metadata.get_ids(where=where)

    if exclude_ids is True:
        ids_to_keep = set(data.index) - set(ids_to_keep)
    filtered = data[data.index.isin(ids_to_keep)]
    if filtered.empty is True:
        raise ValueError('All features were filtered out of the data.')
    return filtered


def filter_features_conditionally(table: biom.Table,
                                  abundance: float,
                                  prevalence: float,
                                  ) -> biom.Table:
    """
    A function to perform joint filtering because it makes life better
    """
    num_observations, num_samples = table.shape
    prevalence = prevalence * num_samples

    # Calculates the filtering parameters on the original table
    def _filter_f(values, id_, metadata):
        return (values >= abundance).sum() >= prevalence

    # Normalized the table to get the prevalance
    # Copy is because biom really wants to normalize the original table. By
    # copying and not using inplace, the original table is preserved.
    # Redundant, but better safe that sorry.
    table_norm = table.copy().norm(axis='sample', inplace=False)
    table_norm.filter(_filter_f, axis='observation', inplace=True)
    filter_ids = table_norm.ids(axis='observation')

    new_table = table.filter(filter_ids, axis='observation', inplace=False)

    return new_table
