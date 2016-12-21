# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import sqlite3

import biom
import qiime2
import numpy as np


def _ids_where(metadata, where):
    df = metadata.to_dataframe()

    conn = sqlite3.connect(':memory:')
    conn.row_factory = lambda cursor, row: row[0]

    df.to_sql('metadata', conn)
    id_column = df.index.name

    c = conn.cursor()

    # In general we wouldn't want to format our query in this way because
    # it leaves us open to sql injection, but it seems acceptable here for
    # a few reasons:
    # 1) This is a throw-away database which we're just creating to have
    #    access to the query language, so any malicious behavior wouldn't
    #    impact any data that isn't temporary
    # 2) The substitution syntax recommended in the docs doesn't allow
    #    us to specify complex where statements, which is what we need to do
    #    here. For example, we need to specify things like:
    #        WHERE Subject='subject-1' AND SampleType='gut'
    #    but their qmark/named-style syntaxes only supports substition of
    #    variables, such as:
    #        WHERE Subject=?
    # 3) sqlite3.Cursor.execute will only execute a single statement so
    #    inserting multiple statements (e.g., "Subject='subject-1'; DROP...")
    #    will result in an OperationalError being raised.
    query = ('SELECT "{0}" FROM metadata WHERE {1} GROUP BY "{0}" '
             'ORDER BY "{0}";'.format(id_column, where))

    try:
        c.execute(query)
    except sqlite3.OperationalError:
        conn.close()
        raise ValueError("Selection of ids failed with query:\n %s"
                         % query)

    ids = c.fetchall()
    conn.close()
    return ids


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


def _filter(table, min_frequency, max_frequency, min_nonzero, max_nonzero,
            metadata, where, axis):
    if min_frequency == 0 and max_frequency is None and min_nonzero == 0 and\
       max_nonzero is None and metadata is None and where is None:
        raise ValueError("No filtering was requested.")
    if metadata is None and where is not None:
        raise ValueError("Metadata must be provided if 'where' is "
                         "specified.")

    if where is not None:
        ids_to_keep = _ids_where(metadata, where)
    elif metadata is not None:
        ids_to_keep = metadata.to_dataframe().index
    else:
        ids_to_keep = table.ids(axis=axis)

    filter_fn1 = _get_biom_filter_function(
        ids_to_keep, min_frequency, max_frequency, min_nonzero, max_nonzero)
    table.filter(filter_fn1, axis=axis, inplace=True)

    # filter on the opposite axis to remove any entities that now have a
    # frequency of zero
    filter_fn2 = _get_biom_filter_function(
        ids_to_keep=table.ids(axis=_other_axis_map[axis]), min_frequency=1,
        max_frequency=None, min_nonzero=0, max_nonzero=None)
    table.filter(filter_fn2, axis=_other_axis_map[axis], inplace=True)


def filter_samples(table: biom.Table, min_frequency: int=0,
                   max_frequency: int=None, min_features: int=0,
                   max_features: int=None,
                   sample_metadata: qiime2.Metadata=None, where: str=None)\
                  -> biom.Table:
    _filter(table=table, min_frequency=min_frequency,
            max_frequency=max_frequency, min_nonzero=min_features,
            max_nonzero=max_features, metadata=sample_metadata,
            where=where, axis='sample')

    return table


def filter_features(table: biom.Table, min_frequency: int=0,
                    max_frequency: int=None, min_samples: int=0,
                    max_samples: int=None,
                    feature_metadata: qiime2.Metadata=None, where: str=None)\
                   -> biom.Table:
    _filter(table=table, min_frequency=min_frequency,
            max_frequency=max_frequency, min_nonzero=min_samples,
            max_nonzero=max_samples, metadata=feature_metadata,
            where=where, axis='observation')

    return table
