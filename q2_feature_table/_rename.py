# ----------------------------------------------------------------------------
# Copyright (c) 2016-2020, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import warnings

import biom
import qiime2


def _generate_new_names(old_ids, rename, strict, verbose=False):
    """
    Checks the list of ids to be renamed
    """
    # all old ids need to be unique
    if (rename.value_counts() > 1).any():
        raise ValueError('All new ids must be unique.\n'
                         'Try the "group" method in this plugin if you want '
                         'to combine multiple ids in the same table.')
    old_disjoint = list(set(old_ids) - set(rename.index))
    new_disjoint = list(set(rename.index) - set(old_ids))

    if (len(old_disjoint) > 0) & strict:
        missing_rename = '\n    '.join(old_disjoint)
        raise ValueError(
            'There are ids in the table which do not have new names.\n'
            'Either turn off strict mode or provide new names for all ids.'
            '\n'
            'The following ids are missing new names:\n'
            '    %s' % missing_rename
            )
    elif (len(old_disjoint) > 0) & verbose:
        missing_rename = '\n    '.join(old_disjoint)
        warnings.warn(UserWarning(
            'There are ids in the original table which do not have new names.'
            '\n'
            'The following ids will not be included:\n'
            '   %s' % missing_rename)
        )
    if (len(new_disjoint) > 0) & verbose:
        missing_rename = '\n    '.join(new_disjoint)
        warnings.warn(UserWarning(
            'There are ids supplied for renaming that are not in the table.\n'
            'The following ids will not be mapped:\n'
            '   %s' % missing_rename)
        )

    new_ids = {id_: rename.to_dict().get(id_, id_) for id_ in old_ids}
    return new_ids


def rename_ids(table: biom.Table,
               metadata: qiime2.CategoricalMetadataColumn,
               axis: str = 'sample',
               strict: bool = False)\
               -> biom.Table:

    rename = metadata.to_series()
    if axis == 'feature':
        axis = 'observation'
    old_ids = table.ids(axis=axis)

    new_ids = _generate_new_names(old_ids, rename, strict, False)

    updated = table.update_ids(new_ids, axis=axis, inplace=False)

    return updated
