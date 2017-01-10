# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pkg_resources
import shutil

import pandas as pd
import qiime2
import q2templates

TEMPLATES = pkg_resources.resource_filename('q2_feature_table',
                                            '_explore_rarefaction')


def explore_rarefaction(output_dir: str, table: pd.DataFrame,
                        sample_metadata: qiime2.Metadata=None) -> None:

    with open(os.path.join(output_dir, 'data.jsonp'), 'w') as fh:
        fh.write("load(")
        table.to_json(fh)
        fh.write(', ')
        if sample_metadata:
            sample_metadata.to_dataframe().to_json(fh)
        else:
            fh.write('{}')
        fh.write(', ')
        table.sum(axis=1).to_json(fh)
        fh.write(');')

    context = {'min_count': 0,
               'max_count': table.sum(axis=1).max()}

    index = os.path.join(TEMPLATES, 'assets', 'index.html')
    q2templates.render(index, output_dir, context=context)

    shutil.copytree(os.path.join(TEMPLATES, 'assets', 'app'),
                    os.path.join(output_dir, 'app'))
