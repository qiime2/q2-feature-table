# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import json

from qiime.sdk.type import Type
import biom
import pandas as pd

from feature_table import __version__


class FeatureTable(Type, variant_of=(Type.Artifact, Type.Metadata),
                   fields='Content'):
    class Content:
        pass

    def load(self, data_reader):
        fh = data_reader.get_file('feature-table.biom')
        return biom.Table.from_json(json.load(fh))

    def save(self, data, data_writer):
        fh = data_writer.create_file('feature-table.biom')
        fh.write(data.to_json(generated_by='feature-table %s' % __version__))

    def get_columns(self, data):
        columns = data.ids(axis='observation')
        return pd.Series([None] * len(columns), index=columns)

    def get_values(self, data, column):
        return pd.Series(data.data(column, axis='observation'),
                         index=data.ids(axis='sample'))

class Frequency(Type, variant_of=FeatureTable.Content):
    pass


class RelativeFrequency(Type, variant_of=FeatureTable.Content):
    pass


class PresenceAbsence(Type, variant_of=FeatureTable.Content):
    pass
