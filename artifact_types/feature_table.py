# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime.sdk.type import Type
import biom
from feature_table import __version__


class FeatureTable(Type, variant_of=Type.Artifact, fields='Content'):
    class Content:
        pass

    def load(self, data_reader):
        fh = data_reader.get_file('feature-table.biom')
        return Table.from_json(json.load(fh))

    def save(self, data, data_writer):
        fh = data_writer.create_file('feature-table.biom')
        fh.write(data.to_json(generated_by='feature-table %s' % __version__))


class Frequency(Type, variant_of=FeatureTable.Content):
    pass


class RelativeFrequency(Type, variant_of=FeatureTable.Content):
    pass


class PresenceAbsence(Type, variant_of=FeatureTable.Content):
    pass
