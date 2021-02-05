# ----------------------------------------------------------------------------
# Copyright (c) 2016-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
import versioneer

setup(
    name="q2-feature-table",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    package_data={'q2_feature_table': ['citations.bib'],
                  'q2_feature_table._summarize': [
                      'summarize_assets/*.html',
                      'summarize_assets/licenses/*',
                      'summarize_assets/vega/licenses/*',
                      'summarize_assets/vega/js/*',
                      'summarize_assets/vega/css/*',
                      'tabulate_seqs_assets/js/*',
                      'tabulate_seqs_assets/index.html'],
                  'q2_feature_table._core_features': [
                      'core_features_assets/index.html'],
                  'q2_feature_table._heatmap': [
                      'assets/index.html'],
                  },
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Functionality for working with sample by feature tables.",
    license='BSD-3-Clause',
    url="https://qiime2.org",
    entry_points={
        'qiime2.plugins':
        ['q2-feature-table=q2_feature_table.plugin_setup:plugin']
    },
    zip_safe=False,
)
