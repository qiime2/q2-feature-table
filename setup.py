# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-feature-table",
    # TODO stop duplicating version string
    version='0.0.6',
    packages=find_packages(),
    install_requires=['qiime >= 2.0.6', 'q2-types >= 0.0.6',
                      'biom-format >= 2.1.5, < 2.2.0', 'seaborn',
                      'scikit-bio', 'q2templates >= 0.0.6', 'numpy'],
    package_data={'q2_feature_table': ['workflows/*md'],
                  'q2_feature_table._summarize': [
                        'summarize_assets/*.html',
                        'tabulate_seqs_assets/index.html'
                  ]},
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Functionality for working with sample by feature tables.",
    license='BSD-3-Clause',
    url="http://www.qiime.org",
    entry_points={
        'qiime.plugins':
        ['q2-feature-table=q2_feature_table.plugin_setup:plugin']
    }
)
