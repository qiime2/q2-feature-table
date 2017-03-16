# ----------------------------------------------------------------------------
# Copyright (c) 2016-2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

setup(
    name="q2-feature-table",
    version='2017.3.0.dev',
    packages=find_packages(),
    install_requires=['qiime2 == 2017.3.*', 'q2-types == 2017.3.*',
                      'q2templates == 2017.3.*', 'seaborn', 'numpy',
                      'biom-format >= 2.1.5, < 2.2.0', 'scikit-bio',
                      # `ipywidgets` included to avoid ShimWarning from
                      # `seaborn` imports:
                      #  https://github.com/mwaskom/seaborn/issues/874
                      'ipywidgets'],
    package_data={'q2_feature_table._summarize': [
                        'summarize_assets/*.html',
                        'summarize_assets/app/*.js',
                        'summarize_assets/app/js/*.js',
                        'summarize_assets/app/css/*.css',
                        'tabulate_seqs_assets/js/*',
                        'tabulate_seqs_assets/index.html'
                  ]},
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
