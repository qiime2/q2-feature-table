# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

from feature_table import __version__

setup(
    name="fetaure_table",
    version=__version__,
    packages=find_packages(),
    install_requires=['biom-format >= 2.1.5, < 2.2.0', 'scipy', 'IPython',
                      'ipywidgets', 'seaborn'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Functionality for working with sample x feature tables.",
    license="BSD",
    url="http://www.qiime.org",
)
