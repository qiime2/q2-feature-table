# ----------------------------------------------------------------------------
# Copyright (c) 2016--, QIIME development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

from normalize import __version__

setup(
    name="normalize",
    version=__version__,
    packages=find_packages(),
    install_requires=['biom-format >= 2.1.5, < 2.2.0'],
    author="Greg Caporaso",
    author_email="gregcaporaso@gmail.com",
    description="Normalization of sample by feature tables.",
    license="BSD",
    keywords="biom normalize ",
    url="http://www.qiime.org",
)
