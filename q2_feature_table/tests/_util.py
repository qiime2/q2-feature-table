# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.type import CaptureHolder


def fake_CaptureHolder_factory(value: any):
    """
    Creates a bad and lame CaptureHolder purely to replace passing a raw value
    into raw unwrapped functions when testing.

    Paramters
    ---------
    value : any
        Whatever value you are trying to pass into your function that is
        expecting a CaptureHolder
    """
    return CaptureHolder('', value, None, None)
