# ----------------------------------------------------------------------------
# Copyright (c) 2016-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


class FakeCaptureHolder:
    """
    Mimic just enough of the CaptureHolder that this can be used in tests that
    call unwrapped functions that expect a CaptureHolder
    """
    def __init__(self, value=None):
        self._value = value
        self._set = False

    @property
    def value(self):
        return self._value

    def set_value(self, value):
        self._value = value
