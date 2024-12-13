# Copyright 2024 ByteDance and/or its affiliates.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


class DefaultNoneWithType(object):
    def __init__(self, dtype):
        self.dtype = dtype


class ValueMaybeNone(object):
    def __init__(self, value):
        assert value is not None
        self.dtype = type(value)
        self.value = value


class GlobalConfigValue(object):
    def __init__(self, global_key):
        self.global_key = global_key


class RequiredValue(object):
    def __init__(self, dtype):
        self.dtype = dtype


class ListValue(object):
    def __init__(self, value, dtype=None):
        if value is not None:
            self.value = value
            self.dtype = type(value[0])
        else:
            self.value = None
            self.dtype = dtype


def get_bool_value(bool_str: str):
    bool_str_lower = bool_str.lower()
    if bool_str_lower in ("false", "f", "no", "n", "0"):
        return False
    elif bool_str_lower in ("true", "t", "yes", "y", "1"):
        return True
    else:
        raise ValueError(f"Cannot interpret {bool_str} as bool")
