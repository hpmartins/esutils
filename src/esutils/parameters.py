# coding: utf-8
# /*##########################################################################
#
# Copyright (c) 2016-2018 European Synchrotron Radiation Facility
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/

from __future__ import absolute_import, division, unicode_literals

import copy
import datetime
import gzip
import json
import numpy as np
import os
try:
    import cPickle as pickle
except ImportError:
    import pickle
import re
import subprocess
import sys

import pkg_resources

import collections
class odict(collections.OrderedDict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

class SlaterIntegrals(object):
    path = pkg_resources.resource_filename('esutils', 'data/slater_integrals.json')
    with open(path, 'r') as rf:
        table = json.load(rf)
        
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        


class QuantyCalculation(object):
    # Make the parameters a class attribute. This speeds up the creation
    # of a new calculation object; significantly.
    path = pkg_resources.resource_filename('esutils', 'data/parameters.json.gz')

    with gzip.open(path, 'rb') as p:
        parameters = json.loads(p.read().decode('utf-8'), object_pairs_hook=odict)

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

        parameters = self.parameters