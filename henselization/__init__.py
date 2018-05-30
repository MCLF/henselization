# -*- coding: utf-8 -*-
r"""
Monkey patching to enable Henselizations in the Sage library.

EXAMPLES::

    sage: from henselization import *
    sage: ZZ.henselization(2)
    Henselization of Integer Ring with respect to 2-adic valuation

"""
#*****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import
from . import sage as monkey_sage
import sage.all

from recursive_monkey_patch import monkey_patch
monkey_patch(monkey_sage, sage)

from sage.rings.padics.henselization.henselization import Henselization
sage.all.Henselization = Henselization
