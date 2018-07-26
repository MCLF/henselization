# -*- coding: utf-8 -*-
r"""
Benchmarks for construction and working with Henselizations.
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
from sage.all import QQ

class Rational:
    r"""
    Benchmarks creation of relatively trivial Henselizations.

    TESTS::

        sage: from henselization.benchmarks.henselizations import Rational
        sage: Rational()

    """
    def setup(self):
        r"""
        Load the henselization monkey patches before running the benchmarks.

        EXAMPLES::

            sage: from henselization.benchmarks.henselizations import Rational
            sage: Rational().setup()

        """
        import henselization
        henselization; # silence pyflakes warning about an unused import

    def time_create(self):
        r"""
        Time the creation of trivial henselizations.

        TESTS::
    
            sage: import henselization
            sage: from henselization.benchmarks.henselizations import Rational
            sage: Rational().time_create()
    
        """
        QQ.henselization(2)
