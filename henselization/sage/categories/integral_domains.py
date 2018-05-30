# -*- coding: utf-8 -*-
r"""
Adds a ``henselization`` method to many rings in Sage.

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

class IntegralDomains:
    class ParentMethods:
        def henselization(self, prime):
            r"""
            Return the henselization of this ring at ``prime``.

            .. NOTE:

                IntegralDomains() is probably not the right place to put this
                generic method. It could be wherever a valuation(prime) method
                is present. So, if this package ever gets into Sage, this
                method should be added manually to all the rings that implement
                valuation(prime), i.e., number fields, orders, function fields,
                (and padics, where it's trivial.)

            """
            from sage.all import Henselization
            return Henselization(self, self.valuation(prime))
