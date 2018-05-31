# -*- coding: utf-8 -*-

#*****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .util import AbstractMonkey

class Monkey(AbstractMonkey):
    _trac = "https://trac.sagemath.org/ticket/22983"

    def _test(self):
        from sage.all import QQ, PolynomialRing
        R = PolynomialRing(QQ, 'x')
        x = R.gen()
        if R.quo(x) is not R.quo(x):
            raise Exception("#22983 has not been fixed")
    
    def _patch(self):
        import sage.rings.polynomial.polynomial_quotient_ring
        sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRingFactory = PolynomialQuotientRingFactory
        sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing = PolynomialQuotientRingFactory("sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing")
        del sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_generic.__reduce__
        del sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_domain.__reduce__
        del sage.rings.polynomial.polynomial_quotient_ring.PolynomialQuotientRing_field.__reduce__

from sage.structure.factory import UniqueFactory
class PolynomialQuotientRingFactory(UniqueFactory):
    def create_key(self, ring, polynomial, names=None):
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_commutative
        if not isinstance(ring, PolynomialRing_commutative):
            raise TypeError("ring must be a polynomial ring")
        from sage.rings.polynomial.polynomial_element import Polynomial
        if not isinstance(polynomial, Polynomial):
            raise TypeError("must be a polynomial")
        if not polynomial.parent() is ring:
            raise TypeError("polynomial must be in ring")

        c = polynomial.leading_coefficient()
        if not c.is_unit():
            raise TypeError("polynomial must have unit leading coefficient")

        if names is None:
            names = tuple([x + 'bar' for x in ring.variable_names()])
        else:
            from sage.rings.polynomial.polynomial_quotient_ring import normalize_names
            names = normalize_names(ring.ngens(), names)

        return ring, polynomial, names

    def create_object(self, version, key):
        ring, polynomial, names = key

        R = ring.base_ring()
        from sage.categories.all import IntegralDomains
        if R in IntegralDomains():
            try:
                is_irreducible = polynomial.is_irreducible()
            except NotImplementedError: # is_irreducible sometimes not implemented
                pass
            else:
                if is_irreducible:
                    from sage.categories.all import Fields
                    if R in Fields():
                        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_field
                        return PolynomialQuotientRing_field(ring, polynomial, names)
                    else:
                        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_domain
                        return PolynomialQuotientRing_domain(ring, polynomial, names)
        from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
        return PolynomialQuotientRing_generic(ring, polynomial, names)

Monkey().patch()
