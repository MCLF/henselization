# -*- coding: utf-8 -*-
r"""
Exact completions of rings

AUTHORS:

- Julian Rüth (2016-11-15): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the completion/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.path.dirname(os.getcwd()))

from sage.rings.ring import Ring, IntegralDomain, Field
from sage.structure.factory import UniqueFactory
from sage.misc.lazy_attribute import lazy_attribute

class CompletionFactory(UniqueFactory):
    r"""
    The completion of ``R`` with respect to ``v``, a discrete non-trivial
    valuation on ``R``.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: Completion(QQ, v)
        Completion of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, R, v):
        r"""
        Create a key which uniquely identifies this completion.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: Completion(QQ, v) is Completion(QQ, v) # indirect doctest
            True

        """
        from sage.categories.all import IntegralDomains
        if R not in IntegralDomains():
            raise ValueError("R must be a ring")
        from mac_lane import DiscretePseudoValuationSpace
        if v not in DiscretePseudoValuationSpace(R) or not v.is_discrete_valuation():
            raise ValueError("v must be a discrete valuation on R")
        if v.is_trivial():
            raise ValueError("v must not be trivial")

        return R, v

    def create_object(self, version, key):
        r"""
        Create the completion identified by ``key``.

        TESTS::
            
            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: Completion(QQ, v) # indirect doctest
            Completion of Rational Field with respect to 5-adic valuation

        """
        R, v = key
        from sage.categories.all import Fields
        if R in Fields():
            return CompleteField(R, v)
        else:
            return CompleteDomain(R, v)

Completion = CompletionFactory("Completion")

class CompleteRing_base(Ring):
    r"""
    Abstract base class for the completion of ``R`` with respect to ``v``.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: Completion(ZZ, v)
        Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, R, v, category = None):
        r"""
        TESTS::

            sage: from completion import *
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x)
            sage: C = Completion(K, v)
            sage: isinstance(C, CompleteRing_base)
            True
            sage: TestSuite(C).run() # long time
    
        """
        # TODO: This should be R.category().Completion() rather, however, the
        # CompleteDiscreteValuationRing/Field() category is currently not set
        # up correctly. First, precision_absolute() and precision_relative()
        # are implementation details. Then CDVFields() should have CDVRings()
        # as a super category.
        from sage.categories.all import Fields, IntegralDomains
        category = category or (Fields() if R in Fields() else IntegralDomains())
        Ring.__init__(self, R, category=category)

        self._v = v
        from sage.rings.all import PolynomialRing
        self._polynomial_ring = PolynomialRing(self, 'x')

    def characteristic(self):
        r"""
        Return the characteristic of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.characteristic()
            0

        """
        return self.base().characteristic()

    def is_finite(self):
        r"""
        Return whether this ring is finite.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.is_finite()
            False

        """
        # since the underlying valuation is non-trivial, there must be
        # infinitely many elements
        return False

    def _coerce_map_from_(self, other):
        r"""
        Return a coercion from ``other`` to this ring if one exists.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.has_coerce_map_from(QQ) # indirect doctest
            True

        """
        if self.base().has_coerce_map_from(other):
            return True

    def _an_element_(self):
        r"""
        Return an element of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.an_element()
            0

        """
        return self(0)

    def _element_constructor_(self, x):
        r"""
        Create an element in this parent from ``x``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0) # indirect doctest
            sage: x.parent() is K
            True

        """
        if x in self.base():
            x = self.base()(x)
            return self._base_element_class(self, x)
        if isinstance(x, tuple):
            if len(x) == 2:
                v, i = x
                from sage.rings.all import NN
                from mac_lane.limit_valuation import MacLaneLimitValuation
                if isinstance(v, MacLaneLimitValuation) and v.domain() is self._polynomial_ring and i in NN:
                    return self._mac_lane_element_class(self, v, i)
        raise ValueError("Can not convert %r to an element in %r"%(x, self))

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the base ring.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0) # indirect doctest
            sage: isinstance(x, K._base_element_class)
            True

        """
        from base_element import BaseElement
        return self.__make_element_class__(BaseElement)

    @lazy_attribute
    def _mac_lane_element_class(self):
        r"""
        The class for elements which are given as coefficients of key
        polynomial of limit valuations.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: F = f.factor() # indirect doctest
            sage: isinstance(F[0][0][0], K._mac_lane_element_class)
            True

        """
        from mac_lane_element import MacLaneElement
        return self.__make_element_class__(MacLaneElement)

    def _repr_(self):
        r"""
        Return a printable representation of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: Completion(QQ, v)
            Completion of Rational Field with respect to 2-adic valuation

        """
        return "Completion of %r with respect to %r"%(self.base_ring(), self._v)

    def some_elements(self):
        r"""
        Return some typical elements of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: len(K.some_elements())
            100

        """
        return map(self, self.base_ring().some_elements())

    def valuation(self):
        r"""
        Return the valuation on this complete ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation()
            2-adic valuation

        """
        from .valuation import Valuation
        return Valuation(self)

    def _factor_univariate_polynomial(self, f):
        r"""
        Return the factorization of ``f`` over this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x + 1
            sage: f.factor() # indirect doctest
            x + 1
            sage: f = x^2 + 1
            sage: f.factor()
            x^2 + 1

        A non-trivial example::

            sage: G = GaussianIntegers().fraction_field()
            sage: v = pAdicValuation(G, 2)
            sage: K = Completion(G, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: f.factor() # long time
            (x + (-I - 8) + O(?)) * (x + (5*I - 4) + O(?))

        Another non-trivial example::

            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: f.factor() # long time
            (x + 57 + O(?)) * (x + 68 + O(?))

        """
        if f.is_constant():
            raise NotImplementedError
        if not f.is_monic():
            raise NotImplementedError
        if not f.is_squarefree():
            # TODO: piece together the factorization of the factors of the squarefree decomposition
            raise NotImplementedError

        from sage.structure.factorization import Factorization
        G = self._polynomial_ring(f)
        approximants = self.valuation().mac_lane_approximants(G)
        if len(approximants) == 1:
            return Factorization([(G, 1)])
        approximants = sum([approximant.mac_lane_step(G) for approximant in approximants], [])
        factors = []
        for approximant in approximants:
            # we only need to perform a Mac Lane step when there is ramification
            # introduced in the last step (to get the degrees right), anyway, it
            # does not hurt to do another step
            approximant = approximant.mac_lane_step(G)
            assert len(approximant)==1
            approximant = approximant[0]

            from sage.rings.all import infinity
            if approximant(approximant.phi()) == infinity:
                factors.append(approximant.phi())
                continue

            degree = approximant.phi().degree()
            from mac_lane.limit_valuation import LimitValuation
            limit = LimitValuation(approximant, G)
            coefficients = [self((limit, d)) for d in range(degree + 1)]
            if f.is_monic():
                coefficients[-1] = self(1)
            factor = f.parent()(coefficients)
            factors.append(factor)
        return Factorization([(factor, 1) for factor in factors], unit=self.one(), sort=False)

class CompleteDomain(CompleteRing_base, IntegralDomain):
    r"""
    Base class for the completion of the integral domain ``R`` with respect
    to ``v``.

    EXAMPLES::

        sage: from completion import *
        sage: R.<x> = QQ[]
        sage: v = pAdicValuation(QQ, 2)
        sage: w = GaussValuation(R, v)
        sage: Completion(R, w)
        Completion of Univariate Polynomial Ring in x over Rational Field with respect to Gauss valuation induced by 2-adic valuation

    """
    def __init__(self, R, v, category = None):
        r"""
        TESTS::

            sage: from completion import *
            sage: R.<x> = QQ[]
            sage: v = pAdicValuation(QQ, 2)
            sage: w = GaussValuation(R, v)
            sage: C = Completion(R, w)
            sage: isinstance(C, CompleteDomain)
            True
            sage: TestSuite(C).run()
            
        """
        CompleteRing_base.__init__(self, R, v, category)
        IntegralDomain.__init__(self, R, category=self.category())

class CompleteField(CompleteRing_base, Field):
    r"""
    Abstract base class for the completion of the field ``K`` with respect to
    ``v``.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: Completion(QQ, v)
        Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, K, v, category = None):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: isinstance(K, CompleteField)
            True
            sage: TestSuite(K).run()

        """
        CompleteRing_base.__init__(self, K, v, category)
        Field.__init__(self, K, category=self.category())
