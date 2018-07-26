# -*- coding: utf-8 -*-
r"""
Elements in a Henselization which are described by a limit valuation
====================================================================

AUTHORS:

- Julian Rüth (2016-11-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016-2018 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from .henselization_element import HenselizationElement_base, HenselizationElement_Ring, HenselizationElement_Field

class MacLaneElement_base(HenselizationElement_base):
    r"""
    Element class for elements of :class:`CompleteRing_base` which are given by
    the limit of the coefficients in a specific degree of the key polynomials
    of a :class:`MacLaneLimitValuation`.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(5)
        sage: R.<x> = K[]
        sage: F = (x^2 + 1).factor()

    The coefficients of the factors are the limits of the coefficients of the
    key polynomials::

        sage: F
        (x + 2 + O(5)) * (x + 3 + O(5))

    """
    def __init__(self, parent, limit_valuation, degree):
        r"""
        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: F = (x^2 + 1).factor()
            sage: a = F[0][0][0]
            sage: from sage.rings.padics.henselization.mac_lane_element import MacLaneElement_base
            sage: isinstance(a, MacLaneElement_base)
            True
            sage: TestSuite(a).run() # long time

        """
        super(MacLaneElement_base, self).__init__(parent)
        self._limit_valuation = limit_valuation
        self._degree = degree

    def _repr_(self):
        r"""
        Return a printable representation of this element.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: (x^2 + 1).factor() # indirect doctest
            (x + 2 + O(5)) * (x + 3 + O(5))

        """
        approximation = self._limit_valuation._initial_approximation.phi()[self._degree]
        if approximation:
            approximation = repr(approximation)
        precision = self._precision()
        uniformizer = self._limit_valuation.uniformizer()
        error = repr(uniformizer)
        if precision != 1:
            if any([op in error for op in '+-*/']):
                error = "(%s)"%(error,)
            power = repr(precision)
            if any([op in power for op in '+-*/']):
                power = "(%s)"%(power,)
            error = "%s^%s"%(error, power)
        error = "O(%s)"%(error,)
        if approximation:
            return "%s + %s"%(approximation, error)
        else:
            return error

    def _cache_key(self):
        r"""
        Return a key which identifies this element for caching.

        EXAMPLES:

        Mac Lane elements are not hashable as it seems too hard to implement a
        hash function that is fast and consistent with equality::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: F = (x^2 + 1).factor()
            sage: a = F[0][0][0]
            sage: hash(a)
            Traceback (most recent call last):
            ...
            TypeError: ... is not hashable

        Therefore, they implement :meth:`_cache_key`, so that they can be used
        in caches::

            sage: a._cache_key()
            ([ Gauss valuation induced by 5-adic valuation, v(x + 2) = 1 , … ], 0)
            
        """
        return self._limit_valuation, self._degree

    def _richcmp_(self, other, op):
        r"""
        Return whether this element relates to ``other`` with respect to
        ``op``.

        EXAMPLES:

        We can sometimes compare elements that came from the same factorization::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: F = (x^2 + 1).factor()
            sage: f = F[0][0] # x + 2 + …
            sage: g = F[1][0] # x - 2 - …
            sage: f[0] == g[0]
            False
            sage: f[0] == f[0]
            True

        In some cases, we can also compare to exact elements::

            sage: G = GaussianIntegers().fraction_field()
            sage: K = G.henselization(2)
            sage: R.<x> = K[]
            sage: F = (x^2 + 1).factor()
            sage: f = F[0][0] # x + I
            sage: g = F[1][0] # x - I
            sage: f[0] == g[0]
            False
            sage: f[0] == f[0]
            True
            sage: f[0] == 0
            False

        """
        if op == 2:
            from base_element import BaseElement_base
            if isinstance(other, BaseElement_base):
                if other == 0 and self._degree == 0:
                    # if the constant coefficient were zero, the polynomial would be reducible
                    return False
                if (other - self._limit_valuation._approximation.phi()[self._degree]).valuation() < self._precision():
                    return False
                if self._degree == 0:
                    phi = self._limit_valuation._approximation.phi()
                    if phi.degree() == 1:
                        # other is the constant of the defining polynomial of this
                        # element if x + other has infinite valuation
                        x = phi.parent().gen()
                        from sage.rings.all import infinity
                        return self._limit_valuation(x + other) is infinity

                # We could try to push the approximation indefinitely (but this won't work if other is actually equal)
                raise NotImplementedError("comparison of %s which is the coefficient in degree %s of %s to base element %s"%(self, self._degree, self._limit_valuation, other))
            if isinstance(other, MacLaneElement_base):
                if self._limit_valuation.parent() is other._limit_valuation.parent():
                    if self._limit_valuation._G == other._limit_valuation._G:
                        if self._limit_valuation == other._limit_valuation:
                            if self._degree == other._degree:
                                # the same coefficient of the same factor
                                return True
                            else:
                                raise NotImplementedError("comparison of coefficients of an approximate factor")
                        else:
                            if self._degree == other._degree and self._limit_valuation._approximation.phi().degree() == other._limit_valuation._approximation.phi().degree() == 1:
                                # the constant coefficients of two linear
                                # factors of the same polynomial can not be
                                # identical if their factors are different
                                return False
                            else:
                                raise NotImplementedError("comparison of coefficients of different factors")
                    else:
                        raise NotImplementedError("comparison of Mac Lane elements that come from different factorizations")
                raise NotImplementedError("comparison of Mac Lane elements that come from valuations on different rings")
        elif op == 3:
            return not (self == other)
        raise NotImplementedError

    def _precision(self):
        r"""
        Return the precision (in terms of valuation) to which the element is
        known.

        ALGORITHM:

        Suppose that `g` is an (exact) factor of the polynomial `G` and let
        `\phi` be `v_n = [v_0,\dots,v_n(\phi_n)=\mu_n]` be an approximant of
        the valuation corresponding to `g`. Suppose that `\theta` is a root of
        `g` and write `\Delta = \phi-g`.
        To measure the precision of the factors of `\phi`, we are
        interested in `v_0(\Delta)` where `v_0` denotes the Gauss valuation.
        We have `v_0(\Delta) > v(\Delta(\theta)) - v(g_{\deg\Delta}(\theta)) =
        v(\phi(\theta)) - v(\prod \phi_i^{j_i})=\mu_n - \sum j_i \mu_i` where
        `g_{\deg\Delta}` is the best approximation to `g` of degree
        `\deg\Delta`, i.e., a monic polynomial with maximal valuation at
        `\theta`. (cf. Lemma 4.5 in [GNP2012])

        REFERENCES:

        .. [GNP2012] Jordi Guàrdia, Enric Nart, Sebastian Pauli
                     "Single-factor lifting and factorization of polynomials over local fields"
                     Journal of Symbolic Computation 47 (2012) 1318-1346

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: a._precision()
            1

        """
        w = self._limit_valuation._approximation
        e = reversed([v.E()/v._base_valuation.E() for v in w.augmentation_chain()[:-1]])
        h = reversed([v.mu() - v._base_valuation(v.phi()) for v in w.augmentation_chain()[:-1]])

        v0 = 0
        e0 = 1
        for hi,ei in zip(h,e):
            e0 *= ei
            v0 += hi/e0

        valuations = w.valuations(self._limit_valuation._G)
        h_phi = valuations.next() - valuations.next()
        return v0 + h_phi / e0

    def valuation(self):
        r"""
        Return the valuation of this element.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: a.valuation()
            0

        """
        ret = self._limit_valuation._approximation.phi()[self._degree].valuation()
        if self._precision() > ret:
            return ret
        # we could try to push the approximation indefinitely (but this won't work if this element is actually zero)
        raise NotImplementedError

    def reduction(self):
        r"""
        Return the reduction of this element module the element of positive
        :meth:`valuation`.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: a.reduction()
            2

        """
        if self._precision() > 0:
            return self._limit_valuation._approximation.phi()[self._degree].reduction()
        # we could try to push the approximation indefinitely (and it would actually work)
        raise NotImplementedError

    def approximation(self, precision):
        r"""
        Return an approximation to this element which is known to differ from
        the actual by at most ``precision``.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: a
            2 + O(5)
            sage: a.approximation(precision=10)
            6139557

        """
        while self._precision() < precision:
            self._limit_valuation._improve_approximation()
        return self.parent()(self._limit_valuation._approximation.phi()[self._degree])

class MacLaneElement_Ring(MacLaneElement_base, HenselizationElement_Ring):
    r"""
    A :class:`MacLaneElement_base` that lives in a ring that is not a field.

    EXAMPLES::

        sage: from henselization import *
        sage: S = ZZ.henselization(5)
        sage: R.<x> = S[]
        sage: F = (x^2 + 1).factor()
        sage: a = F[0][0][0]; a
        2 + O(5)

    TESTS::

        sage: from sage.rings.padics.henselization.mac_lane_element import MacLaneElement_Ring
        sage: isinstance(a, MacLaneElement_Ring)
        True
        sage: TestSuite(a).run() # long time

    """

class MacLaneElement_Field(MacLaneElement_base, HenselizationElement_Field):
    r"""
    A :class:`MacLaneElement_base` that lives in a field.

    EXAMPLES::

        sage: from henselization import *
        sage: R.<x> = QQ.henselization(5)[]
        sage: F = (x^2 + 1).factor()
        sage: a = F[0][0][0]; a
        2 + O(5^10)

    TESTS::

        sage: from sage.rings.padics.henselization.mac_lane_element import MacLaneElement_Field
        sage: isinstance(a, MacLaneElement_Field)
        True
        sage: TestSuite(a).run() # long time

    """
