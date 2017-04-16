# -*- coding: utf-8 -*-
r"""
Elements in a completion which are described by a limit valuation.

AUTHORS:

- Julian Rüth (2016-11-16): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.element import IntegralDomainElement

class MacLaneElement(IntegralDomainElement):
    r"""
    Element class for elements of :class:`CompleteRing_base` which are given by
    the limit of the coefficients in a specific degree of the key polynomials
    of a :class:`MacLaneLimitValuation`.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: C = Completion(QQ, v)
        sage: R.<x> = C[]
        sage: F = (x^2 + 1).factor()

    The coefficients of the factors are the limits of the coefficients of the
    key polynomials::

        sage: F
        (x + 2 + O(?)) * (x + 3 + O(?))

    """
    def __init__(self, parent, limit_valuation, degree):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
            sage: F = (x^2 + 1).factor()
            sage: a = F[0][0][0]
            sage: isinstance(a, MacLaneElement)
            True
            sage: TestSuite(a).run() # long time

        """
        super(MacLaneElement, self).__init__(parent)
        self._limit_valuation = limit_valuation
        self._degree = degree

    def _repr_(self):
        r"""
        Return a printable representation of this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
            sage: (x^2 + 1).factor() # indirect doctest
            (x + 2 + O(?)) * (x + 3 + O(?))

        """
        approximation = repr(self._limit_valuation._initial_approximation.phi()[self._degree])
        if ' ' in approximation:
            approximation = "(" + approximation + ")"
        return approximation + " + O(?)"

    def _richcmp_(self, other, op):
        r"""
        Return whether this element relates to ``other`` with respect to
        ``op``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: b = (x^2 + 1).factor()[0][0][0]
            sage: a == b
            True

        """
        if op == 2:
            from base_element import BaseElement_base
            if isinstance(other, BaseElement_base):
                if (other - self._limit_valuation._approximation.phi()[self._degree]).valuation() < self._precision():
                    return False
                # we could try to push the approximation indefinitely (but this won't work if other is actually equal)
                raise NotImplementedError("comparison to base elements")
            if isinstance(other, MacLaneElement):
                if self._limit_valuation.parent() is other._limit_valuation.parent():
                    # we currently only handle the trivial case here, i.e., the
                    # elements are indistinguishable
                    return (self._limit_valuation == other._limit_valuation
                            and self._degree == other._degree)
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: a
            2 + O(?)
            sage: a.approximation(precision=10)
            6139557

        """
        while self._precision() < precision:
            self._limit_valuation._improve_approximation()
        return self.parent()(self._limit_valuation._approximation.phi()[self._degree])
