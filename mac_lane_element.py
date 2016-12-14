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

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the completion/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.path.dirname(os.getcwd()))

from sage.structure.element import IntegralDomainElement

class MacLaneElement(IntegralDomainElement):
    r"""
    Element class for elements of :class:`CompleteRing_base` which are given by
    the limit of the coefficients in a specific degree of the key polynomials
    of a :class:`MacLaneLimitValuation`.

    EXAMPLES::

        sage: from completion import *
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

            sage: from completion import *
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

            sage: from completion import *
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

            sage: from completion import *
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

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
            sage: a = (x^2 + 1).factor()[0][0][0]
            sage: a._precision()
            1

        """
        # Let G(x) be the factor that the key polynomials of
        # self._limit_valuation approximate.
        # Write v_n=[v_1(x)=mu_1, ..., v_n(phi_n)=mu_n] for an
        # approximant of self._limit_valuation.
        v_n = self._limit_valuation._approximation
        # We can write G = phi_n + sum (a_i - b_i) phi_{n-1}^i in the
        # phi_n-adic expansion with
        # sum b_i phi_{n-1}^i, sum a_i phi_{n-1}^i the phi_{n-1}-adic
        # expansions of phi_n and G respectively.
        # Since G is a key polynomial for v_n, v_n(G) = mu_n
        # = v_n(sum(a_i - b_i) phi_{n-1}^i)
        # = min v_{n-1}(a_i - b_i) + i mu_{n-1}.
        # Thus v_{n-1}(a_i - b_i) >= mu_n - i mu_{n-1}.
        if v_n._base_valuation.phi() == v_n.domain().gen():
            # When phi_{n-1}=x, the bound on a_i-b_i is a bound on the coefficient of G in degree i
            return v_n._mu - self._degree * v_n._base_valuation(v_n._base_valuation.phi())
        if v_n._base_valuation.phi().degree() == 1:
            assert v_n._base_valuation.phi().is_monic()
            # When phi_{n-1}=x+d, we can write
            # G = sum b_i (x+d)^i = sum b_i sum bin(i,j) x^j d^(i-j)
            #   = sum x^j sum bin(i,j) b_i d^(i-j)
            # Since we know a lower bound on v(a_i-b_i), we get a lower bound
            # on the error of the j-th coefficient we introduce if we replace
            # b_i in a_i in the expansion of G (i.e., if we go from G to
            # phi_n.)
            from sage.rings.all import infinity
            from sage.all import binomial
            error = infinity
            j = self._degree
            for i in range(j, v_n.phi().degree() + 1):
                error_ai_bi = v_n._mu - i * v_n._base_valuation(v_n._base_valuation.phi())
                error = min(error, v_n._base_valuation(binomial(i,j)) + error_ai_bi + (i-j) * v_n._base_valuation(v_n._base_valuation.phi()[0]))
            return error
            
        # It should be possible to generalize this to other cases, but
        # it has not been done yet.
        raise NotImplementedError

    def valuation(self):
        r"""
        Return the valuation of this element.

        EXAMPLES::

            sage: from completion import *
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

            sage: from completion import *
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

            sage: from completion import *
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
