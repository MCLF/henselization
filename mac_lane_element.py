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
        (x + 57 + O(?)) * (x + 68 + O(?))

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
        IntegralDomainElement.__init__(self, parent)
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
            (x + 57 + O(?)) * (x + 68 + O(?))

        """
        approximation = repr(self._limit_valuation._initial_approximation.phi()[self._degree])
        if ' ' in approximation:
            approximation = "(" + approximation + ")"
        return approximation + " + O(?)"

    def _cmp_(self, other):
        r"""
        Compare ``self`` and ``other``.

        """
        raise NotImplementedError

    def _richcmp(self, other, op):
        r"""
        """
        print self,other,op
        raise NotImplementedError

    def valuation(self):
        r"""
        Return the valuation of this element.

        """
        raise NotImplementedError

    def reduction(self):
        r"""
        Return the reduction of this element module the element of positive
        :meth:`valuation`.

        """
        raise NotImplementedError
