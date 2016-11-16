# -*- coding: utf-8 -*-
r"""
Elements of the completion that come from the uncompleted ring

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

from sage.structure.element import IntegralDomainElement

class BaseElement(IntegralDomainElement):
    r"""
    Element class for elements of :class:`CompleteRing_base` which are in its
    base ring.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: K(0)
        0

    """
    def __init__(self, parent, x):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0)
            sage: isinstance(x, BaseElement)
            True
            sage: TestSuite(x).run()

        """
        IntegralDomainElement.__init__(self, parent)
        self._x = x

    def _repr_(self):
        r"""
        Return a printable representation of this element.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(1) # indirect doctest
            1
            
        """
        return repr(self._x)

    def _add_(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(1) + K(2) # indirect doctest
            3

        """
        if isinstance(other, BaseElement):
            return self.parent()(self._x + other._x)
        raise NotImplementedError

    def _sub_(self, other):
        r"""
        Return the difference of ``self`` and ``other``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(1) - K(2) # indirect doctest
            -1

        """
        if isinstance(other, BaseElement):
            return self.parent()(self._x - other._x)
        raise NotImplementedError

    def _mul_(self, other):
        r"""
        Return the product of ``self`` and ``other``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(3) * K(2) # indirect doctest
            6

        """
        if isinstance(other, BaseElement):
            return self.parent()(self._x * other._x)
        raise NotImplementedError

    def _div_(self, other):
        r"""
        Return the quotient of ``self`` and ``other``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(3) / K(2) # indirect doctest
            3/2

        """
        if isinstance(other, BaseElement):
            return self.parent()(self._x / other._x)
        raise NotImplementedError

    def _cmp_(self, other):
        r"""
        Compare ``self`` and ``other``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: cmp(K(3), K(2)) # indirect doctest
            1

        """
        if isinstance(other, BaseElement):
            return cmp(self._x, other._x)
        raise NotImplementedError

    def valuation(self):
        r"""
        Return the valuation of this element.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(1/4)
            sage: x.valuation()
            -2

        """
        return self.parent()._v(self._x)

    def reduction(self):
        r"""
        Return the reduction of this element module the element of positive
        :meth:`valuation`.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(4)
            sage: x.reduction()
            0

        """
        return self.parent()._v.reduce(self._x)
