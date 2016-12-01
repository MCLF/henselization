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

from sage.structure.element import Element, FieldElement, IntegralDomainElement, coerce_binop

class BaseElement_base(IntegralDomainElement):
    r"""
    Abstract base class for elements of :class:`Completion_base` which are in
    one of its uncompleted ``base`` fields.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: x = K(0); x
        0

    """
    def __init__(self, parent, base, valuation, x):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0)
            sage: isinstance(x, BaseElement_base)
            True
            sage: TestSuite(x).run()

        """
        super(BaseElement_base, self).__init__(parent)
        self._x = x
        self._base = base
        self._valuation = valuation

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
        if isinstance(other, BaseElement_base):
            if self._base is other._base and self._valuation is other._valuation:
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
        if isinstance(other, BaseElement_base):
            if self._base is other._base and self._valuation is other._valuation:
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
        if isinstance(other, BaseElement_base):
            if self._base is other._base and self._valuation is other._valuation:
                return self.parent()(self._x * other._x)
        raise NotImplementedError

    def _div_(self, other):
        r"""
        Return the quotient of ``self`` and ``other``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: K = Completion(ZZ, v)
            sage: K(2) / K(3) # indirect doctest
            2/3

        """
        if isinstance(other, BaseElement_base):
            if self._base is other._base and self._valuation is other._valuation:
                return self.parent()(self._x / other._x)
        raise NotImplementedError

    def _richcmp_(self, other, op):
        r"""
        Return the result of comparing this element to ``other`` with respect
        to ``op``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(3) == K(2) # indirect doctest
            False

        """
        if op == 2: # ==
            if isinstance(other, BaseElement_base):
                if self._base is other._base and self._valuation is other._valuation:
                    return  self._x == other._x
        if op == 3: # !=
            return not (self == other)
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
        return self._valuation(self._x)

    def reduction(self):
        r"""
        Return the reduction of this element modulo the elements of positive
        :meth:`valuation`.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(4)
            sage: x.reduction()
            0

        """
        return self._valuation.reduce(self._x)


class BaseElement_Ring(BaseElement_base):
    r"""
    An element of :class:`Completion_Ring` which is in one of its uncompleted
    ``base`` fields.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: R = Completion(ZZ, v)
        sage: x = R(0); x
        0

    TESTS::

        sage: isinstance(R(0), BaseElement_Ring)
        True

    """
    def _floordiv_(self, other):
        r"""
        Return the quotient of ``self`` and ``other``.

        We implement the variation (3) as given in
        :meth:`sage.rings.padics.pAdicGenericElement._mod_`.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: R(3) // R(2) # indirect doctest
            1

        """
        if isinstance(other, BaseElement_Ring):
            if self._base is other._base and self._valuation is other._valuation:
                from sage.rings.all import ZZ
                other_unit_part = self.parent()(self.parent()._base_valuation.shift(other._x, -other.valuation()))
                ret = self.parent()._base_valuation.shift((self / other_unit_part)._x, - other.valuation())
                return self.parent()(ret)
        raise NotImplementedError

    def _mod_(self, other):
        r"""
        Return the remainder of division of this element by ``other``.

        We implement the variation (3) as given in
        :meth:`sage.rings.padics.pAdicGenericElement._mod_`.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: R(3) % R(2) # indirect doctest
            1

        """
        return self - (self // other) * other


class BaseElement_Field(BaseElement_base, FieldElement):
    r"""
    An element of :class:`Completion_Field` which is in one of its uncompleted
    ``base`` fields.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: R = Completion(ZZ, v)
        sage: R(0)
        0

    """
    def __init__(self, parent, base, valuation, x):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: isinstance(K(0), BaseElement_Field)
            True

        """
        super(BaseElement_Field, self).__init__(parent, base, valuation, x)
