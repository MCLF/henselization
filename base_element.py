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
from sage.misc.cachefunc import cached_method

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

    @cached_method
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

    def _vector_(self, base=None):
        r"""
        Return the coefficients of this element over the power basis this
        elements parents over ``base``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K(1)._vector_()
            (1,)

        """
        if base is None:
            base = self.parent()
        if base is self.parent():
            return (self,)
        else:
            return sum((self.parent().base_ring()(c)._vector_(base) for c in list(self._x)), ())

    def matrix(self, base=None):
        r"""
        Return the matrix of this element over ``base``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.gen().matrix(base=K)
            [ 0  1]
            [-1 -1]

        """
        V,f,g = self.parent().vector_space(base=base)
        return V.hom([g(f(b)*self) for b in V.basis()]).matrix()

    def approximation(self, precision=None):
        r"""
        Return an approximation to this element which is know to differ from
        the actual element by at most ``precision``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.gen().approximation(123)
            x

        """
        return self
        

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
                other_unit_part = self.parent()(self.parent()._base_fraction_field_valuation.shift(other._x, -other.valuation()))
                x = (self / other_unit_part)._x
                if x in self.parent()._base_valuation.domain():
                    x = self.parent()._base_valuation.domain()(x)
                    ret = self.parent()._base_valuation.shift(x, -other.valuation())
                else:
                    raise NotImplementedError("shift of element with integral denominator")
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
