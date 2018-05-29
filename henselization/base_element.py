# -*- coding: utf-8 -*-
r"""
Elements of the Henselization that come from the uncompleted ring

AUTHORS:

- Julian Rüth (2016-11-15): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016-2018 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.cachefunc import cached_method
from henselization_element import HenselizationElement_Field, HenselizationElement_Ring, HenselizationElement_base
from sage.all import ZZ, QQ

class BaseElement_base(HenselizationElement_base):
    r"""
    Abstract base class for elements of :class:`Henselization_base` which are in
    one of its uncompleted ``base`` fields.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: x = K(0); x
        0

    """
    def __init__(self, parent, base, valuation, x):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: x = K(0)
            sage: isinstance(x, BaseElement_base)
            True
            sage: TestSuite(x).run()

        """
        super(BaseElement_base, self).__init__(parent)
        self._x = x
        self._base = base
        self._valuation = valuation

    def _integer_(self):
        return ZZ(self._x)

    def _rational_(self):
        return QQ(self._x)

    def __hash__(self):
        r"""
        Return a hash value for this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: hash(K(1)) # indirect doctest
            1

        """
        return hash(self._x)

    def _repr_(self):
        r"""
        Return a printable representation of this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(1) # indirect doctest
            1
            
        """
        return repr(self._x)

    def _add_(self, other):
        r"""
        Return the sum of ``self`` and ``other``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = ZZ.valuation(2)
            sage: K = Henselization(ZZ, v)
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(3) == K(2) # indirect doctest
            False

        """
        if op == 2: # ==
            from mac_lane_element import MacLaneElement_base
            if isinstance(other, MacLaneElement_base):
                return other._richcmp_(self, op)
            elif isinstance(other, BaseElement_base):
                from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
                if (self._base is other._base or
                    # polynomial quotient rings are not unique parents yet so
                    # we need to work around the failing "is"
                    (is_PolynomialQuotientRing(self._base) and is_PolynomialQuotientRing(other._base) and self._base == other._base)):
                    if self._valuation is other._valuation:
                        return self._x == other._x
            elif other in self.parent().base():
                return self._x == other
        if op == 3: # !=
            return not (self == other)
        raise NotImplementedError

    @cached_method
    def valuation(self):
        r"""
        Return the valuation of this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: x = K(4)
            sage: x.reduction()
            0

        """
        return self._valuation.reduce(self._x)

    def _vector_(self, base=None):
        r"""
        Return the coefficients of this element when written as a linear
        combination over a basis of its parent over ``base``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(1)._vector_()
            (1,)

        In a simple extension, the coefficients are taken with respect to the
        internal number field representation of the elements::

            sage: R.<x> = K[]
            sage: L.<x> = K.extension(x^2 - 2)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 - x)
            sage: y = M(M.base().gen())
            sage: y._vector_(K)
            (0, 1, 0, 0)
            sage: (y^4)._vector_(K)
            (-62, 0, 0, 0)

        """
        if base is None:
            base = self.parent()
        if base is self.parent():
            return (self,)
        elif self.parent().base().base_ring() is base.base():
            return tuple(self._x.list())
        elif self.parent().base_ring().base() is self.parent().base().base_ring():
            return sum((self.parent().base_ring()(c)._vector_(base) for c in self._x.list()), ())
        else:
            raise NotImplementedError("Vector space representation of an element of %s over %s"%(self.parent(), base))

    def matrix(self, base=None):
        r"""
        Return the matrix of this element over ``base``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.gen().matrix(base=K)
            [ 0  1]
            [-1 -1]

        """
        V,f,g = self.parent().module(base=base)
        return V.hom([g(f(b)*self) for b in V.basis()]).matrix()

    def approximation(self, precision=None):
        r"""
        Return an approximation to this element which is know to differ from
        the actual element by at most ``precision``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.gen().approximation(123)
            x

        """
        return self

    def _relative_size(self):
        r"""
        Return an estimate on the coefficient size of this element.

        The number returned is an estimate on the factor between the number of
        Bits used by this element and the minimal number of bits used by an
        element congruent to it.

        This is used by :meth:`simplify` to decide whether simplification of
        coefficients is going to lead to a significant shrinking of the
        coefficients of this elements.

        EXAMPLES:: 

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(1024)._relative_size()
            6

        """
        return self._valuation._relative_size(self._x)

    def simplify(self, error=None, force=False):
        r"""
        Return a simplified version of this element.

        Produce an element which differs from this element by an element of
        valuation strictly greater than the valuation of this element  (or
        strictly greater than ``error`` if set.)

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(1025).simplify(force=True)
            1

        """
        return self.parent()(self._valuation.simplify(self._x, error=error, force=force))

    def _upper_bound(self):
        r"""
        Return an upper bound of the valuation of this element.

        Use this method to get an approximation of the valuation when speed is
        more important than accuracy.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(1025)._upper_bound()
            0

        """
        return self._valuation.upper_bound(self._x)

    def _lower_bound(self):
        r"""
        Return an lower bound of the valuation of this element.

        Use this method to get an approximation of the valuation when speed is
        more important than accuracy.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: K(1025)._lower_bound()
            0

        """
        return self._valuation.lower_bound(self._x)
        

class BaseElement_Ring(BaseElement_base, HenselizationElement_Ring):
    r"""
    An element of :class:`Henselization_Ring` which is in one of its uncompleted
    ``base`` fields.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = ZZ.valuation(2)
        sage: R = Henselization(ZZ, v)
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = ZZ.valuation(2)
            sage: R = Henselization(ZZ, v)
            sage: R(3) // R(2) # indirect doctest
            1

        """
        if isinstance(other, BaseElement_Ring):
            if self._base is other._base and self._valuation is other._valuation:
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

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = ZZ.valuation(2)
            sage: R = Henselization(ZZ, v)
            sage: R(3) % R(2) # indirect doctest
            1

        """
        return self - (self // other) * other


class BaseElement_Field(BaseElement_base, HenselizationElement_Field):
    r"""
    An element of :class:`Henselization_Field` which is in one of its uncompleted
    ``base`` fields.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = ZZ.valuation(2)
        sage: R = Henselization(ZZ, v)
        sage: R(0)
        0

    """
    def __init__(self, parent, base, valuation, x):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: isinstance(K(0), BaseElement_Field)
            True

        """
        super(BaseElement_Field, self).__init__(parent, base, valuation, x)
