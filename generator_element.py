# -*- coding: utf-8 -*-
r"""
Symbolic elements which generate an extension of a Henselization

AUTHORS:

- Julian Rüth (2017-04-30): initial version

"""
#*****************************************************************************
#       Copyright (C) 2017-2018 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.element import IntegralDomainElement
from sage.misc.cachefunc import cached_method
from henselization_element import HenselizationElement_base

class GeneratorElement(HenselizationElement_base):
    r"""
    Element class for generators of :class:`HenselizationExtension_base` which
    merely exist symbolically as the root of a certain polynomial but which
    support essentially no arithmetic.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(3)
        sage: K = Henselization(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 - 3)
        sage: R.<y> = L[]
        sage: M = L.extension(y^3 - 3)
        sage: M.gen()
        y

    In towers of extensions, these elements are also used to make explicit the
    choice of embedding (however, this is not implemented yet)::

        sage: M(L.gen()) == -M(L.gen())
        Traceback (most recent call last):
        ...
        NotImplementedError: Selection of approximate root of x^2 - 3 in Extension defined by y^3 - 3 of Extension defined by x^2 - 3 of Henselization of Rational Field with respect to 3-adic valuation

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(3)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 - 3)
            sage: R.<y> = L[]
            sage: M = L.extension(y^3 - 3)
            sage: y = M.gen()
            sage: isinstance(y, GeneratorElement)
            True
            sage: TestSuite(y).run()

        """
        super(GeneratorElement, self).__init__(parent)

    def _repr_(self):
        r"""
        Return a printable representation of this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(3)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 - 3)
            sage: R.<y> = L[]
            sage: M = L.extension(y^3 - 3)
            sage: M.gen()
            y

        """
        return self.parent().variable_name()

    def _richcmp_(self, other, op):
        r"""
        Compare this elment to ``other`` with respect to ``op``.

        TESTS::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(3)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 - 3)
            sage: R.<y> = L[]
            sage: M = L.extension(y^3 - 3)
            sage: y = M.gen()

        Generators can not be elements in the base ring::

            sage: y == 0
            False

        """
        if op == 2:
            from base_element import BaseElement_base
            if isinstance(other, BaseElement_base):
                return False
            if isinstance(other, GeneratorElement) and other.parent() == self.parent():
                return True
            raise NotImplementedError("comparison of generator %r to %r"%(self, other))
        elif op == 3:
            return not (self == other)
        raise NotImplementedError

    @cached_method
    def valuation(self):
        r"""
        Return the valuation of this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(3)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 - 3)
            sage: R.<y> = L[]
            sage: M = L.extension(y^3 - 3)
            sage: y = M.gen()
            sage: y.valuation()
            1/3

        """
        return self.parent()._polynomial[0].valuation() / self.parent()._polynomial.degree()

    @cached_method
    def reduction(self):
        r"""
        Return the reduction of this element in the :meth:`residue_ring`.

        EXAMPLES:

        This is not implemented yet::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(3)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 - 3)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 - y - 1)
            sage: y = M.gen()
            sage: y.reduction()
            Traceback (most recent call last):
            ...
            NotImplementedError: approximations of generators

        """
        return self.approximation(self.parent().valuation().value_group().gen()).reduction()

    @cached_method
    def approximation(self, precision):
        r"""
        Return an approximation of this element which only differs by an
        element of valuation at most ``precision``.

        EXAMPLES:

        This is not implemented yet::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(3)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 - 3)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 - y - 1)
            sage: y = M.gen()
            sage: y.approximation(precision=0)
            Traceback (most recent call last):
            ...
            NotImplementedError: approximations of generators

        """
        raise NotImplementedError("approximations of generators")
