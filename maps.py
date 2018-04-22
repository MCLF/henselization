# -*- coding: utf-8 -*-
r"""
Maps from and to Henselizations of rings

AUTHORS:

- Julian Rüth (2016-12-01): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.morphism import Morphism
from sage.structure.unique_representation import UniqueRepresentation

class ConvertMap_generic(Morphism):
    r"""
    Conversion map for codomains which can handle elements in the fraction
    field of the base of a Henselization.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(5)
        sage: K = Henselization(QQ, v)
        sage: QQ.convert_map_from(K)
        Generic morphism:
            From: Henselization of Rational Field with respect to 5-adic valuation
            To:   Rational Field

    """
    def _call_(self, x):
        r"""
        Evaluate this map at ``x``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(5)
            sage: K = Henselization(QQ, v)
            sage: f = QQ.convert_map_from(K)
            sage: f(K(0))
            0

        """
        from base_element import BaseElement_base
        if isinstance(x, BaseElement_base):
            return self.codomain()(x._x)
        raise NotImplementedError("Conversion of %s into %s"%(x, self.codomain()))


class ExtensionCoercion_generic(ConvertMap_generic):
    r"""
    Coercion map from a Henselization to an algebraic extension.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: L.coerce_map_from(K)
        Generic morphism:
            From: Henselization of Rational Field with respect to 2-adic valuation
            To:   Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    """
    def is_injective(self):
        r"""
        Return whether this coercion is injective.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = ZZ.valuation(2)
            sage: S = Henselization(ZZ, v)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: T.coerce_map_from(S).is_injective()
            True

        """
        return True


class VectorSpaceHenselizationIsomorphism(Morphism):
    r"""
    Base class for isomorphisms of Henselizations and vector spaces.

    TESTS::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: f = L.module()[1]
        sage: isinstance(f, VectorSpaceHenselizationIsomorphism)
        True

    """
    def _repr_type(self):
        r"""
        Return the type of this morphism.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.module()[1]; f # indirect doctest
            Isomorphism morphism:
              From: Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
              To:   Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

        """
        return "Isomorphism"

    def is_injective(self):
        r"""
        Return that this isomorphism is injective.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.module()[1]
            sage: f.is_injective()
            True

        """
        return True

    def is_surjective(self):
        r"""
        Return that this isomorphism is surjective.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.module()[1]
            sage: f.is_surjective()
            True

        """
        return True

    def __nonzero__(self):
        r"""
        Return that this isomorphism is not the zero morphism.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.module()[1]
            sage: bool(f) # indirect doctest
            True

        """
        return True


class VectorSpaceToHenselization(VectorSpaceHenselizationIsomorphism, UniqueRepresentation):
    r"""
    An isomorphism from a vector space to a complete field.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: f = L.module()[1]; f
        Isomorphism morphism:
          From: Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
          To:   Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    TESTS:

    Run the test suite but skip the tests that require the zero element of the
    hom space (not implemented)::

        # TODO: for some strange reason pickling does not work - I suppose that
        # this is somehow an artifact of running this standalone.
        sage: TestSuite(f).run(skip=("_test_nonzero_equal", "_test_pickling"))

    """
    from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
    __metaclass__ = InheritComparisonClasscallMetaclass

    def __init__(self, parent, basis):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.module()[1]
            sage: isinstance(f, VectorSpaceToHenselization)
            True

        """
        Morphism.__init__(self, parent)
        UniqueRepresentation.__init__(self)

        self._basis = basis

    def _call_(self, v):
        r"""
        Evaluate this morphism at ``v``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: V,f,_ = L.module(base=K)
            sage: f(V((1,2)))
            2*x + 1

        """
        return sum(self.codomain()(x)*b for x,b in zip(v,self._basis) if x != 0)


class HenselizationToVectorSpace(VectorSpaceHenselizationIsomorphism, UniqueRepresentation):
    r"""
    An isomorphism from a complete field to a vector space.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: f = L.module()[2]; f
        Isomorphism morphism:
          From: Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
          To:   Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    TESTS:

    Run the test suite but skip the tests that require the zero element of the
    hom space (not implemented)::

        # TODO: for some strange reason pickling does not work - I suppose that
        # this is somehow an artifact of running this standalone.
        sage: TestSuite(f).run(skip=("_test_nonzero_equal", "_test_pickling"))

    """
    from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
    __metaclass__ = InheritComparisonClasscallMetaclass

    def __init__(self, parent, base):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.module()[2]
            sage: isinstance(f, HenselizationToVectorSpace)
            True

        """
        Morphism.__init__(self, parent)
        UniqueRepresentation.__init__(self)

        self._base = base

    def _call_(self, x):
        r"""
        Evaluate this morphism at ``x``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: _,_,f = L.module(base=K)
            sage: f(L.gen() + 2)
            (2, 1)

        """
        return self.codomain()(x._vector_(base=self._base))


class RelativeExtensionCoercion_generic(Morphism):
    r"""
    A coercion between extensions of complete rings which extend each other.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: M = L.extension(y^2 + 2)
        sage: f = M.coerce_map_from(L); f
        Generic morphism:
          From: Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
          To:   Extension defined by y^2 + 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    TESTS::

        sage: isinstance(f, RelativeExtensionCoercion_generic)
        True
        sage: TestSuite(f).run()

    """
    def is_injective(self):
        r"""
        Return whether this coercion is injective, which is the case since it
        is the embedding of a ring extension.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 + 2)
            sage: M.coerce_map_from(L).is_injective()
            True

        """
        return True

    def is_surjective(self):
        r"""
        Return whether this coercion is surjective, which is only the case for
        trivial extensions::

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 + 2)
            sage: M.coerce_map_from(L).is_surjective()
            False

            sage: R.<z> = M[]
            sage: N = M.extension(z - 1)
            sage: N.coerce_map_from(M).is_surjective()
            True

        ^"""
        if self.domain() is self.codomain().base_ring():
            return self.codomain().degree() == 1
        raise NotImplementedError("surjectivity only for simple extensions")

    def _call_(self, x):
        r"""
        Evaluate this map at ``x``.

        EXAMPLES:

        This is currently not implemented for non-trivial cases (to implement
        this, we would need to make sure that we are choosing
        the roots of the defining polynomial consistently.)::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L.<x> = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 + 2)
            sage: M.coerce(x)
            Traceback (most recent call last):
            ...
            NotImplementedError: Selection of approximate root of x^2 + x + 1 in Extension defined by y^2 + 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

        """
        from base_element import BaseElement_base
        if isinstance(x, BaseElement_base):
            if self.domain()._base.is_subring(self.codomain()._base):
                return self.codomain()(x._x)
            minpoly = x._x.minpoly()
            if self.codomain().has_coerce_map_from(minpoly.base_ring()):
                if minpoly.degree() == 1:
                    return self.codomain()(-minpoly[0])
                raise NotImplementedError("Selection of approximate root of %s in %s"%(minpoly, self.codomain()))
        raise NotImplementedError("Coercion of %s into %s"%(x, self.codomain()))

    def _hash_(self):
        r"""
        Return a hash value for this morphism.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 + 2)
            sage: f = M.coerce_map_from(L)
            sage: hash(f) == hash(loads(dumps(f)))
            True

        """
        return hash((self.domain(), self.codomain()))

    def _richcmp_(self, other, op):
        r"""
        Compare this element to ``other``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = QQ.valuation(2)
            sage: K = Henselization(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 + 2)
            sage: f = M.coerce_map_from(L)
            sage: f == loads(dumps(f))
            True

        """
        from sage.structure.richcmp import op_EQ, op_NE
        if type(self) is not type(other):
            raise NotImplementedError("Comparison operators only implemented for morphisms with the same underlying implementation.")

        if op == op_EQ:
            return self.domain() is other.domain() and self.codomain() is other.codomain()
        if op == op_NE:
            return not self == other
        raise NotImplementedError("Operator not implemented for this morphism")

class QuotientConversion_generic(Morphism):
    r"""
    A conversion between two quotients that is induced by a conversion between
    their bases.
    """
    def _call_(self, x):
        r"""
        Evaluate this morphism at ``x``.
        """
        return x._x.lift().map_coefficients(self.codomain().base_ring().convert_map_from(self.domain().base_ring()))(self.codomain().gen())
