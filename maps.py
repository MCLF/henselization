# -*- coding: utf-8 -*-
r"""
Maps from and to completions of rings

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

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the completion/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.path.dirname(os.getcwd()))

from sage.categories.morphism import Morphism
from sage.structure.unique_representation import UniqueRepresentation

class ConvertMap_generic(Morphism):
    r"""
    Conversion map for codomains which can handle elements in the fraction
    field of the base of a completion.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: QQ.convert_map_from(K)
        Generic morphism:
            From: Completion of Rational Field with respect to 5-adic valuation
            To:   Rational Field

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: f = QQ.convert_map_from(K)
            sage: isinstance(f, ConvertMap_generic)
            True

        """
        super(ConvertMap_generic, self).__init__(parent)

    def _call_(self, x):
        r"""
        Evaluate this map at ``x``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: f = QQ.convert_map_from(K)
            sage: f(K(0))
            0

        """
        from base_element import BaseElement_base
        if isinstance(x, BaseElement_base):
            return self.codomain()(x._x)
        raise NotImplementedError

class ExtensionCoercion_generic(ConvertMap_generic):
    r"""
    Coercion map from a completion to an algebraic extension.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: L.coerce_map_from(K)
        Generic morphism:
            From: Completion of Rational Field with respect to 2-adic valuation
            To:   Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

    """
    def is_injective(self):
        r"""
        Return whether this coercion is injective.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: S = Completion(ZZ, v)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: T.coerce_map_from(S).is_injective()
            True

        """
        return True

class VectorSpaceCompletionIsomorphism(Morphism):
    r"""
    Base class for isomorphisms of completions and vector spaces.

    TESTS::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: f = L.vector_space()[1]
        sage: isinstance(f, VectorSpaceCompletionIsomorphism)
        True

    """
    def _repr_type(self):
        r"""
        Return the type of this morphism.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.vector_space()[1]; f # indirect doctest
            Isomorphism morphism:
              From: Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation
              To:   Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

        """
        return "Isomorphism"

    def is_injective(self):
        r"""
        Return that this isomorphism is injective.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.vector_space()[1]
            sage: f.is_injective()
            True

        """
        return True

    def is_surjective(self):
        r"""
        Return that this isomorphism is surjective.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.vector_space()[1]
            sage: f.is_surjective()
            True

        """
        return True

    def __nonzero__(self):
        r"""
        Return that this isomorphism is not the zero morphism.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.vector_space()[1]
            sage: bool(f) # indirect doctest
            True

        """
        return True


class VectorSpaceToCompletion(VectorSpaceCompletionIsomorphism, UniqueRepresentation):
    r"""
    An isomorphism from a vector space to a complete field.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: f = L.vector_space()[1]; f
        Isomorphism morphism:
          From: Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation
          To:   Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

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

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.vector_space()[1]
            sage: isinstance(f, VectorSpaceToCompletion)
            True

        """
        Morphism.__init__(self, parent)
        UniqueRepresentation.__init__(self)

        self._basis = basis

    def _call_(self, v):
        r"""
        Evaluate this morphism at ``v``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: V,f,_ = L.vector_space(base=K)
            sage: f(V((1,2)))
            2*x + 1

        """
        return sum(x*b for x,b in zip(v,self._basis))


class CompletionToVectorSpace(VectorSpaceCompletionIsomorphism, UniqueRepresentation):
    r"""
    An isomorphism from a complete field to a vector space.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: f = L.vector_space()[2]; f
        Isomorphism morphism:
          From: Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation
          To:   Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

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

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: f = L.vector_space()[2]
            sage: isinstance(f, CompletionToVectorSpace)
            True

        """
        Morphism.__init__(self, parent)
        UniqueRepresentation.__init__(self)

        self._base = base

    def _call_(self, x):
        r"""
        Evaluate this morphism at ``x``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: _,_,f = L.vector_space(base=K)
            sage: f(L.gen() + 2)
            (2, 1)

        """
        return self.codomain()(x._vector_(base=self._base))
