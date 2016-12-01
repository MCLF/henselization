# -*- coding: utf-8 -*-
r"""
Exact completions of rings

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

from sage.rings.ring import Ring, PrincipalIdealDomain, Field
from sage.structure.factory import UniqueFactory
from sage.misc.lazy_attribute import lazy_attribute

class CompletionFactory(UniqueFactory):
    r"""
    The completion of ``R`` with respect to ``v``.
    
    INPUT:
    
    - ``R`` -- a field or an excellent integral domain whose localization at
      the maximal ideal of ``v`` is the discrete valuation ring of ``v``.

    - ``v`` -- a non-trivial discrete valuation on ``R``.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: Completion(QQ, v)
        Completion of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, R, v):
        r"""
        Create a key which uniquely identifies this completion.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: Completion(QQ, v) is Completion(QQ, v) # indirect doctest
            True

        """
        from sage.categories.all import IntegralDomains
        if R not in IntegralDomains():
            raise ValueError("R must be a ring")
        from mac_lane import DiscretePseudoValuationSpace
        if v not in DiscretePseudoValuationSpace(R) or not v.is_discrete_valuation():
            raise ValueError("v must be a discrete valuation on R")
        if v.is_trivial():
            raise ValueError("v must not be trivial")

        return R, v

    def create_object(self, version, key):
        r"""
        Create the completion identified by ``key``.

        TESTS::
            
            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: Completion(QQ, v) # indirect doctest
            Completion of Rational Field with respect to 5-adic valuation

        """
        R, v = key
        from sage.categories.all import Fields
        if v.value_semigroup().is_group():
            return CompleteField(R, v)
        else:
            return CompleteDomain(R, v)

Completion = CompletionFactory("Completion")

class ExtensionFactory(UniqueFactory):
    r"""
    Return the algebraic extension which adjoins to ``base`` a root of
    ``polynomial``.

    Do not call this factory directly, but call ``base.extension()`` instead.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: K.extension(x^2 - 5) # indirect doctest
        Extension defined by x^2 - 5 of Completion of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, base, polynomial, name):
        r"""
        Return a key that uniquely defines this extension.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: K.extension(x^2 - 5) is K.extension(x^2 - 5) # indirect doctest
            True

        """
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        if not is_Polynomial(polynomial):
            raise TypeError("polynomial must be a polynomial")
        if len(polynomial.parent().gens()) != 1 or polynomial.parent().base() is not base:
            raise ValueError("polynomial must be an element of a univariate polynomial ring over %r but %r is not"%(base, polynomial.parent()))
        if not polynomial.is_irreducible():
            raise ValueError("polynomial must be irreducible but %r is not"%(polynomial,))
        polynomial.change_variable_name(name)

        return base, polynomial

    def create_object(self, version, key):
        r"""
        Return the extension defined by ``key``.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: K.extension(x^2 - 5) # indirect doctest
            Extension defined by x^2 - 5 of Completion of Rational Field with respect to 5-adic valuation

        """
        base, polynomial = key

        # Currently, we only support extensions which are defined by
        # polynomials which already live over the exact base ring of the
        # completion.
        base_polynomial = polynomial.map_coefficients(base.base().convert_map_from(base))
        base_extension = base.base().extension(base_polynomial, names=(polynomial.variable_name(),))
        base_valuation = base._original_valuation.extension(base_extension)

        from sage.categories.all import Fields, IntegralDomains
        if base in Fields():
            return CompleteExtensionField(base, base_extension, base_valuation, polynomial, base_polynomial, name=polynomial.variable_name())
        elif base in IntegralDomains():
            return CompleteExtensionDomain(base, base_extension, base_valuation, polynomial, base_polynomial, name=polynomial.variable_name())
        raise NotImplementedError

Extension = ExtensionFactory("Extension")

class CompleteRing_base(Ring):
    r"""
    Abstract base class for the completion of ``R`` with respect to ``v``.

    INPUT:

    - ``R`` -- an excellent integral domain whose localization at the elements
      of positive valuation is the valuation ring of ``v``.

    - ``v`` -- a non-trivial discrete valuation on ``R``

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: Completion(ZZ, v)
        Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, R, v, category):
        r"""
        TESTS::

            sage: from completion import *
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x)
            sage: C = Completion(K, v)
            sage: isinstance(C, CompleteRing_base)
            True
            sage: TestSuite(C).run() # long time
    
        """
        Ring.__init__(self, R, category=category)

        self._original_valuation = v

        # The completion contains the inverses of elements of valuation zero,
        # therefore we extend the valuation to the field of fractions.
        fraction_field = R.fraction_field()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        if is_PolynomialRing(R) and R.ngens() == 1:
            # For polynomial rings, the support for valuations on function
            # fields is better than the support for naive fraction fields of
            # polynomial rings; the two are virtually identical, so we use
            # function fields.
            from sage.rings.all import FunctionField
            fraction_field = FunctionField(R.base().fraction_field(), names=(R.variable_name(),))
        self._base_valuation = v.extension(fraction_field)

        from sage.rings.all import PolynomialRing
        self._polynomial_ring = PolynomialRing(self, 'x')

        # monkey patch the broken _gcd_univariate_polynomial from UniqueFactorizationDomains
        self._gcd_univariate_polynomial_original = self._gcd_univariate_polynomial
        self._gcd_univariate_polynomial = lambda f, g: f.parent()(self._gcd_univariate_polynomial_original(f,g))

        from sage.categories.homset import Hom
        from sage.categories.morphism import SetMorphism
        R.register_conversion(SetMorphism(Hom(self,R), lambda x: R(x._x)))
        if R is not fraction_field:
            fraction_field.register_conversion(SetMorphism(Hom(self,fraction_field), lambda x: fraction_field(x._x)))

    def characteristic(self):
        r"""
        Return the characteristic of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.characteristic()
            0

        """
        return self.base().characteristic()

    def is_finite(self):
        r"""
        Return whether this ring is finite.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.is_finite()
            False

        """
        # since the underlying valuation is non-trivial, there must be
        # infinitely many elements
        return False

    def uniformizer(self):
        r"""
        Return a uniformizing element for :meth:`valuation`.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.uniformizer()
            2

        """
        return self(self._base_valuation.uniformizer())

    def _coerce_map_from_(self, other):
        r"""
        Return a coercion from ``other`` to this ring if one exists.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.has_coerce_map_from(QQ) # indirect doctest
            True

        """
        if self.base().has_coerce_map_from(other):
            return True
        from sage.structure.parent import is_Parent
        if is_Parent(other):
            if other.base().is_subring(self.base()):
                if isinstance(other, CompleteRing_base):
                    try:
                        extension = other._base_valuation.extension(self.base())
                    except NotImplementedError:
                        pass
                    else:
                        if self._base_valuation == other._base_valuation.extension(self.base()):
                        from maps import BaseExtensionCoercion
                        return BaseExtensionCoercion(other.Hom(self))
        return super(CompleteRing_base, self)._coerce_map_from_(other)


    def _an_element_(self):
        r"""
        Return an element of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.an_element()
            2

        """
        return self(self._base_valuation.uniformizer())

    def _element_constructor_(self, x):
        r"""
        Create an element in this parent from ``x``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0) # indirect doctest
            sage: x.parent() is K
            True

        """
        from base_element import BaseElement_base
        if isinstance(x, BaseElement_base):
            x = x._x
        if x in self.base():
            x = self.base()(x)
            return self._base_element_class(self, x)
        if x in self.base().fraction_field():
            x = self.base().fraction_field()(x)
            from sage.categories.fields import Fields
            if self in Fields() or self._base_valuation(x) >= 0:
                return self._base_element_class(self, x)
        if isinstance(x, tuple):
            if len(x) == 2:
                v, i = x
                from sage.rings.all import NN
                from mac_lane.limit_valuation import MacLaneLimitValuation
                if isinstance(v, MacLaneLimitValuation) and v.domain() is self._polynomial_ring and i in NN:
                    return self._mac_lane_element_class(self, v, i)
        raise ValueError("Can not convert %r to an element in %r"%(x, self))

    @lazy_attribute
    def _mac_lane_element_class(self):
        r"""
        The class for elements which are given as coefficients of key
        polynomial of limit valuations.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: F = f.factor() # indirect doctest
            sage: isinstance(F[0][0][0], K._mac_lane_element_class)
            True

        """
        from mac_lane_element import MacLaneElement
        return self.__make_element_class__(MacLaneElement)

    def _repr_(self):
        r"""
        Return a printable representation of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: Completion(QQ, v)
            Completion of Rational Field with respect to 2-adic valuation

        """
        return "Completion of %r with respect to %r"%(self.base_ring(), self._original_valuation)

    def some_elements(self):
        r"""
        Return some typical elements of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: len(K.some_elements())
            100

        """
        return map(self, self.base_ring().some_elements())

    def valuation(self):
        r"""
        Return the valuation on this complete ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation()
            2-adic valuation

        """
        from .valuation import Valuation
        return Valuation(self)

    def residue_field(self):
        r"""
        Return the residue field, i.e., the elements of non-negative valuation
        module elements of positive valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.residue_field()
            Finite Field of size 2

        """
        return self.valuation().residue_field()

    def _factor_univariate_polynomial(self, f):
        r"""
        Return the factorization of ``f`` over this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x + 1
            sage: f.factor() # indirect doctest
            x + 1
            sage: f = x^2 + 1
            sage: f.factor()
            x^2 + 1

        A non-trivial example::

            sage: G = GaussianIntegers().fraction_field()
            sage: v = pAdicValuation(G, 2)
            sage: K = Completion(G, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: f.factor() # long time
            (x + (-I - 8) + O(?)) * (x + (5*I - 4) + O(?))

        Another non-trivial example::

            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: f.factor() # long time
            (x + 57 + O(?)) * (x + 68 + O(?))

        """
        if f.is_constant():
            raise NotImplementedError
        if not f.is_monic():
            raise NotImplementedError
        if not f.is_squarefree():
            # TODO: piece together the factorization of the factors of the squarefree decomposition
            raise NotImplementedError

        from sage.structure.factorization import Factorization
        G = self._polynomial_ring(f)
        approximants = self.valuation().mac_lane_approximants(G)
        if len(approximants) == 1:
            return Factorization([(G, 1)])
        approximants = sum([approximant.mac_lane_step(G) for approximant in approximants], [])
        factors = []
        for approximant in approximants:
            # we only need to perform a Mac Lane step when there is ramification
            # introduced in the last step (to get the degrees right), anyway, it
            # does not hurt to do another step
            approximant = approximant.mac_lane_step(G)
            assert len(approximant)==1
            approximant = approximant[0]

            from sage.rings.all import infinity
            if approximant(approximant.phi()) == infinity:
                factors.append(approximant.phi())
                continue

            degree = approximant.phi().degree()
            from mac_lane.limit_valuation import LimitValuation
            limit = LimitValuation(approximant, G)
            coefficients = [self((limit, d)) for d in range(degree + 1)]
            if f.is_monic():
                coefficients[-1] = self(1)
            factor = f.parent()(coefficients)
            factors.append(factor)
        return Factorization([(factor, 1) for factor in factors], unit=self.one(), sort=False)

    def extension(self, poly, names=None, name=None):
        r"""
        Return the algebraic extension of this ring by adjoining a root of
        the irreducible polynomial ``poly``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: K.extension(x^2 + x + 1)
            Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

        """
        if not names is None:
            name = names
        if isinstance(name, tuple):
            name = name[0]
        if name is None:
            name = poly.parent().variable_name()

        return Extension(self, poly, name)

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        EXAMPLES:

        This ring is generated by its one element::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.ngens()
            1
            sage: K.gens()
            (1,)

        """
        return 1

    def gen(self, i):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES:

        This ring is generated by its one element::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.gen(0)
            1
            sage: K.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: ring has only one generator

        """
        if i == 0:
            return self.one()
        raise ValueError("ring has only one generator")
        

class CompleteDomain(CompleteRing_base, PrincipalIdealDomain):
    r"""
    The completion of the integral domain ``R`` with respect to ``v``.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: Completion(ZZ, v)
        Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, R, v, category = None):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: C = Completion(ZZ, v)
            sage: isinstance(C, CompleteDomain)
            True
            sage: TestSuite(C).run() # long time
            
        """
        if category is None:
            from completions import CompleteDiscreteValuationRings
            category = CompleteDiscreteValuationRings()
        CompleteRing_base.__init__(self, R, v, category)
        PrincipalIdealDomain.__init__(self, R, category=self.category())
        from sage.categories.fields import Fields
        if R in Fields():
            raise TypeError("R must not be a field")

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the base ring.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: x = R(0) # indirect doctest
            sage: isinstance(x, R._base_element_class)
            True

        """
        from base_element import BaseElementRing
        return self.__make_element_class__(BaseElementRing)

    def is_field(self, *args, **kwargs):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: R.is_field()
            False

        """
        return False

    def fraction_field(self):
        r"""
        Return the fraction field of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: R.fraction_field()
            Completion of Rational Field with respect to 2-adic valuation

        """
        return Completion(self._base_valuation.domain(), self._base_valuation)


class CompleteField(CompleteRing_base, Field):
    r"""
    The completion of the field ``K`` with respect to ``v``.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: Completion(QQ, v)
        Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, K, v, category = None):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: isinstance(K, CompleteField)
            True
            sage: TestSuite(K).run() # long time

        """
        if category is None:
            from completions import CompleteDiscreteValuationFields
            category = CompleteDiscreteValuationFields()
        CompleteRing_base.__init__(self, K, v, category)
        Field.__init__(self, K, category=self.category())

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the base ring.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0) # indirect doctest
            sage: isinstance(x, K._base_element_class)
            True

        """
        from base_element import BaseElementField
        return self.__make_element_class__(BaseElementField)


class CompleteExtension_base(Ring):
    r"""
    Abstract base class for the extension by adjunction of a root of
    ``polynomial`` to the complete ring ``base`.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: K.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base, polynomial, name, category=None):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: isinstance(L, CompleteExtension_base)
            True
            sage: TestSuite(L).run() # long time

        """
        Ring.__init__(self, base, category=category or base.category())

        self._polynomial = polynomial
        self._name = name
        self._base_ring = base

    def base_ring(self):
        r"""
        Return the base ring of this ring, i.e., the ring from which this
        algebraic extension has been created.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.base_ring() is K
            True

        """
        return self._base_ring

    def _repr_(self):
        r"""
        Return a printable representation of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: K.extension(x^2 + x + 1) # indirect doctest
            Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

        """
        return "Extension defined by %r of %r"%(self._polynomial, self.base_ring())

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        EXAMPLES:

        This extension is generated by a root of its defining polynomial::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.ngens()
            1
            sage: L.gens()
            (x,)

        """
        return 1

    def gen(self, i):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.gen(0)
            x
            sage: L.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: ring has only one generator

        """
        if i == 0:
            return self.base().gen(0)
        raise ValueError("ring has only one generator")
        

class CompleteExtensionDomain(CompleteExtension_base, CompleteDomain):
    r"""
    Extension by adjunction of a root of ``polynomial`` to the complete
    integral domain ``base`.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: S = Completion(ZZ, v)
        sage: R.<x> = S[]
        sage: S.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base, base_extension, valuation, polynomial, base_polynomial, name, category=None):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: S = Completion(ZZ, v)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: isinstance(T, CompleteExtensionDomain)
            True
            sage: TestSuite(T).run() # long time

        """
        CompleteExtension_base.__init__(self, base, polynomial, category)
        CompleteDomain.__init__(self, base_extension, valuation, category=self.category())


class CompleteExtensionField(CompleteExtension_base, CompleteField):
    r"""
    Extension by adjunction of a root of ``polynomial`` to the complete field
    ``base`.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: K.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base, base_extension, valuation, polynomial, base_polynomial, name, category=None):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: isinstance(L, CompleteExtensionField)
            True
            sage: TestSuite(L).run() # long time

        """
        CompleteExtension_base.__init__(self, base, polynomial, category)
        CompleteField.__init__(self, base_extension, valuation, category=self.category())
