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
from sage.rings.ring import CommutativeRing, Field
from sage.structure.factory import UniqueFactory
from sage.misc.lazy_attribute import lazy_attribute

from sage.categories.fields import Fields
from sage.structure.element import is_Element

class CompletionFactory(UniqueFactory):
    r"""
    Return the completion of ``R`` with respect to ``v``.
    
    INPUT:
    
    - ``R`` -- a field or an excellent integral domain whose localization at
      the maximal ideal of ``v`` is the discrete valuation ring of ``v``. (Most
      of these conditions are not verified by this factory.)

    - ``v`` -- a non-trivial discrete valuation on ``R``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: Completion(QQ, v)
        Completion of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, R, v):
        r"""
        Create a key which uniquely identifies this completion.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: Completion(QQ, v) is Completion(QQ, v) # indirect doctest
            True

        """
        from sage.categories.all import IntegralDomains
        if R not in IntegralDomains():
            raise ValueError("R must be an integral domain")
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
            
            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: Completion(QQ, v) # indirect doctest
            Completion of Rational Field with respect to 5-adic valuation

        """
        R, v = key

        if v.value_semigroup().is_group():
            return Completion_Field(R, v)
        else:
            return Completion_Ring(R, v)

Completion = CompletionFactory("Completion")


class ExtensionFactory(UniqueFactory):
    r"""
    Return the algebraic extension which adjoins to ``base`` a root of
    ``polynomial``.

    Do not call this factory directly, but call ``base.extension()`` instead.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: K.extension(x^2 - 5) # indirect doctest
        Extension defined by x^2 - 5 of Completion of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, base, polynomial, name, check=True):
        r"""
        Return a key that uniquely defines this extension.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
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
        if polynomial.is_constant():
            raise ValueError("polynomial must not be constant")
        if any(c not in base.base() for c in polynomial.coefficients(sparse=False)):
            raise NotImplementedError("polynomial must have coefficients in %r"%(base.base(),))

        polynomial = polynomial.change_variable_name(name)

        if check and not polynomial.is_squarefree():
            # We only check squarefreeness here. Irreducibility is checked
            # automatically, when the extensions of the valuations on base to
            # the ring are constructed. (If there is more than one extension,
            # i.e., the polynomial is not irreducible, extension() is going to
            # complain.)
            raise ValueError("polynomial must be irreducible but %r is not"%(polynomial,))

        return base, polynomial

    def create_object(self, version, key):
        r"""
        Return the extension defined by ``key``.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: K.extension(x^2 - 5) # indirect doctest
            Extension defined by x^2 - 5 of Completion of Rational Field with respect to 5-adic valuation

        """
        base, polynomial = key

        from sage.categories.all import Fields, IntegralDomains
        if base.valuation().value_semigroup().is_group():
            return CompletionExtension_Field(base, polynomial)
        else:
            return CompletionExtension_Ring(base, polynomial)

Extension = ExtensionFactory("Extension")


class Completion_base(CommutativeRing):
    r"""
    Abstract base class for the completion of ``base`` with respect to
    ``base_valuation``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: Completion(ZZ, v)
        Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base, base_valuation, category):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: K.<x> = FunctionField(QQ)
            sage: v = FunctionFieldValuation(K, x)
            sage: C = Completion(K, v)
            sage: isinstance(C, Completion_base)
            True
            sage: TestSuite(C).run() # long time
    
        """
        super(Completion_base, self).__init__(base_ring=base, category=category)

        self._base_valuation = base_valuation

        # The completion contains not only the elements of base but many more
        # elements of the fraction field of base, namely the elements of
        # valuation zero. We therefore, extend base_valuation to that field of
        # fractions.
        self._base_fraction_field = self.base().fraction_field()
        from sage.rings.polynomial.polynomial_ring import is_PolynomialRing
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
        if is_PolynomialRing(base) and base.ngens() == 1:
            # For polynomial rings, the support for valuations on function
            # fields is better than the support for naive fraction fields of
            # polynomial rings; the two are virtually identical, so we use
            # function fields.
            from sage.rings.all import FunctionField
            self._base_fraction_field = FunctionField(base.base().fraction_field(), names=(base.variable_name(),))
        elif is_PolynomialQuotientRing(base) and base.ngens() == 1:
            # We could rewrite quotient rings to field extensions here.
            # However, they often have sever performance penalties (e.g. in the
            # case of number fields where the construction of a relative number
            # field can take a long time because it falls back on the
            # construction of the absolute number field.)
            # We therefore only rewrite the base field but keep the quotient.
            if base.base_ring().fraction_field() is not base.base_ring():
                self._base_fraction_field = base.base().change_ring(base.base_ring().fraction_field()).quo(base.modulus())
        self._base_fraction_field_valuation = self._base_valuation.extension(self._base_fraction_field)

        # monkey patch the broken _gcd_univariate_polynomial from UniqueFactorizationDomains
        if not hasattr(self, '_gcd_univariate_polynomial_original'):
            self._gcd_univariate_polynomial_original = self._gcd_univariate_polynomial
            self._gcd_univariate_polynomial = lambda f, g: f.parent()(self._gcd_univariate_polynomial_original(f,g))

        # provide conversions from the completion back to its base and its field of fractions
        from maps import ConvertMap_generic
        from sage.categories.all import SetsWithPartialMaps
        # TODO: use weak references
        homspace = self.Hom(self._base, category=SetsWithPartialMaps())
        self._base.register_conversion(homspace.__make_element_class__(ConvertMap_generic)(homspace))
        if self._base is not self._base_fraction_field:
            homspace = self.Hom(self._base_fraction_field, category=SetsWithPartialMaps())
            self._base_fraction_field.register_conversion(homspace.__make_element_class__(ConvertMap_generic)(homspace))

    def _gcd_univariate_polynomial(self, f, g):
        r"""
        Return the greatest common divisor of ``f`` and ``g``.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: (x^2 - 1).gcd(x - 1) # indirect doctest
            x - 1

        """
        if all(c in self._base_fraction_field for c in f.coefficients(sparse=False)) and all(c in self._base_fraction_field for c in g.coefficients(sparse=False)):
            # Since the gcd is independent of the base field, we can compute it
            # over the base field which is usually much faster (using a
            # multimodular algorithm for example) than the naive implementation
            return f.change_ring(self._base_fraction_field).gcd(g.change_ring(self._base_fraction_field)).change_ring(self)
        else:
            raise NotImplementedError

    def _xgcd_univariate_polynomial(self, f, g):
        r"""
        Return the xgcd of ``f`` and  ``g``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: x.xgcd(x) # indirect doctest
            (x, 0, 1)

        """
        if all(c in self._base_fraction_field for c in f.coefficients(sparse=False)) and all(c in self._base_fraction_field for c in g.coefficients(sparse=False)):
            # Since the gcd is independent of the base field, we can compute it
            # over the base field which is usually much faster (using a
            # multimodular algorithm for example) than the naive implementation
            h,s,t = f.change_ring(self._base_fraction_field).xgcd(g.change_ring(self._base_fraction_field))
            h = h.change_ring(self)
            s = s.change_ring(self)
            t = t.change_ring(self)
            return h,s,t
        else:
            raise NotImplementedError
    
    def base_ring(self):
        r"""
        Return the base ring of this ring.

        There are two different base rings for completions. The :meth:`base`
        which is the ring that the completion completes and the
        :meth:`base_ring` which is the ring from which the completion has been
        constructed.

        EXAMPLES:

        For completions that are not algebraical extensions of another
        completion, these two concepts coincide::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.base_ring() is K.base() is QQ
            True

        """
        return self.base()

    def vector_space(self, base=None):
        r"""
        Return a ``base``-vector space isomorphic to this field together with
        isomorphisms to and from this vector space.

        INPUT:

        - ``base`` -- a field of which this field is a finite extension
          (default: the field itself)

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.vector_space(base=L)
            (Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation
               To:   Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation
               To:   Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation)
            sage: L.vector_space(base=K)
            (Vector space of dimension 2 over Completion of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Vector space of dimension 2 over Completion of Rational Field with respect to 2-adic valuation
               To:   Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation
               To:   Vector space of dimension 2 over Completion of Rational Field with respect to 2-adic valuation)

        """
        if base is None:
            base = self
        basis = self._vector_space_basis(base=base)
        V = base**len(basis)
        from sage.all import Hom
        to_self_parent = Hom(V, self)
        from maps import VectorSpaceToCompletion, CompletionToVectorSpace
        to_self = to_self_parent.__make_element_class__(VectorSpaceToCompletion)(to_self_parent, basis)
        from_self_parent = Hom(self, V)
        from_self = from_self_parent.__make_element_class__(CompletionToVectorSpace)(from_self_parent, base)
        return (V, to_self, from_self)

    def _vector_space_basis(self, base):
        r"""
        Return a basis of this field as a vector space over ``base``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K._vector_space_basis(K)
            (1,)

        """
        if base is None:
            base = self
        if base is self:
            return (self(1),)
        raise NotImplementedError

    def characteristic(self):
        r"""
        Return the characteristic of this ring.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.has_coerce_map_from(QQ) # indirect doctest
            True

        """
        if self.base().has_coerce_map_from(other):
            return True
        
        from maps import ExtensionCoercion_generic
        if self.base_ring().has_coerce_map_from(other):
            homspace = other.Hom(self)
            return homspace.__make_element_class__(ExtensionCoercion_generic)(homspace)

        if isinstance(other, Completion_base):
            if other.base().is_subring(self.base()):
                try:
                    extension = other._base_valuation.extension(self.base())
                except NotImplementedError:
                    pass
                else:
                    if self._base_valuation == extension:
                        homspace = other.Hom(self)
                        return homspace.__make_element_class__(ExtensionCoercion_generic)(homspace)
        return super(Completion_base, self)._coerce_map_from_(other)

    def _an_element_(self):
        r"""
        Return an element of this ring.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.an_element()
            2

        """
        return self.uniformizer()

    def _element_constructor_(self, x):
        r"""
        Create an element in this parent from ``x``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0) # indirect doctest
            sage: x.parent() is K
            True

        """
        if x in self.base():
            x = self._base_fraction_field(x)
            return self._base_element_class(self, self._base_fraction_field, self._base_fraction_field_valuation, x)

        if x in self._base_fraction_field:
            x = self._base_fraction_field(x)
            if self in Fields() or self._base_fraction_field_valuation(x) >= 0:
                return self._base_element_class(self, self._base_fraction_field, self._base_fraction_field_valuation, x)

        if isinstance(x, tuple):
            if len(x) == 2:
                v, i = x
                from sage.rings.all import NN
                from mac_lane.limit_valuation import MacLaneLimitValuation
                if isinstance(v, MacLaneLimitValuation) and v.domain().base() is self and i in NN:
                    return self._mac_lane_element_class(self, v, i)
        raise ValueError("Can not convert %r to an element in %r"%(x, self))

    @lazy_attribute
    def _mac_lane_element_class(self):
        r"""
        The class for elements which are given as coefficients of key
        polynomial of limit valuations.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: Completion(QQ, v)
            Completion of Rational Field with respect to 2-adic valuation

        """
        return "Completion of %r with respect to %r"%(self.base_ring(), self._base_valuation)

    def some_elements(self):
        r"""
        Return some typical elements of this ring.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: len(K.some_elements())
            100

        """
        # At the moment we do not return any MacLaneElement elements. They do
        # not support any basic arithmetic (yet) and so make lots of TestSuite
        # tests fail.
        return map(self, self.base().some_elements())

    def valuation(self):
        r"""
        Return the valuation on this complete ring.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.residue_field()
            Finite Field of size 2

        """
        return self.valuation().residue_field()

    def extension(self, f, names=None, name=None, check=True):
        r"""
        Return the algebraic extension of this ring obtained by adjoining a
        root of the irreducible polynomial ``f``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: K.extension(x^2 + x + 1)
            Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

        """
        if names is not None:
            name = names
        if isinstance(name, tuple):
            name = name[0]
        if name is None:
            name = f.parent().variable_name()

        return Extension(self, f, name, check=check)

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        EXAMPLES:

        This ring is generated by its one element::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.ngens()
            1
            sage: K.gens()
            (1,)

        """
        return 1

    def gen(self, i=0):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES:

        This ring is generated by its one element::

            sage: sys.path.append(os.getcwd()); from completion import *
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

    def _factor_univariate_polynomial(self, f):
        r"""
        Return the factorization of ``f`` over this ring.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
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
            sage: f.factor()
            (x + I + O(?)) * (x + I)

        Another non-trivial example::

            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: f.factor()
            (x + 2 + O(?)) * (x + 3 + O(?))

        """
        from sage.misc.misc import verbose
        verbose("Factoring %r over %r"%(f, self))
        if f.is_constant():
            raise NotImplementedError("factorization of constant polynomials")
        if not f.is_monic():
            raise NotImplementedError("factorization of non-monic polynomials")
        if not f.is_squarefree():
            F = f.squarefree_decomposition()
            from sage.structure.factorization import Factorization
            factors = []
            unit = F.unit()
            for squarefree_factor,e in F:
                G = squarefree_factor.factor()
                unit *= G.unit()**e
                for factor, ee in G:
                    factors.append((factor, ee*e))
            return Factorization(factors, unit=unit, simplify=False, sort=False)

        from sage.structure.factorization import Factorization
        approximants = self.valuation().mac_lane_approximants(f, require_maximal_degree=True)
        if len(approximants) == 1:
            return Factorization([(f, 1)], sort=False)
        factors = []
        for approximant in approximants:
            from sage.rings.all import infinity
            if approximant(approximant.phi()) == infinity:
                factors.append(approximant.phi())
                continue

            degree = approximant.phi().degree()
            from mac_lane.limit_valuation import LimitValuation
            limit = LimitValuation(approximant, f)
            coefficients = [self((limit, d)) for d in range(degree + 1)]
            if f.is_monic():
                coefficients[-1] = self(1)
            factor = f.parent()(coefficients)
            factors.append(factor)
        return Factorization([(factor, 1) for factor in factors], unit=self.one(), sort=False, simplify=False)

    def ideal(self, *args, **kwds):
        r"""
        Essentially a monkey-patched version of ``ideal`` which only reduces
        the generators through gcd if the ring is a PrincipalIdealDomain (i.e.,
        it has that type). We are a principal ideal domain (by category) but do
        not inherit from that deprecated type.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: R.ideal(1, 1)
            Principal ideal (1) of Completion of Integer Ring with respect to 2-adic valuation

        """
        I = super(Completion_base, self).ideal(*args, **kwds)
        gens = I.gens()
        # Use GCD algorithm to obtain a principal ideal
        g = gens[0]
        if len(gens) == 1:
            try:
                g = g.gcd(g) # note: we set g = gcd(g, g) to "canonicalize" the generator: make polynomials monic, etc.
            except (AttributeError, NotImplementedError):
                pass
        else:
            for h in gens[1:]:
                g = g.gcd(h)
        gens = [g]
        return super(Completion_base, self).ideal(gens)

    def _krasner_bound(self, minpoly, assume_irreducible=False):
        r"""
        Return modulo which valuation ``minpoly`` defines a unique extension.

        More specifically, return a bound for each coefficient such that adding
        error terms of a *higher* valuation to ``minpoly`` does not change the
        isomorphism class of the extension obtained by adjoining a root of
        ``minpoly`` to this ring.

        ALGORITHM:

        Write `f` for the ``minpoly``, and let `\alpha` be a root of `f`.
        Let `f(x)=\sum a_i x^i` and write `g(x)=\sum b_i x^i` for an
        irreducible polynomial that is such that the coefficients `b_i` are
        close to the `a_i`. Let `\beta` be a root of `g`.
        By Krasner's Lemma `\beta` and `\alpha` generate isomorphic fields if
        `v(\alpha-\beta) > \max v(\alpha-\alpha')` for `\alpha'\neq\alpha`
        conjugates of `\alpha`.

        Note that
        `v(\beta - \alpha) + \sum v(\beta-\alpha')=v(f(\beta))=v(\sum (a_i-b_i)\beta^i)`.
        The right hand side is bounded from below by the minimum of the
        `v(a_i-b_i) + v(\beta^i)`.

        Our aim is now to determine from the above term `v(\sum(a_i-b_i)\beta^i)`
        whether the condition of the Krasner Lemma
        `v(\alpha-\beta) > \max v(\alpha-\alpha')` holds.
        It suffices that
        `v(\alpha-\beta) + \sum v(\alpha'-\beta) > \max v(\alpha-\alpha') + \sum v(\alpha-\alpha')`.
        because in the sum
        `v(\alpha-\alpha')\ge \max\{v(\alpha-\beta),v(\beta-\alpha')\}`.

        To compute the right hand side of the the above inequality, consider
        the polynomial `f(x-\alpha)` which can be written in its Taylor
        expansion as `f(x-\alpha)=\sum f^{(i)}(\alpha)\frac{x^i}{i!}`.  Note
        that the slopes of the Newton polygon of `f(x-\alpha)` provide the
        distances of `\alpha` from its conjugates, i.e., `v(\alpha-\alpha')`.

        To sum this up: we need to choose all the bounds, i.e., the
        `v(a_i-b_i)` such that the `v(a_i-b_i)+v(\beta^i)` exceed
        `\max v(\alpha-\alpha') + \sum v(\alpha-\alpha')`.

        INPUT:

        - ``minpoly`` -- a monic irreducible separable polynomial defined over
          this ring

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: C = Completion(QQ, v)
            sage: R.<x> = C[]
            sage: C._krasner_bound(x^2 + x + 1)
            [0, 0, +Infinity]
            sage: C._krasner_bound(x^2 + 2)
            [3, 5/2, +Infinity]

        """
        if not minpoly.is_monic():
            raise ValueError("minpoly must be monic")
        if minpoly.is_constant():
            raise ValueError("minpoly must not be constant")
        if self.characteristic() != 0:
            if minpoly.gcd(minpoly.derivative()):
                raise ValueError("Krasner's Lemma only applies to separable polynomials")
        if not assume_irreducible and not minpoly.is_irreducible():
            raise ValueError("minpoly must be irreducible")

        f = minpoly
        # determine the coefficients of the Taylor expansion to compute the distances v(alpha-alpha')
        derivatives = [f]
        while derivatives[-1]:
            derivatives.append(derivatives[-1].derivative())
        derivatives.pop()
        # we compute valuations such as v(\alpha) in the formal ring R=K[x]/(f)
        ext = self.extension(f)
        ext_valuation = self.valuation().extension(ext)
        from sage.all import ZZ
        taylor = [d(ext.gen())/ZZ(i).factorial() for i,d in enumerate(derivatives)]
        valuations = [ext_valuation(c) for c in taylor]
        # the distances v(\alpha-\alpha') are the slopes of the Newton polygon of the Taylor expansion
        from sage.geometry.newton_polygon import NewtonPolygon
        distances = NewtonPolygon(enumerate(valuations)).principal_part().slopes()
        distances = [-d for d in distances]

        # For Krasner's lemma to hold, another polynomial g with root beta must satisfy
        # v(alpha-beta) + sum v(alpha'-beta) > max v(alpha-alpha') + sum v(alpha-alpha').
        # The right hand side of this expression is:
        krasner_bound = max(distances + [0]) + sum(distances)

        # First, we need the bounds to be such that the roots of g have the
        # same valuations as the roots of f, i.e., g needs to have the same
        # Newton polygon (which must have only one segment)
        newton_slope = f[0].valuation()/f.degree()
        from sage.all import infinity
        newton_bound = [(f.degree() - i)*newton_slope for i in range(f.degree())] + [infinity]
        # We also need v(a_i-b_i) + v(beta^i) > krasner_bound:
        coefficient_bound = [krasner_bound - i*f[0].valuation()/f.degree() for i in range(f.degree()+1)]
        return [max(*b) for b in zip(newton_bound, coefficient_bound)]


class Completion_Ring(Completion_base):
    r"""
    The completion of ``base`` with respect to ``base_valuation``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: Completion(ZZ, v)
        Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base, base_valuation, category = None):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: C = Completion(ZZ, v)
            sage: isinstance(C, Completion_Ring)
            True
            sage: TestSuite(C).run() # long time
            
        """
        from sage.categories.fields import Fields
        if base in Fields():
            raise TypeError("base must not be a field")

        if category is None:
            from completions import CompleteDiscreteValuationRings
            category = CompleteDiscreteValuationRings()

        super(Completion_Ring, self).__init__(base=base, base_valuation=base_valuation, category=category)

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the fraction field
        of :meth:`base`.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: x = R(0) # indirect doctest
            sage: isinstance(x, R._base_element_class)
            True

        """
        from base_element import BaseElement_Ring
        return self.__make_element_class__(BaseElement_Ring)

    def is_field(self, *args, **kwargs):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: R = Completion(ZZ, v)
            sage: R.fraction_field()
            Completion of Rational Field with respect to 2-adic valuation

        """
        return Completion(self._base_fraction_field, self._base_fraction_field_valuation)


class Completion_Field(Completion_base, Field):
    r"""
    The completion of the field ``base`` with respect to ``base_valuation``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: Completion(QQ, v)
        Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base, base_valuation, category = None):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: isinstance(K, Completion_Field)
            True
            sage: TestSuite(K).run() # long time

        """
        if category is None:
            from completions import CompleteDiscreteValuationFields
            category = CompleteDiscreteValuationFields()

        super(Completion_Field, self).__init__(base=base, base_valuation=base_valuation, category=category)

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the base ring.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: x = K(0) # indirect doctest
            sage: isinstance(x, K._base_element_class)
            True

        """
        from base_element import BaseElement_Field
        return self.__make_element_class__(BaseElement_Field)


class CompletionExtension_base(Completion_base):
    r"""
    Abstract base class for the extension of a completion by adjunction of a
    root of ``polynomial`` to the complete ring ``base_ring``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: L.<a> = K.extension(x^2 + x + 1); L
        Extension defined by a^2 + a + 1 of Completion of Rational Field with respect to 2-adic valuation

        sage: R.<x> = L[]
        sage: M.<b> = L.extension(x^12 - 4*x^11 + 2*x^10 + 13*x^8 - 16*x^7 - 36*x^6 + 168*x^5 - 209*x^4 + 52*x^3 + 26*x^2 + 8*x - 13); M
        Extension defined by b^12 - 4*b^11 + 2*b^10 + 13*b^8 - 16*b^7 - 36*b^6 + 168*b^5 - 209*b^4 + 52*b^3 + 26*b^2 + 8*b - 13 of Extension defined by a^2 + a + 1 of Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, category=None):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: isinstance(L, CompletionExtension_base)
            True
            sage: TestSuite(L).run() # long time

        """
        self._base_ring = base_ring
        self._polynomial = polynomial
        self._name = polynomial.variable_name()
        self._assign_names((self._name,))

        coefficient_ring = base_ring._base
        base_extension_polynomial = polynomial.map_coefficients(coefficient_ring, coefficient_ring)
        from sage.rings.all import QQ
        if base_ring.base() is QQ:
            base_extension = base_ring.base().extension(base_extension_polynomial, names=(self._name,))
        else:
            base_extension = base_extension_polynomial.parent().quo(base_extension_polynomial)
            from sage.categories.all import Fields
            if base_extension in Fields():
                # trigger refinement of category of base_extension
                pass
        base_extension_valuation = base_ring._base_valuation.extension(base_extension)

        super(CompletionExtension_base, self).__init__(base=base_extension, base_valuation=base_extension_valuation, category=category or base_ring.category())

    def degree(self):
        r"""
        Return the degree of this extension over its base.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.degree()
            2
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 + y + L.gen())
            sage: M.degree()
            2

        """
        return self._polynomial.degree()

    def base_ring(self):
        r"""
        Return the base ring of this ring, i.e., the ring from which this
        algebraic extension has been created.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<a> = K[]
            sage: L = K.extension(a^2 + a + 1)
            sage: L.ngens()
            1
            sage: L.gens()
            (a,)

        """
        return 1

    def gen(self, i=0):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<a> = K[]
            sage: L = K.extension(a^2 + a + 1)
            sage: L.gen(0)
            a
            sage: L.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: ring has only one generator

        """
        if i == 0:
            return self(self.base().gen(0))
        raise ValueError("ring has only one generator")

    def _simple_model(self):
        r"""
        Return a simple model of this extension.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<a> = K[]
            sage: L = K.extension(a^2 + a + 1)
            sage: L._simple_model()
            (Number Field in a with defining polynomial a^2 + a + 1, 2-adic valuation)

        """
        base = self.base_ring()
        while base.base() is not base.base_ring():
            base = base.base_ring()
        if self.base_ring() is base:
            return self._base, self._base_valuation
        else:
            return self._eisenstein_model(base=base)
        if self.base_ring().base() is self.base_ring().base_ring():
            return self._base, self._base_valuation
        else:

            return self._simple_eisenstein_model()

    def _eisenstein_model(self, base=None):
        r"""
        Return a model of this extension which is generated by a uniformizing
        element.

        Unless specificed otherwise, the model is constructod over the
        ``base_ring`` of this extension.

        EXAMPLES:

        A totally ramified extension::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: C = Completion(QQ, v)
            sage: R.<t> = C[]
            sage: f = t^12 - 4*t^11 + 2*t^10 + 13*t^8 - 16*t^7 - 36*t^6 + 168*t^5 - 209*t^4 + 52*t^3 + 26*t^2 + 8*t - 13
            sage: C12 = C.extension(f)
            sage: C12._eisenstein_model()
            (Number Field in x with defining polynomial x^12 + 4*x^9 + 12*x^8 + 12*x^5 + 6*x^4 + 4*x^3 + 6, 2-adic valuation)

        A totally ramified extension over an unramified extension::

            sage: R.<u> = C[]
            sage: C2 = C.extension(u^2 + u + 1)
            sage: R.<t> = C2[]
            sage: f = t^12 - 4*t^11 + 2*t^10 + 13*t^8 - 16*t^7 - 36*t^6 + 168*t^5 - 209*t^4 + 52*t^3 + 26*t^2 + 8*t - 13
            sage: C24 = C2.extension(f)
            sage: C24._eisenstein_model(base=C)
            (Number Field in x with defining polynomial x^24 + 4*x^20 + 4*x^17 + 10*x^16 + 8*x^13 + 4*x^12 + 8*x^9 + 12*x^8 + 24*x^7 + 24*x^5 + 28*x^4 + 36, 2-adic valuation)

        """
        if base is None:
            base = self.base_ring()
        if self.gen().valuation() == self.valuation().value_group().gen() and self.base_ring() is base:
            return self._base, self._base_valuation

        prototypical_generator = self.valuation().lift(self.valuation().residue_field().gen()) * self.uniformizer()
        perturbation = 0
        while True:
            generator = prototypical_generator*(1 + self.uniformizer()*perturbation)
            assert generator.valuation() == self.uniformizer().valuation()
            charpoly = generator.matrix(base).charpoly()
            if not charpoly.is_squarefree():
                if self.domain().characteristic() != 0:
                    raise NotImplementedError("iteration over the elements of domains of positive characteristic")
                perturbation += 1
                continue
            error = base._krasner_bound(charpoly, assume_irreducible=True)
            charpoly = charpoly.parent()([c.simplify(e + self.valuation().value_group().gen(), force=True) for c,e in zip(charpoly.coefficients(sparse=False), error)])
            isomorphic_extension = base.extension(charpoly)
            return isomorphic_extension._base, isomorphic_extension._base_valuation
        

class CompletionExtension_Ring(CompletionExtension_base, Completion_Ring):
    r"""
    Extension by adjunction of a root of ``polynomial`` to the complete
    valuation ring ``base_ring``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(ZZ, 2)
        sage: S = Completion(ZZ, v)
        sage: R.<x> = S[]
        sage: S.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Completion of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, category=None):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(ZZ, 2)
            sage: S = Completion(ZZ, v)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: isinstance(T, CompletionExtension_Ring)
            True
            sage: TestSuite(T).run() # long time

        """
        super(CompletionExtension_Ring, self).__init__(base_ring=base_ring, polynomial=polynomial, category=category)


class CompletionExtension_Field(CompletionExtension_base, Completion_Field):
    r"""
    Extension by adjunction of a root of ``polynomial`` to the complete field
    ``base_ring``.

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: R.<x> = K[]
        sage: K.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Completion of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, category=None):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: isinstance(L, CompletionExtension_Field)
            True
            sage: TestSuite(L).run() # long time

        """
        super(CompletionExtension_Field, self).__init__(base_ring=base_ring, polynomial=polynomial, category=category)

    def _vector_space_basis(self, base):
        r"""
        Return a basis of this field as a vector space over ``base``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L._vector_space_basis(L)
            (1,)
            sage: L._vector_space_basis(K)
            (1, x)

        """
        if base is self:
            return (self.one(),)
        basis = tuple(self.gen()**i for i in range(self.degree()))
        return tuple(b*self(p) for b in basis for p in self.base_ring()._vector_space_basis(base))
