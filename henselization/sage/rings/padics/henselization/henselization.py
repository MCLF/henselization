# -*- coding: utf-8 -*-
r"""
Henselizations of rings
=======================

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
from sage.rings.ring import CommutativeRing, Field
from sage.structure.factory import UniqueFactory
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

from sage.categories.fields import Fields

class HenselizationFactory(UniqueFactory):
    r"""
    Return the Henselization of ``R`` with respect to ``v``.
    
    INPUT:
    
    - ``R`` -- a field or an excellent integral domain whose localization at
      the maximal ideal of ``v`` is the discrete valuation ring of ``v``. (Most
      of these conditions are not verified by this factory.)

    - ``v`` -- a non-trivial discrete valuation on ``R``.

    EXAMPLES::

        sage: from henselization import *
        sage: QQ.henselization(5) # indirect doctest
        Henselization of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, R, v):
        r"""
        Create a key which uniquely identifies this Henselization.

        TESTS::

            sage: from henselization import *
            sage: QQ.henselization(5) is QQ.henselization(5) # indirect doctest
            True

        """
        from sage.categories.all import IntegralDomains
        if R not in IntegralDomains():
            raise ValueError("R must be an integral domain")
        from sage.rings.valuation.valuation_space import DiscretePseudoValuationSpace
        if v not in DiscretePseudoValuationSpace(R) or not v.is_discrete_valuation():
            raise ValueError("v must be a discrete valuation on R")
        if v.is_trivial():
            raise ValueError("v must not be trivial")

        return R, v

    def create_object(self, version, key):
        r"""
        Create the Henselization identified by ``key``.

        TESTS::
            
            sage: from henselization import *
            sage: Henselization(QQ, QQ.valuation(5)) # indirect doctest
            Henselization of Rational Field with respect to 5-adic valuation

        """
        R, v = key

        if v.value_semigroup().is_group():
            return Henselization_Field(R, v)
        else:
            return Henselization_Ring(R, v)

Henselization = HenselizationFactory("Henselization")


class ExtensionFactory(UniqueFactory):
    r"""
    Return the algebraic extension which adjoins to ``base`` a root of
    ``polynomial``.

    Do not call this factory directly, but call ``base.extension()`` instead.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(5)
        sage: R.<x> = K[]
        sage: K.extension(x^2 - 5) # indirect doctest
        Extension defined by x^2 - 5 of Henselization of Rational Field with respect to 5-adic valuation

    """
    def create_key(self, base, polynomial, name, check=True):
        r"""
        Return a key that uniquely defines this extension.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
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

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: K.extension(x^2 - 5) # indirect doctest
            Extension defined by x^2 - 5 of Henselization of Rational Field with respect to 5-adic valuation

        """
        base, polynomial = key

        if isinstance(base, HenselizationExtension):
            if polynomial.base_ring()._is_monic_mac_lane_polynomial(polynomial) and polynomial.base_ring() is base:
                iterated_extension = self._get_isomorphic_approximation(polynomial)
            else:
                iterated_extension = Extension._create_extension(base, polynomial)
            model, model_valuation = iterated_extension._eisenstein_model(iterated_extension._absolute_base_ring())
            return Extension._create_extension(base, polynomial, model, model_valuation)
        else:
            return Extension._create_extension(base, polynomial)

    def _get_isomorphic_approximation(self, polynomial):
        r"""
        Return a Henselization backed by a quotient by a polynomial that has
        the same splitting field as the approximate ``polynomial`` over
        ``base``.

        EXAMPLES:

        A trivial case, a factor of degree one::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: F = (x^2 + 1).factor()
            sage: F = sorted(F, key=str) # remove randomness in the output
            sage: F[0][0], F[1][0]
            (x + 2 + O(5), x + 3 + O(5))
            sage: from sage.rings.padics.henselization.henselization import Extension
            sage: Extension._get_isomorphic_approximation(F[0][0])._base
            Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x + 2

        Non-trivial factors::

            sage: F = (x^6 + 3).factor()
            sage: F = sorted(F, key=str) # remove randomness in the output
            sage: F[0][0], F[1][0], F[2][0]
            (x^2 + (1 + O(5))*x + 2 + O(5),
             x^2 + (4 + O(5))*x + 2 + O(5),
             x^2 + O(5)*x + 2 + O(5))
            sage: Extension._get_isomorphic_approximation(F[0][0])._base
            Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + x + 2
            sage: Extension._get_isomorphic_approximation(F[1][0])._base
            Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + 4*x + 2
            sage: Extension._get_isomorphic_approximation(F[2][0])._base
            Univariate Quotient Polynomial Ring in xbar over Rational Field with modulus x^2 + 2

        A complex case where the initial approximation is not sufficient::

            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: f = x^12 - 4*x^11 + 2*x^10 + 13*x^8 - 16*x^7 - 36*x^6 + 168*x^5 - 209*x^4 + 52*x^3 + 26*x^2 + 8*x - 13
            sage: L = K.extension(f.change_variable_name('a'))
            sage: F = f.change_ring(L).factor()
            sage: F = sorted(F, key=str) # remove randomness in the output
            sage: g = F[3][0] # a factor of degree 8

        Here, the initial approximation of the factor is not sufficient to
        single out the right splitting field::

            sage: F[3][0][0]
            57/2*a^11 + 30*a^10 - 3/2*a^9 + 35/2*a^8 + 17*a^7 + 13*a^6 + a^5 + 19*a^4 + 33/2*a^3 + 9*a^2 + 37/2*a - 5/6 + O(...)
            sage: Extension._get_isomorphic_approximation(F[3][0])._base
            Univariate Quotient Polynomial Ring in xbar over Number Field in a with defining polynomial a^12 - 4*a^11 + 2*a^10 + 13*a^8 - 16*a^7 - 36*a^6 + 168*a^5 - 209*a^4 + 52*a^3 + 26*a^2 + 8*a - 13 with modulus x^8 + (5/6*a^11 + 9/2*a^10 + 3/2*a^9 + 1/10*a^8 + 5*a^7 + 5*a^6 + 9*a^5 + 7*a^4 + 19/2*a^3 + 17/2*a^2 + 15/2*a + 1/2)*x^7 + (53/2*a^11 - 1/6*a^10 + 19/2*a^9 + 15*a^7 + 8*a^6 + 17*a^5 + 5*a^4 + 9/2*a^3 + 35/2*a^2 + 47/2*a + 9)*x^6 + (15*a^11 + 29*a^10 + 26*a^9 + 11*a^8 + 28*a^7 + 22*a^6 + 10*a^5 + 28*a^4 + 9*a^3 + 21*a^2 + 12*a + 21)*x^5 + (31/2*a^11 + 5*a^10 + 55/2*a^9 + 15/2*a^8 + 13*a^7 + a^6 + 29*a^5 + 31*a^4 - 1/10*a^3 + 26*a^2 - 1/10*a - 1/2)*x^4 + (29/2*a^11 + 7/2*a^10 + 49/2*a^9 + 7/2*a^8 + 17*a^7 + 25*a^6 + 17*a^5 + 31*a^4 + 37/2*a^3 + 11/2*a^2 + 17/2*a + 7/2)*x^3 + (27/2*a^11 + 3/2*a^10 - 3/10*a^9 + 23*a^8 + 25*a^7 + 22*a^6 + 11*a^5 + 13*a^4 + 55/2*a^3 - 5/6*a^2 - 3/2*a + 22)*x^2 + (12*a^11 + 11*a^9 + 16*a^8 + 4*a^7 + 10*a^6 + 14*a^5 + 8*a^4 + 20*a^3 + 22*a^2 + 15*a + 12)*x + 33/2*a^11 + 28*a^10 + 57/2*a^9 + 11/2*a^8 + 3*a^7 + 23*a^6 + 23*a^5 + 9*a^4 - 3/2*a^3 + 29*a^2 + 37/2*a + 45/2

        """
        coefficient_ring = polynomial.base_ring()._base

        if polynomial.degree() < 1:
            raise ValueError("constants have no splitting field")

        limit = polynomial[0]._limit_valuation
        G = limit._G.change_ring(coefficient_ring)
        # We extend coefficient_ring by a sufficiently precise
        # approximation of the polynomial that limit approaches.
        # We need a root of g = limit._approximation.phi() to be closer to
        # a root of G than to any other root of g (Krasner's Lemma.)
        while True:
            g = limit._approximation.phi().change_ring(coefficient_ring)

            # Let z be a root of g and consider the Newton polygons of
            # G(T+z) and g(T+z)/T in the quotient model = g.parent() mod (g)
            I = g.parent().ideal(g)
            if hasattr(I.gen().is_irreducible, 'set_cache'):
                I.gen().is_irreducible.set_cache(True)
            model = g.parent().quo(I)
            z = model.gen()
            R = model['T']; T = R.gen()

            # We could now determine the valuation on model by running
            # model_valuation = polynomial.base_ring().extension(model)
            # However, this is not necessary. We already know what this
            # valuation is going to look like, namely it is
            # limit._approximation with the last step set to v(g)=infinity.
            from sage.all import infinity
            model_valuation = limit._approximation._base_valuation.change_domain(g.parent()).augmentation(g, infinity, check=False).change_domain(model)

            from sage.rings.valuation.all import GaussValuation
            w = GaussValuation(R, model_valuation)

            if g.degree() == 1:
                break

            g_shift = g(T + z)
            G_shift = G(T + z)
            # Consider the following Newton polygons:
            # GNP = w.newton_polygon(G_shift)
            # gNP = w.newton_polygon(g_shift.shift(-1))
            #
            # If the first slope of GNP is more negative than the first
            # slope of gNP, then z is closer to a root of G than it is to
            # another root of g.
            # if (GNP.slopes()[0] < gNP.slopes()[0]
            #   # or z is an actual root of G
            #   or GNP.vertices()[0][0] > 1):
            #     return g

            # The above idea can be implemented more efficiently:
            # Let ζ be the root of G closest to z, the root of g.
            # Since g is an approximate factor of G, z approaches ζ as we
            # perform limit._improve_approximation(), i.e., as we improve the
            # quality of the approximate factor g. In other words, the absolute
            # value of the first slope of the Newton polygon of  G(T + z) goes
            # to infinity (but the other slopes don't).
            # Therefore, the difference in valuation of the constant
            # coefficient and the linear coefficient of G(T + z) eventually is
            # the first slope of the Newton polygon of G(T + z)
            # If it is not the first slope, it can only be smaller in absolute
            # value, so we have a bound that becomes eventually an equality:
            # GNP.slopes()[0] ≤ w(G(T + z)[1]) - w(G(T + z)[0])
            gNP = w.newton_polygon(g_shift.shift(-1))
            if model_valuation(G_shift[1]) - model_valuation(G_shift[0]) < gNP.slopes()[0]:
                break

            limit._improve_approximation()

        return Extension._create_extension(polynomial.base_ring(), polynomial, model=model, model_valuation=model_valuation)


    def _create_extension(self, base, polynomial, model = None, model_valuation = None):
        r"""
        Return the extension of ``base`` by adjoining a root of ``polynomial``.
    
        TESTS:
    
        A ``model`` and a ``model_valuation`` can be specified to specify which
        number field should back the implementation of the Henselization; if the
        ``model_valuation`` is omitted it is determined automatically::
    
            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = QQ[]
            sage: M.<x> = QQ.extension(x^2 + x + 3)
            sage: R.<x> = K[]
            sage: from sage.rings.padics.henselization.henselization import Extension
            sage: L = Extension._create_extension(K, x^2 + x + 1, model = M)
            sage: L
            Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
            sage: L.base()
            Number Field in x with defining polynomial x^2 + x + 3
    
        If no ``model`` is specified, then an iterated extension is backed by a
        quotient that just divides out ``polynomial``::
    
            sage: R.<y> = L[]
            sage: M = Extension._create_extension(L, y^2 - 2)
            sage: M.base()
            Univariate Quotient Polynomial Ring in ybar over Number Field in x with defining polynomial x^2 + x + 3 with modulus y^2 - 2
    
        """
        return Quotient(base, polynomial, model, model_valuation)

class QuotientFactory(UniqueFactory):
    def create_key(self, base, polynomial, model = None, model_valuation = None):
        if model is None:
            from base_element import BaseElement_base
            if all([isinstance(c, BaseElement_base) for c in polynomial.coefficients(sparse=False)]):
                model_polynomial = polynomial.map_coefficients(base._base, base._base)
            else:
                model_polynomial = Extension._get_isomorphic_approximation(polynomial)._base.modulus()

            if isinstance(base, HenselizationExtension):
                model = model_polynomial.parent().quo(model_polynomial)
            else:
                model = base.base().extension(model_polynomial, names=(model_polynomial.variable_name(),))
        if model_valuation is None:
            model_valuation = base._base_valuation.extension(model)

        return base, polynomial, model, model_valuation

    def create_object(self, version, key):
        base, polynomial, model, model_valuation = key

        from sage.categories.all import Fields
        if model in Fields():
            # triggers refinement of category of model
            pass
    
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
        if not isinstance(base, HenselizationExtension):
            if base in Fields():
                clazz = HenselizationExtensionSimple_Field
            else:
                clazz = HenselizationExtensionSimple_Ring
        elif is_PolynomialQuotientRing(model):
            if base in Fields():
                clazz = HenselizationExtensionIteratedQuotient_Field
            else:
                clazz = HenselizationExtensionIteratedQuotient_Ring
        else:
            if base in Fields():
                clazz = HenselizationExtensionIteratedAbsolute_Field
            else:
                clazz = HenselizationExtensionIteratedAbsolute_Ring
        
        return clazz(base, polynomial, model, model_valuation)

Extension = ExtensionFactory("sage.rings.padics.henselization.henselization.Extension")
Quotient = QuotientFactory("sage.rings.padics.henselization.henselization.Quotient")


class Henselization_base(CommutativeRing):
    r"""
    Abstract base class for the Henselization of ``base`` with respect to
    ``base_valuation``.

    EXAMPLES::

        sage: from henselization import *
        sage: Henselization(ZZ, ZZ.valuation(2))
        Henselization of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base, base_valuation, category):
        r"""
        TESTS::

            sage: from henselization import *
            sage: K.<x> = FunctionField(QQ)
            sage: C = K.henselization(x)
            sage: from sage.rings.padics.henselization.henselization import Henselization_base
            sage: isinstance(C, Henselization_base)
            True
            sage: TestSuite(C).run() # long time

        """
        super(Henselization_base, self).__init__(base_ring=base, category=category)

        self._base_valuation = base_valuation

        # The Henselization contains not only the elements of base but many more
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

        # provide conversions from the Henselization back to its base and its field of fractions
        from maps import ConvertMap_generic
        from sage.categories.all import SetsWithPartialMaps
        # TODO: use weak references
        homspace = self.Hom(self._base, category=SetsWithPartialMaps())
        self._base.register_conversion(homspace.__make_element_class__(ConvertMap_generic)(homspace))
        if self._base is not self._base_fraction_field:
            homspace = self.Hom(self._base_fraction_field, category=SetsWithPartialMaps())
            self._base_fraction_field.register_conversion(homspace.__make_element_class__(ConvertMap_generic)(homspace))

    def _is_squarefree_univariate_polynomial(self, f):
        r"""
        Return whether ``f`` is squarefree.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: F = (x^2 + 4).factor()
            sage: F[0][0].is_squarefree() # indirect doctest
            True
            sage: F[1][0].is_squarefree()
            True

        """
        if self._is_monic_mac_lane_polynomial(f):
            return True
        return f._is_squarefree_generic()

    def _gcd_univariate_polynomial(self, f, g):
        r"""
        Return the greatest common divisor of ``f`` and ``g``.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: (x^2 - 1).gcd(x - 1) # indirect doctest
            x - 1

        """
        # Since the gcd is independent of the base, we can compute it over
        # our base ring which is usually much faster (using a multimodular
        # algorithm for example) than the naive implementation
        return f.change_ring(self._base_fraction_field).gcd(g.change_ring(self._base_fraction_field)).change_ring(self)

    def _xgcd_univariate_polynomial(self, f, g):
        r"""
        Return the xgcd of ``f`` and  ``g``.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

        There are two different base rings for Henselizations. The :meth:`base`
        which is the ring that the Henselization completes and the
        :meth:`base_ring` which is the ring from which the Henselization has been
        constructed.

        EXAMPLES:

        For Henselizations that are not algebraical extensions of another
        Henselization, these two concepts coincide::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K.base_ring() is K.base() is QQ
            True

        """
        return self.base()

    def characteristic(self):
        r"""
        Return the characteristic of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K.characteristic()
            0

        """
        return self.base().characteristic()

    def is_finite(self):
        r"""
        Return whether this ring is finite.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K.uniformizer()
            2

        """
        return self(self._base_valuation.uniformizer())

    def _coerce_map_from_(self, other):
        r"""
        Return a coercion from ``other`` to this ring if one exists.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K.has_coerce_map_from(QQ) # indirect doctest
            True

        """
        if self.base().has_coerce_map_from(other):
            return True
        
        from maps import ExtensionCoercion_generic
        if self.base_ring().has_coerce_map_from(other):
            if isinstance(other, HenselizationExtension):
                if self.base().has_coerce_map_from(other.base()):
                    homspace = other.Hom(self)
                    return homspace.__make_element_class__(ExtensionCoercion_generic)(homspace)

        if isinstance(other, Henselization_base):
            if other.base().is_subring(self.base()):
                if self._base_valuation.restriction(other.base()) == other._base_valuation:
                    homspace = other.Hom(self)
                    return homspace.__make_element_class__(ExtensionCoercion_generic)(homspace)
        return super(Henselization_base, self)._coerce_map_from_(other)

    def _an_element_(self):
        r"""
        Return an element of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K.an_element()
            2

        """
        return self.uniformizer()

    def _element_constructor_(self, x):
        r"""
        Create an element in this parent from ``x``.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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
                from sage.rings.valuation.limit_valuation import MacLaneLimitValuation
                if isinstance(v, MacLaneLimitValuation) and v.domain().base() is self and i in NN:
                    return self._mac_lane_element_class(self, v, i)

        raise ValueError("Can not convert %r to an element in %r"%(x, self))

    def _repr_(self):
        r"""
        Return a printable representation of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: Henselization(QQ, QQ.valuation(2))
            Henselization of Rational Field with respect to 2-adic valuation

        """
        return "Henselization of %r with respect to %r"%(self.base_ring(), self._base_valuation)

    def some_elements(self):
        r"""
        Return some typical elements of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K.residue_field()
            Finite Field of size 2

        """
        return self.valuation().residue_field()

    def extension(self, f, names=None, name=None, check=True):
        r"""
        Return the algebraic extension of this ring obtained by adjoining a
        root of the irreducible polynomial ``f``.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: K.extension(x^2 + x + 1)
            Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: f = x + 1
            sage: f.factor() # indirect doctest
            x + 1
            sage: f = x^2 + 1
            sage: f.factor()
            x^2 + 1

        A non-trivial example::

            sage: G = GaussianIntegers().fraction_field()
            sage: I = G.gen()
            sage: K = G.henselization(2)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: F = f.factor(); F # random output
            (x + 3*I + O(?)) * (x + I)
            sage: (x + I, 1) in F
            True
            sage: (x - I, 1) in F
            True

        Another non-trivial example::

            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: f.factor()
            (x + 2 + O(5)) * (x + 3 + O(5))

        """
        if f.is_constant():
            raise NotImplementedError("factorization of constant polynomials")
        if not f.is_monic():
            raise NotImplementedError("factorization of non-monic polynomials")
        if any([c.valuation() < 0 for c in f.coefficients()]):
            raise NotImplementedError("factorization of non-integral polynomials")
        if self._is_monic_mac_lane_polynomial(f):
            from sage.structure.factorization import Factorization
            return Factorization([(f,1)])
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
        # We have to require_incomparability so we can unique improve the
        # resulting approximations to factors
        approximants = self.valuation().mac_lane_approximants(f, require_incomparability=True, require_maximal_degree=True)
        if len(approximants) == 1:
            return Factorization([(f, 1)], sort=False)
        factors = []
        for approximant in approximants:
            from sage.rings.all import infinity
            if approximant(approximant.phi()) == infinity:
                factors.append(approximant.phi())
                continue

            degree = approximant.phi().degree()
            from sage.rings.valuation.limit_valuation import LimitValuation
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

            sage: from henselization import *
            sage: R = ZZ.henselization(2)
            sage: R.ideal(1, 1)
            Principal ideal (1) of Henselization of Integer Ring with respect to 2-adic valuation

        """
        I = super(Henselization_base, self).ideal(*args, **kwds)
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
        return super(Henselization_base, self).ideal(gens)

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

            sage: from henselization import *
            sage: C = QQ.henselization(2)
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
        # we compute valuations such as v(\alpha) in the formal ring R=K[x]/(f)
        ext = self.base()[f.variable_name()].quo(f.change_ring(self.base()))
        ext_valuation = self._base_valuation.extension(ext)

        from sage.all import ZZ
        taylor = [(ext(d),ZZ(i).factorial()) for i,d in enumerate(derivatives) if i != 0]

        valuations = [ext_valuation(n)-self._base_valuation(d) for n,d in taylor]
        # the distances v(\alpha-\alpha') are the slopes of the Newton polygon of the Taylor expansion
        from sage.geometry.newton_polygon import NewtonPolygon
        distances = NewtonPolygon(list(enumerate(valuations))).slopes()
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

    def _is_monic_mac_lane_polynomial(self, polynomial):
        r"""
        Return whether ``polynomial`` comes from an approximate factorization
        of a monic polynomial over this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: F = (x^2 + 4).factor()
            sage: K._is_monic_mac_lane_polynomial(F[0][0])
            True
            sage: K._is_monic_mac_lane_polynomial(x)
            False

        """
        if polynomial.base_ring() is not self:
            return False
        if polynomial.degree() < 1:
            return False
        if not polynomial.is_monic():
            return False
        from mac_lane_element import MacLaneElement_base
        return polynomial.is_monic() and all([isinstance(c, MacLaneElement_base)
                                              and c._degree == d
                                              and c._limit_valuation == polynomial[0]._limit_valuation
                                              for (d,c) in enumerate(polynomial.coefficients(sparse=False)[:-1])])

    def module(self, base=None):
        r"""
        Return a free ``base``-module isomorphic to this ring together with
        isomorphisms to and from this module.

        INPUT:

        - ``base`` -- a ring of which this ring is a finite extension
          (default: the ring itself)

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.module(base=L)
            (Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
               To:   Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
               To:   Vector space of dimension 1 over Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation)
            sage: L.module(base=K)
            (Vector space of dimension 2 over Henselization of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Vector space of dimension 2 over Henselization of Rational Field with respect to 2-adic valuation
               To:   Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation,
             Isomorphism morphism:
               From: Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation
               To:   Vector space of dimension 2 over Henselization of Rational Field with respect to 2-adic valuation)

        """
        if base is None:
            base = self
        basis = self._module_basis(base=base)
        V = base**len(basis)
        from sage.all import Hom
        to_self_parent = Hom(V, self)
        from maps import VectorSpaceToHenselization, HenselizationToVectorSpace
        to_self = to_self_parent.__make_element_class__(VectorSpaceToHenselization)(to_self_parent, basis)
        from_self_parent = Hom(self, V)
        from_self = from_self_parent.__make_element_class__(HenselizationToVectorSpace)(from_self_parent, base)
        return (V, to_self, from_self)

    def _module_basis(self, base):
        r"""
        Return a basis of this ring as a free module over ``base``.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K._module_basis(K)
            (1,)

        """
        if base is None:
            base = self
        if base is self:
            return (self(1),)
        raise NotImplementedError("Basis with respect to %s"%base)

    def _test_module_basis(self, **options):
        r"""
        Test that :meth:`_module_basis` works correctly.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: K._test_module_basis()

        """
        tester = self._tester(**options)

        basis = self._module_basis(self)
        tester.assertTrue(all(c.parent() is self for c in basis))
        tester.assertEqual(len(basis), 1)


class Henselization_Ring(Henselization_base):
    r"""
    The Henselization of ``base`` with respect to ``base_valuation``.

    EXAMPLES::

        sage: from henselization import *
        sage: Henselization(ZZ, ZZ.valuation(2))
        Henselization of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base, base_valuation, category = None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: C = ZZ.henselization(2)
            sage: from sage.rings.padics.henselization.henselization import Henselization_Ring
            sage: isinstance(C, Henselization_Ring)
            True
            sage: TestSuite(C).run() # long time
            
        """
        from sage.categories.fields import Fields
        if base in Fields():
            raise TypeError("base must not be a field")

        if category is None:
            from sage.categories.discrete_valuation import HenselianDiscreteValuationRings
            category = HenselianDiscreteValuationRings()

        super(Henselization_Ring, self).__init__(base=base, base_valuation=base_valuation, category=category)

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the fraction field
        of :meth:`base`.

        TESTS::

            sage: from henselization import *
            sage: R = ZZ.henselization(2)
            sage: x = R(0) # indirect doctest
            sage: isinstance(x, R._base_element_class)
            True

        """
        from base_element import BaseElement_Ring
        return self.__make_element_class__(BaseElement_Ring)

    @lazy_attribute
    def _mac_lane_element_class(self):
        r"""
        The class for elements which are given as coefficients of key
        polynomial of limit valuations.

        TESTS::

            sage: from henselization import *
            sage: S = ZZ.henselization(5)
            sage: R.<x> = S[]
            sage: f = x^2 + 1
            sage: F = f.factor() # indirect doctest
            sage: isinstance(F[0][0][0], S._mac_lane_element_class)
            True

        """
        from mac_lane_element import MacLaneElement_Ring
        return self.__make_element_class__(MacLaneElement_Ring)

    def is_field(self, *args, **kwargs):
        r"""
        Return whether this ring is a field.

        EXAMPLES::

            sage: from henselization import *
            sage: R = ZZ.henselization(2)
            sage: R.is_field()
            False

        """
        return False

    def fraction_field(self):
        r"""
        Return the fraction field of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: R = ZZ.henselization(2)
            sage: R.fraction_field()
            Henselization of Rational Field with respect to 2-adic valuation

        """
        return Henselization(self._base_fraction_field, self._base_fraction_field_valuation)


class Henselization_Field(Henselization_base, Field):
    r"""
    The Henselization of the field ``base`` with respect to ``base_valuation``.

    EXAMPLES::

        sage: from henselization import *
        sage: Henselization(QQ, QQ.valuation(2))
        Henselization of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base, base_valuation, category = None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: from sage.rings.padics.henselization.henselization import Henselization_Field
            sage: isinstance(K, Henselization_Field)
            True
            sage: TestSuite(K).run() # long time

        """
        if category is None:
            from sage.categories.discrete_valuation import HenselianDiscreteValuationFields
            category = HenselianDiscreteValuationFields()

        super(Henselization_Field, self).__init__(base=base, base_valuation=base_valuation, category=category)

    @lazy_attribute
    def _base_element_class(self):
        r"""
        The class for elements which are already elements of the base ring.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: x = K(0) # indirect doctest
            sage: isinstance(x, K._base_element_class)
            True

        """
        from base_element import BaseElement_Field
        return self.__make_element_class__(BaseElement_Field)

    @lazy_attribute
    def _mac_lane_element_class(self):
        r"""
        The class for elements which are given as coefficients of key
        polynomial of limit valuations.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(5)
            sage: R.<x> = K[]
            sage: f = x^2 + 1
            sage: F = f.factor() # indirect doctest
            sage: isinstance(F[0][0][0], K._mac_lane_element_class)
            True

        """
        from mac_lane_element import MacLaneElement_Field
        return self.__make_element_class__(MacLaneElement_Field)


class HenselizationExtension(Henselization_base):
    r"""
    Abstract base class for the extension of a Henselization by adjunction of a
    root of ``polynomial`` to the complete ring ``base_ring``.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: R.<x> = K[]
        sage: L.<a> = K.extension(x^2 + x + 1); L
        Extension defined by a^2 + a + 1 of Henselization of Rational Field with respect to 2-adic valuation

        sage: R.<x> = L[]
        sage: M.<b> = L.extension(x^12 - 4*x^11 + 2*x^10 + 13*x^8 - 16*x^7 - 36*x^6 + 168*x^5 - 209*x^4 + 52*x^3 + 26*x^2 + 8*x - 13); M
        Extension defined by b^12 - 4*b^11 + 2*b^10 + 13*b^8 - 16*b^7 - 36*b^6 + 168*b^5 - 209*b^4 + 52*b^3 + 26*b^2 + 8*b - 13 of Extension defined by a^2 + a + 1 of Henselization of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtension
            sage: isinstance(L, HenselizationExtension)
            True
            sage: TestSuite(L).run() # long time

        We can handle defining polynomials whose coefficients are Mac Lane
        elements::

            sage: R.<T> = K[]
            sage: f = T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13
            sage: L.<a12> = K.extension(f)
            sage: F = f.change_ring(L).factor()
            sage: g = [g for (g,e) in F if g.degree() == 2][0]
            sage: M.<a24> = L.extension(g); M
            Extension defined by a24^2 + (5*a12^11 + 15/2*a12^10 + 7*a12^9 + 27/4*a12^8 + 6*a12^7 + 7/2*a12^4 + 5*a12^3 + 13/2*a12^2 + 3*a12 + 19/4 + O((1/4*a12^7 + 3/4*a12^6 + 5/4*a12^5 + 7/4*a12^4 + 9/4*a12^3 + 11/4*a12^2 + 9/4*a12 + 3/4)^(5/6)))*a24 + 15/2*a12^11 + 1/2*a12^10 + 11/4*a12^9 + 2*a12^8 + 6*a12^7 + 2*a12^6 + 15/2*a12^5 + 2*a12^4 + 1/2*a12^3 + 1/2*a12^2 + 19/4*a12 + 2 + O((1/4*a12^7 + 3/4*a12^6 + 5/4*a12^5 + 7/4*a12^4 + 9/4*a12^3 + 11/4*a12^2 + 9/4*a12 + 3/4)^(5/6)) of Extension defined by a12^12 - 4*a12^11 + 2*a12^10 + 13*a12^8 - 16*a12^7 - 36*a12^6 + 168*a12^5 - 209*a12^4 + 52*a12^3 + 26*a12^2 + 8*a12 - 13 of Henselization of Rational Field with respect to 2-adic valuation
    
        """
        if not isinstance(base_ring, Henselization_base):
            raise TypeError("base_ring must be a Henselization")
        if polynomial.base_ring() is not base_ring:
            raise ValueError("polynomial must be defined over base_ring")
        if model_valuation.domain() is not model:
            raise ValueError("model_valuation must be defined on model")

        self._base_ring = base_ring
        self._polynomial = polynomial
        self._model = model
        self._name = polynomial.variable_name()
        self._assign_names((self._name,))

        super(HenselizationExtension, self).__init__(base=model, base_valuation=model_valuation, category=category or base_ring.category())

        from sage.categories.all import SetsWithPartialMaps
        homspace = base_ring.Hom(self, category=SetsWithPartialMaps())
        from maps import RelativeExtensionCoercion_generic
        self.register_coercion(homspace.__make_element_class__(RelativeExtensionCoercion_generic)(homspace))

    def degree(self):
        r"""
        Return the degree of this extension over its base.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
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

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: L.base_ring() is K
            True

        """
        return self._base_ring

    def _absolute_base_ring(self):
        r"""
        Return the :class:`Henselization_base` that lies at the bottom of this
        tower of complete extensions.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 - 2)
            sage: M._absolute_base_ring() is K
            True

        """
        ret = self.base_ring()
        while isinstance(ret, HenselizationExtension):
            ret = ret.base_ring()
        return ret

    def _absolute_degree(self):
        r"""
        Return the degree this extension has over the :class:`Henselization_base`
        that lies at the bottom of this tower of complete extensions.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 - 2)
            sage: M._absolute_degree()
            4

        """
        from sage.all import ZZ
        ret = ZZ(1)
        base = self
        while isinstance(base, HenselizationExtension):
            ret *= base.degree()
            base = base.base_ring()
        return ret

    @cached_method
    def _repr_(self):
        r"""
        Return a printable representation of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: K.extension(x^2 + x + 1) # indirect doctest
            Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

        """
        return "Extension defined by %r of %r"%(self._polynomial, self.base_ring())

    def ngens(self):
        r"""
        Return the number of generators of this ring.

        EXAMPLES:

        This extension is generated by a root of its defining polynomial::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<a> = K[]
            sage: L = K.extension(a^2 + a + 1)
            sage: L.ngens()
            1
            sage: L.gens()
            (a,)

        """
        return 1

    def _eisenstein_model(self, base, name=None):
        r"""
        Return a model of this extension which is simple over``base.base()``
        and generated by a uniformizing element.

        EXAMPLES:

        A totally ramified extension::

            sage: from henselization import *
            sage: C = QQ.henselization(2)
            sage: R.<t> = C[]
            sage: f = t^12 - 4*t^11 + 2*t^10 + 13*t^8 - 16*t^7 - 36*t^6 + 168*t^5 - 209*t^4 + 52*t^3 + 26*t^2 + 8*t - 13
            sage: C12 = C.extension(f)
            sage: C12._eisenstein_model(base=C)
            (Number Field in e12 with defining polynomial e12^12 + 4*e12^9 + 12*e12^8 + 12*e12^5 + 6*e12^4 + 4*e12^3 + 6, 2-adic valuation)

        A totally ramified extension over an unramified extension::

            sage: R.<u> = C[]
            sage: C2 = C.extension(u^2 + u + 1)
            sage: R.<t> = C2[]
            sage: f = t^12 - 4*t^11 + 2*t^10 + 13*t^8 - 16*t^7 - 36*t^6 + 168*t^5 - 209*t^4 + 52*t^3 + 26*t^2 + 8*t - 13
            sage: from sage.rings.padics.henselization.henselization import Extension
            sage: C24 = Extension._create_extension(C2, f)
            sage: C24._eisenstein_model(base=C)
            (Number Field in e12 with defining polynomial e12^24 + 4*e12^20 + 4*e12^17 + 10*e12^16 + 8*e12^13 + 4*e12^12 + 8*e12^9 + 12*e12^8 + 24*e12^7 + 24*e12^5 + 28*e12^4 + 36, 2-adic valuation)

        """
        if base is self.base_ring() and self.gen().valuation() == self.valuation().value_group().gen():
            if name is None or name == self.variable_name():
                return self._base, self._base_valuation
            else:
                base = self._base.change_variable_name(name)
                base_valuation = self._base_valuation.change_ring(base)
                return base, base_valuation

        if name is None:
            name = 'e%s'%(self.valuation().value_group().index(base.valuation().value_group()),)

        prototypical_generator = self.valuation().lift(self.valuation().residue_field().gen()) * self.uniformizer()
        perturbation = 0
        while True:
            generator = prototypical_generator*(1 + self.uniformizer()*perturbation)
            assert generator.valuation() == self.uniformizer().valuation()
            M = generator.matrix(base)
            M = M.apply_map(base.base().convert_map_from(base))
            charpoly = M.charpoly()
            charpoly = charpoly.change_ring(base)
            charpoly = charpoly.change_variable_name(name)
            if not charpoly.is_squarefree():
                if self.domain().characteristic() != 0:
                    raise NotImplementedError("iteration over the elements of domains of positive characteristic")
                perturbation += 1
                continue
            error = base._krasner_bound(charpoly, assume_irreducible=True)
            charpoly = charpoly.parent()([c.simplify(e + self.valuation().value_group().gen(), force=True) for c,e in zip(charpoly.coefficients(sparse=False), error)])
            eisenstein_model = base.base().extension(charpoly.change_ring(base.base()), names=(name,), check=False)
            model_valuation = base._base_valuation.extension(eisenstein_model)
            return eisenstein_model, model_valuation


class HenselizationExtension_Field(HenselizationExtension, Henselization_Field):
    r"""
    A :class:`HenselizationExtension` that is a field.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(3)
        sage: R.<x> = K[]
        sage: K.extension(x^2 + 3)
        Extension defined by x^2 + 3 of Henselization of Rational Field with respect to 3-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(3)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + 3)
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtension_Field
            sage: isinstance(L, HenselizationExtension_Field)
            True
            sage: TestSuite(L).run() # long time

        """
        if model not in Fields():
            raise TypeError("model must be a field")
        super(HenselizationExtension_Field, self).__init__(base_ring=base_ring, polynomial=polynomial, model=model, model_valuation=model_valuation, category=category)

    def _test_extension_module_basis(self, **options):
        r"""
        Test that :meth:`_module_basis` works correctly.

        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + 2)
            sage: L._test_extension_module_basis()

        """
        tester = self._tester(**options)

        basis = self._module_basis(self._absolute_base_ring())
        tester.assertTrue(all(c.parent() is self for c in basis))
        tester.assertEqual(len(basis), self._absolute_degree())


class HenselizationExtension_Ring(HenselizationExtension, Henselization_Ring):
    r"""
    A :class:`HenselizationExtension` that is not a field.

    EXAMPLES::

        sage: from henselization import *
        sage: S = ZZ.henselization(3)
        sage: R.<x> = S[]
        sage: T = S.extension(x^2 + 3); T
        Extension defined by x^2 + 3 of Henselization of Integer Ring with respect to 3-adic valuation

    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtension_Ring
        sage: isinstance(T, HenselizationExtension_Ring)
        True
        sage: TestSuite(T).run() # long time

    """
    def fraction_field(self):
        r"""
        Return the fraction field of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: S = ZZ.henselization(3)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + 3)
            sage: T.fraction_field()
            Extension defined by x^2 + 3 of Henselization of Rational Field with respect to 3-adic valuation

        """
        base_ring = self.base_ring().fraction_field()
        return Extension(base_ring, self._polynomial.change_ring(base_ring), name=self.variable_name(), check=False)


class HenselizationExtensionAbsolute(HenselizationExtension):
    r"""
    Extension of a :class:`Henselization` whose model is an absolute
    extension.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: L.<x> = K[]
        sage: K.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: L.<x> = K[]
            sage: M = K.extension(x^2 + x + 1)
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionAbsolute
            sage: isinstance(M, HenselizationExtensionAbsolute)
            True
            sage: TestSuite(M).run() # long time

        """
        super(HenselizationExtensionAbsolute, self).__init__(base_ring=base_ring, polynomial=polynomial, model=model, model_valuation=model_valuation, category=category)

        if not model.base_ring() is self._absolute_base_ring().base():
            raise ValueError("model must be an absolute extension of the base of absolute base_ring")

    def _module_basis(self, base):
        r"""
        Return a basis of this ring as a free module over ``base``.
 
        EXAMPLES::

           sage: from henselization import *
           sage: K = QQ.henselization(2)
           sage: R.<x> = K[]
           sage: L = K.extension(x^2 + x + 1)
           sage: L._module_basis(base = K)
           (1, x)

        """
        if base is self._absolute_base_ring():
            gen = self._base(self._base.fraction_field().gen())
            return tuple(self(gen)**i for i in range(self._absolute_degree()))
        return super(HenselizationExtensionAbsolute, self)._module_basis(base = base)


class HenselizationExtensionAbsolute_Field(HenselizationExtensionAbsolute, HenselizationExtension_Field):
    r"""
    Extension of a field :class:`Henselization` whose model is an absolute
    extension.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: L.<x> = K[]
        sage: L = K.extension(x^2 + 2); L
        Extension defined by x^2 + 2 of Henselization of Rational Field with respect to 2-adic valuation
    
    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionAbsolute_Field
        sage: isinstance(L, HenselizationExtensionAbsolute_Field)
        True
        sage: TestSuite(L).run() # long time

    """
    pass


class HenselizationExtensionSimple(HenselizationExtensionAbsolute):
    r"""
    Simple extension of a :class:`Henselization`.

    EXAMPLES::

        sage: from henselization import *
        sage: S = ZZ.henselization(2)
        sage: R.<x> = S[]
        sage: S.extension(x^2 + x + 1)
        Extension defined by x^2 + x + 1 of Henselization of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: S = ZZ.henselization(2)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionSimple
            sage: isinstance(T, HenselizationExtensionSimple)
            True
            sage: TestSuite(T).run() # long time

        """
        if isinstance(base_ring, HenselizationExtension):
            raise TypeError("base_ring must not be an extension")
        super(HenselizationExtensionSimple, self).__init__(base_ring=base_ring, polynomial=polynomial, model=model, model_valuation=model_valuation, category=category)

    def gen(self, i=0):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<a> = K[]
            sage: L = K.extension(a^2 + a + 1)
            sage: L.gen(0)
            a
            sage: L.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: ring has only one generator

        """
        if i != 0:
            raise ValueError("ring has only one generator")
        return self(self._base_fraction_field.gen())


class HenselizationExtensionSimple_Field(HenselizationExtensionSimple, HenselizationExtensionAbsolute_Field):
    r"""
    Simple extension of a :class:`Henselization_Field` that is not an extension.

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1); L
        Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionSimple_Field
        sage: isinstance(L, HenselizationExtensionSimple_Field)
        True
        sage: TestSuite(L).run() # long time

    """


class HenselizationExtensionSimple_Ring(HenselizationExtensionSimple, HenselizationExtension_Ring):
    r"""
    Simple extension of a :class:`Henselization_Ring` that is not an extension.

    EXAMPLES::

        sage: from henselization import *
        sage: S = ZZ.henselization(2)
        sage: R.<x> = S[]
        sage: T = S.extension(x^2 + x + 1); T
        Extension defined by x^2 + x + 1 of Henselization of Integer Ring with respect to 2-adic valuation

    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionSimple_Ring
        sage: isinstance(T, HenselizationExtensionSimple_Ring)
        True
        sage: TestSuite(T).run() # long time

    """


class HenselizationExtensionIteratedQuotient(HenselizationExtension):
    r"""
    Extension of a :class:`HenselizationExtensionAbsolute` that is realized
    as a quotient over its model.

    EXAMPLES:

    Instances of this class are not exposed through public methods but only
    exist internally when building iterated extensions::

        sage: K = QQ.henselization(2)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: from sage.rings.padics.henselization.henselization import Extension
        sage: Extension._create_extension(L, y^2 - 2)
        Extension defined by y^2 - 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: from sage.rings.padics.henselization.henselization import Extension
            sage: M = Extension._create_extension(L, y^2 - 2)
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionIteratedQuotient
            sage: isinstance(M, HenselizationExtensionIteratedQuotient)
            True
            sage: TestSuite(M).run() # long time

        """
        if not isinstance(base_ring, HenselizationExtension):
            raise TypeError("base_ring must be an extension")
        from sage.rings.polynomial.polynomial_quotient_ring import is_PolynomialQuotientRing
        if not is_PolynomialQuotientRing(model):
            raise TypeError("model must be a quotient")

        super(HenselizationExtensionIteratedQuotient, self).__init__(base_ring=base_ring, polynomial=polynomial, model=model, model_valuation=model_valuation, category=category)

        if not model.base_ring() is base_ring.base():
            raise ValueError("model must be a quotient over the base of base_ring")

    def gen(self, i=0):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: from sage.rings.padics.henselization.henselization import Extension
            sage: M = Extension._create_extension(L, y^2 - 2)
            sage: M.gen(0)
            ybar
            sage: M.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: ring has only one generator

        """
        if i != 0:
            raise ValueError("ring has only one generator")
        return self(self._base_fraction_field.gen())

    @cached_method
    def _module_basis(self, base):
        r"""
        Return a basis of this ring as a free module over ``base``.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: from sage.rings.padics.henselization.henselization import Extension
            sage: M = Extension._create_extension(L, y^4 - 2)
            sage: M._module_basis(L)
            (1, ybar, ybar^2, ybar^3)
            sage: M._module_basis(K)
            (1, x, ybar, x*ybar, ybar^2, x*ybar^2, ybar^3, x*ybar^3)

        """
        if base is self.base_ring():
            return tuple(self(self.base().gen()**i) for i in range(self.base().degree())) 
        if base is self._absolute_base_ring():
            return tuple(b*self(p) for b in self._module_basis(self.base_ring()) for p in self.base_ring()._module_basis(base))
        return super(HenselizationExtensionIteratedQuotient, self)._module_basis(base = base)

    def _coerce_map_from_(self, other):
        if isinstance(other, HenselizationExtensionIteratedQuotient):
            if self.base_ring().has_coerce_map_from(other.base_ring()):
                if other._polynomial.change_ring(self.base_ring()) == self._polynomial:
                    homspace = other.Hom(self)
                    from maps import QuotientConversion_generic
                    return homspace.__make_element_class__(QuotientConversion_generic)(homspace)
                    
        return super(HenselizationExtensionIteratedQuotient, self)._coerce_map_from_(other)


class HenselizationExtensionIteratedQuotient_Field(HenselizationExtensionIteratedQuotient, HenselizationExtension_Field):
    r"""
    Extension of a :class:`HenselizationExtensionAbsolute` field that is realized
    as a quotient over its model.

    EXAMPLES:

    Instances of this class are not exposed through public methods but only
    exist internally when building iterated extensions::

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: from sage.rings.padics.henselization.henselization import Extension
        sage: M = Extension._create_extension(L, y^4 - 2); M
        Extension defined by y^4 - 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionIteratedQuotient_Field
        sage: isinstance(M, HenselizationExtensionIteratedQuotient_Field)
        True
        sage: TestSuite(M).run() # long time

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        super(HenselizationExtensionIteratedQuotient_Field, self).__init__(base_ring, polynomial, model,model_valuation,category)


class HenselizationExtensionIteratedQuotient_Ring(HenselizationExtensionIteratedQuotient, HenselizationExtension_Ring):
    r"""
    Extension of a :class:`HenselizationExtensionAbsolute` domain that is not a
    field that is realized as a quotient over its model.

    EXAMPLES:

    Instances of this class are not exposed through public methods but only
    exist internally when building iterated extensions::

        sage: from henselization import *
        sage: S = ZZ.henselization(2)
        sage: R.<x> = S[]
        sage: T = S.extension(x^2 + x + 1)
        sage: R.<y> = T[]
        sage: U.<y> = T.extension(y^4 - 2); U # known bug, see https://github.com/MCLF/henselization/issues/15
        Extension defined by y^4 - 2 of Extension defined by x^2 + x + 1 of Henselization of Integer Ring with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        TESTS::

            sage: from henselization import *
            sage: S = ZZ.henselization(2)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: R.<y> = T[]
            sage: U.<y> = T.extension(y^4 - 2) # known bug, see https://github.com/MCLF/henselization/issues/15
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionIteratedQuotient_Ring
            sage: isinstance(U, HenselizationExtensionIteratedQuotient_Ring) # known bug, see https://github.com/MCLF/henselization/issues/15
            True
            sage: TestSuite(U).run() # long time, known bug, see https://github.com/MCLF/henselization/issues/15

        """
        super(HenselizationExtensionIteratedQuotient_Ring, self).__init__(base_ring, polynomial, model,model_valuation,category)

    def fraction_field(self):
        r"""
        Return the fraction field of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: S = ZZ.henselization(2)
            sage: R.<x> = S[]
            sage: T = S.extension(x^2 + x + 1)
            sage: R.<y> = T[]
            sage: U.<y> = T.extension(y^4 - 2) # known bug, see https://github.com/MCLF/henselization/issues/15
            sage: U.fraction_field() # known bug, see https://github.com/MCLF/henselization/issues/15
            Extension defined by y^4 - 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

        """
        K = self.base_ring().fraction_field()
        return Extension._create_extension(K, self._polynomial.change_ring(K))


class HenselizationExtensionIteratedAbsolute(HenselizationExtensionAbsolute):
    r"""
    Extension of a :class:`HenselizationExtension` that is realized
    as an absolute extension.

    EXAMPLES:

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: L.extension(y^2 - 2)
        Extension defined by y^2 - 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    """
    def __init__(self, base_ring, polynomial, model, model_valuation, category=None):
        r"""
        TESTS::
    
            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<x> = K[]
            sage: L = K.extension(x^2 + x + 1)
            sage: R.<y> = L[]
            sage: M = L.extension(y^2 - 2)
            sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionIteratedAbsolute
            sage: isinstance(M, HenselizationExtensionIteratedAbsolute)
            True
            sage: TestSuite(M).run() # long time
    
        """
        if not isinstance(base_ring, HenselizationExtension):
            raise TypeError("base_ring must be an extension")
        super(HenselizationExtensionIteratedAbsolute, self).__init__(base_ring=base_ring, polynomial=polynomial, model=model, model_valuation=model_valuation, category=category)

    @cached_method
    def gen(self, i=0):
        r"""
        Return the ``i``-th generator of this ring.

        EXAMPLES::

            sage: from henselization import *
            sage: K = QQ.henselization(2)
            sage: R.<a> = K[]
            sage: L = K.extension(a^2 + a + 1)
            sage: L.gen(0)
            a
            sage: L.gen(1)
            Traceback (most recent call last):
            ...
            ValueError: ring has only one generator

        """
        if i != 0:
            raise ValueError("ring has only one generator")

        from generator_element import GeneratorElement
        return self.__make_element_class__(GeneratorElement)(self)


class HenselizationExtensionIteratedAbsolute_Field(HenselizationExtensionIteratedAbsolute, HenselizationExtensionAbsolute_Field):
    r"""
    Extension of a :class:`HenselizationExtension` field that is realized as an
    absolute extension.

    EXAMPLES:

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: R.<x> = K[]
        sage: L = K.extension(x^2 + x + 1)
        sage: R.<y> = L[]
        sage: M = L.extension(y^4 - 2); M
        Extension defined by y^4 - 2 of Extension defined by x^2 + x + 1 of Henselization of Rational Field with respect to 2-adic valuation

    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionIteratedAbsolute_Field
        sage: isinstance(M, HenselizationExtensionIteratedAbsolute_Field)
        True
        sage: TestSuite(M).run() # long time

    """


class HenselizationExtensionIteratedAbsolute_Ring(HenselizationExtensionIteratedAbsolute, HenselizationExtension_Ring):
    r"""
    Extension of a :class:`HenselizationExtension` that is not a field and
    realized as an absolute extension.

    EXAMPLES::

        sage: from henselization import *
        sage: S = ZZ.henselization(2)
        sage: R.<x> = S[]
        sage: T = S.extension(x^2 + x + 1)
        sage: R.<y> = T[]
        sage: U = T.extension(y^4 - 2); U # known bug, see https://github.com/MCLF/henselization/issues/15
        Extension defined by y^4 - 2 of Extension defined by x^2 + x + 1 of Henselization of Integer Ring with respect to 2-adic valuation

    TESTS::

        sage: from sage.rings.padics.henselization.henselization import HenselizationExtensionIteratedAbsolute_Ring
        sage: isinstance(U, HenselizationExtensionIteratedAbsolute_Ring) # known bug, see https://github.com/MCLF/henselization/issues/15
        True
        sage: TestSuite(U).run() # long time, known bug, see https://github.com/MCLF/henselization/issues/15

    """
