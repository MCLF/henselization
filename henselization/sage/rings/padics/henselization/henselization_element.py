# -*- coding: utf-8 -*-
r"""
Base class for elements of Henselizations
=========================================

Concrete elements of a Henselization can be of several kinds depending on how
they were constructed.

Base Elements
-------------

Since any Henselization is backed by an absolute number field internally, the
elemnts of that absolute number field can be understood as elements of the
Henselization::

    sage: from henselization import *

    sage: H = QQ.henselization(2)
    sage: p = H(2); p
    2

These elements provide the underlying arithmetic of number field elements but
also keep track of other properties such as their valuation::

    sage: p.valuation()
    1

See :module:`sage.rings.padics.henselization.base_element` for details about
this type of element.


Mac Lane Elements
-----------------

A polynomial that is irreducible over a certain number field might factor over
its Henselization. The coefficients of that factorization are not elements of
the underlying absolute number field. Internally, they are described by a
limit valuation, see :module:`sage.rings.valuations.limit_valuation`.

    sage: H = QQ.henselization(5)

    sage: R.<x> = QQ[]
    sage: f = x^2 + 1
    sage: f.is_irreducible()
    True
    sage: f = f.change_ring(H)
    sage: f.is_irreducible()
    False

    sage: f.factor()
    (x + 2 + O(5)) * (x + 3 + O(5))

    sage: a = _[0][0][0]; a
    2 + O(5)

Such elements suppport a number of fundamental operations, however, direct
arithmetic with this kind of element is quite limited currently since we do not
want to perform computations symbolically which would be very slow::

    sage: a.valuation()
    0
    sage: a == 2
    False
    sage: a.reduction()
    2
    sage: a + 1
    Traceback (most recent call last):
    ...
    NotImplementedError: addition not implemented for Henselization of Rational Field with respect to 5-adic valuation

When arithmetic is needed, such elements can be approximated by base elements
to arbitrary precision, i.e., essentially with a certain number of p-adic
digits::

    sage: a.approximation(10)
    -3116/237
    sage: _ + 1
    -2879/237

See :module:`sage.rings.padics.henselization.mac_lane_element` for details
about this type of element.


Generator Elements
------------------

Finally, there are elements that represent the generators of relative extensions::

    sage: H = QQ.henselization(3)
    sage: R.<x> = H[]
    sage: L.<a> = H.extension(x^2 - 3)
    sage: R.<y> = L[]
    sage: M.<b> = L.extension(y^2 - a)

Here, ``a`` is a base element given by the underlying absolute number field of
``L``. The element ``b`` is a generator element since the code could not figure
out a consistent choice of base elements that match the roots of the given
polynomials. Consequently, the arithmetic that can be performed with ``b`` is
currently quite limited:

    sage: b + 1
    Traceback (most recent call last):
    ...
    NotImplementedError: ...

However, in this case, this can be patched up quite easily manually by making
an explicit choice of roots of the defining polynomials::

    sage: R.<x> = M._base[]
    sage: f = x^2 - 3
    sage: f.factor()
    (x - e4^2) * (x + e4^2)
    sage: a = M(-_[0][0][0])

    sage: f = x^2 - a
    sage: f.factor()
    (x - e4) * (x + e4)
    sage: b = M(-_[0][0][0])

Now ``a`` and ``b`` are proper base elements in the field ``M``:

    sage: a + b
    e4^2 + e4
    sage: _.valuation()
    1/4

See :module:`sage.rings.padics.henselization.generator_element` for details
about generator elements.


AUTHORS:

- Julian Rüth (2017-05-04): initial version

"""
#*****************************************************************************
#       Copyright (C) 2017 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.structure.element import IntegralDomainElement

class HenselizationElement_base(IntegralDomainElement):
    r"""
    Abstract base class for elements of :class:`Henselization_base`

    EXAMPLES::


        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: x = K(0); x
        0

    TESTS::

        sage: from sage.rings.padics.henselization.henselization_element import HenselizationElement_base
        sage: isinstance(x, HenselizationElement_base)
        True
        sage: TestSuite(x).run()

    """

class HenselizationElement_Field(HenselizationElement_base):
    r"""
    Abstract base class for elements of :class:`Henselization_Field`

    EXAMPLES::

        sage: from henselization import *
        sage: K = QQ.henselization(2)
        sage: x = K(0); x
        0

    TESTS::

        sage: from sage.rings.padics.henselization.henselization_element import HenselizationElement_Field
        sage: isinstance(x, HenselizationElement_Field)
        True
        sage: TestSuite(x).run()

    """

class HenselizationElement_Ring(HenselizationElement_base):
    r"""
    Abstract base class for elements of :class:`Henselization_Ring`

    EXAMPLES::

        sage: from henselization import *
        sage: S = ZZ.henselization(2)
        sage: x = S(0); x
        0

    TESTS::

        sage: from sage.rings.padics.henselization.henselization_element import HenselizationElement_Ring
        sage: isinstance(x, HenselizationElement_Ring)
        True
        sage: TestSuite(x).run()

    """
    def euclidean_degree(self):
        r"""
        Return an Euclidean degree of this element.

        EXAMPLES::

            sage: from henselization import *
            sage: S = ZZ.henselization(3)
            sage: R.<x> = S[]
            sage: T.<y> = S.extension(x^2 + 3)
            sage: y.euclidean_degree()
            1

        """
        if not self:
            raise ValueError("Euclidean degree of the zero element not defined")
        return self.valuation() / self.parent().valuation().value_group().gen()
