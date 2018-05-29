# -*- coding: utf-8 -*-
r"""
Base class for elements of Henselizations

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

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: x = K(0); x
        0

    TESTS::

        sage: isinstance(x, HenselizationElement_base)
        True
        sage: TestSuite(x).run()

    """

class HenselizationElement_Field(HenselizationElement_base):
    r"""
    Abstract base class for elements of :class:`Henselization_Field`

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = QQ.valuation(2)
        sage: K = Henselization(QQ, v)
        sage: x = K(0); x
        0

    TESTS::

        sage: isinstance(x, HenselizationElement_Field)
        True
        sage: TestSuite(x).run()

    """

class HenselizationElement_Ring(HenselizationElement_base):
    r"""
    Abstract base class for elements of :class:`Henselization_Ring`

    EXAMPLES::

        sage: sys.path.append(os.getcwd()); from henselization import *
        sage: v = ZZ.valuation(2)
        sage: S = Henselization(ZZ, v)
        sage: x = S(0); x
        0

    TESTS::

        sage: isinstance(x, HenselizationElement_Ring)
        True
        sage: TestSuite(x).run()

    """
    def euclidean_degree(self):
        r"""
        Return an Euclidean degree of this element.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from henselization import *
            sage: v = ZZ.valuation(3)
            sage: S = Henselization(ZZ, v)
            sage: R.<x> = S[]
            sage: T.<y> = S.extension(x^2 + 3)
            sage: y.euclidean_degree()
            1

        """
        if not self:
            raise ValueError("Euclidean degree of the zero element not defined")
        return self.valuation() / self.parent().valuation().value_group().gen()
