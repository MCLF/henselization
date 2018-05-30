# -*- coding: utf-8 -*-
r"""
Categories for Henselian discrete valuation rings.

Eventually, there should also be an axiom Henselian on the level of
DiscreteValuationRings.

AUTHORS:

- Julian Rüth (2016-11-23): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016-2018 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.category_singleton import Category_singleton

# TODO: This should be an actual monkey patch. This is not part of the henselization package.
class DiscreteValuationRings:
    class ElementMethods:
        def is_squarefree(self):
            r"""
            Return whether this element is squarefree.
        
            EXAMPLES:
        
            An element is squarefree if its composition into prime factors of its
            parent has no repeated factors. For field elements, this is therefore
            trivial since there are no prime elements::
        
                sage: from henselization import *
                sage: K = QQ.henselization(2)
                sage: x = K(4)
                sage: x.is_squarefree()
                True
        
            Over rings, the only prime is the uniformizing element::
        
                sage: R = ZZ.henselization(2)
                sage: R(4).is_squarefree()
                False
                sage: R(9).is_squarefree()
                True
        
            """
            from sage.categories.fields import Fields
            if self.parent() in Fields:
                return True
            return self.valuation() <= self.parent().valuation().value_group().gen()

class HenselianDiscreteValuationRings(Category_singleton):
    def super_categories(self):
        """
        EXAMPLES::

            sage: from henselization import *
            sage: ZZ.henselization(2).category().super_categories()
            [Category of discrete valuation rings]
        """
        from sage.categories.discrete_valuation import DiscreteValuationRings
        return [DiscreteValuationRings()]


class HenselianDiscreteValuationFields(Category_singleton):
    def super_categories(self):
        """
        EXAMPLES::

            sage: from henselization import *
            sage: QQ.henselization(2).category().super_categories()
            [Category of discrete valuation fields]
        """
        from sage.categories.discrete_valuation import DiscreteValuationFields
        return [DiscreteValuationFields()]
