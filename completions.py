# -*- coding: utf-8 -*-
r"""
Categories for complete discrete valuation rings.

Eventually, this should be merged with the corresponding category that is
already in Sage. Then, there should also be an axiom Complete on the level of
DiscreteValuationRings.

AUTHORS:

- Julian Rüth (2016-11-23): initial version

"""
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.misc.abstract_method import abstract_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_with_axiom import CategoryWithAxiom, all_axioms
from sage.categories.discrete_valuation import DiscreteValuationRings, DiscreteValuationFields

#TODO: This should be a patch to DicsreteValuationRings.ElementMethods
def is_squarefree(self):
    r"""
    Return whether this element is squarefree.

    EXAMPLES:

    An element is squarefree if its composition into prime factors of its
    parent has no repeated factors. For field elements, this is therefore
    trivial since there are no prime elements::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 2)
        sage: K = Completion(QQ, v)
        sage: x = K(4)
        sage: x.is_squarefree()
        True

    Over rings, the only prime is the uniformizing element::

        sage: v = pAdicValuation(ZZ, 2)
        sage: K = Completion(ZZ, v)
        sage: K(4).is_squarefree()
        False
        sage: K(9).is_squarefree()
        True

    """
    from sage.categories.fields import Fields
    if self.parent() in Fields:
        return True
    return self.valuation() <= self.parent().valuation().value_group().gen()
DiscreteValuationRings().element_class.is_squarefree = is_squarefree
from sage.categories.fields import Fields
Fields().element_class.is_squarefree = lambda self: True

class CompleteDiscreteValuationRings(Category_singleton):
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationRings().super_categories()
            [Category of discrete valuation rings]
        """
        return [DiscreteValuationRings()]


class CompleteDiscreteValuationFields(Category_singleton):
    def super_categories(self):
        """
        EXAMPLES::

            sage: CompleteDiscreteValuationFields().super_categories()
            [Category of discrete valuation fields]
        """
        return [DiscreteValuationFields()]
