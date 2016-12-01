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
        Morphism.__init__(self, parent)

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
