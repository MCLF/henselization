# -*- coding: utf-8 -*-
r"""
Valuations on completions of rings

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

from sage.structure.factory import UniqueFactory

import mac_lane.valuation

class CompletionValuationFactory(UniqueFactory):
    r"""
    Return the valuation on the completion ``domain``.

    EXAMPLES:

    Do not call this factory directly, but call the ``valuation`` method of a
    completion::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: K.valuation() # indirect doctest
        5-adic valuation

    """
    def create_key(self, domain):
        r"""
        Return a key that uniquely identifies this valuation.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: K.valuation() is K.valuation() # indirect doctest
            True

        """
        return domain,

    def create_object(self, version, key):
        r"""
        Return the valuation described by ``key``.

        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: K.valuation() # indirect doctest
            5-adic valuation

        """
        domain, = key
        from mac_lane import DiscretePseudoValuationSpace
        parent = DiscretePseudoValuationSpace(domain)
        return parent.__make_element_class__(CompletionValuation)(parent)

Valuation = CompletionValuationFactory("completion.valuation.Valuation")


class CompletionValuation(mac_lane.valuation.DiscreteValuation):
    r"""
    The valuation on a completion.

    This class turns the valuation that was used to create a completion into a
    valuation on the actual completion.

    EXAMPLES::

        sage: from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: K.valuation()
        5-adic valuation

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: v = K.valuation()
            sage: isinstance(v, CompletionValuation)
            True
            sage: TestSuite(v).run() # long time

        """
        mac_lane.valuation.DiscreteValuation.__init__(self, parent)
        self._base_valuation = self.domain()._base_valuation

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: K.valuation() # indirect doctest
            5-adic valuation

        """
        return repr(self._base_valuation)

    def _call_(self, x):
        r"""
        Evaluation this valuation at ``x``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: w = K.valuation()
            sage: w(K(25))
            2

        """
        return x.valuation()

    def uniformizer(self):
        r"""
        Return a uniformizing element for this valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: w = K.valuation()
            sage: w.uniformizer()
            5

        """
        return self.domain()(self._base_valuation.uniformizer())

    def residue_ring(self):
        r"""
        Return the residue ring of this valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: w = K.valuation()
            sage: w.residue_ring()
            Finite Field of size 5

        """
        return self._base_valuation.residue_ring()

    def lift(self, F):
        r"""
        Lift ``F`` from the residue ring to an element in the domain of this
        valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: w = K.valuation()
            sage: w.lift(1)
            1

        """
        F = self.residue_ring().coerce(F)
        return self.domain()(self._base_valuation.lift(F))

    def reduce(self, f):
        r"""
        Reduce ``f`` modulo the elements of positive valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: w = K.valuation()
            sage: w.reduce(6)
            1

        """
        f = self.domain().coerce(f)
        return f.reduction()

    def value_semigroup(self):
        r"""
        Return the value semigroup of this valuation.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(ZZ, 5)
            sage: R = Completion(ZZ, v)
            sage: w = R.valuation()
            sage: w.value_semigroup()
            Additive Abelian Semigroup generated by 1

        """
        return self.domain()._base_valuation.value_semigroup()

    def extensions(self, ring):
        r"""
        Return the extensions of this valuation to ``ring``.

        EXAMPLES::

            sage: from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: w = K.valuation()
            sage: R.<a> = K[]
            sage: L.<a> = K.extension(a^2 + a + 1)
            sage: w.extensions(L)
            [2-adic valuation]

        """
        if self.domain().is_subring(ring):
            from completion import Completion_base
            if isinstance(ring, Completion_base):
                return [ring.valuation()]
        return super(Valuation, self).extensions(ring)

