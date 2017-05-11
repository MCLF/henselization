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
from sage.structure.factory import UniqueFactory

import mac_lane.valuation

class CompletionValuationFactory(UniqueFactory):
    r"""
    Return the valuation on the completion ``domain``.

    EXAMPLES:

    Do not call this factory directly, but call the ``valuation`` method of a
    completion::

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: K.valuation() # indirect doctest
        5-adic valuation

    """
    def create_key(self, domain):
        r"""
        Return a key that uniquely identifies this valuation.

        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

        sage: sys.path.append(os.getcwd()); from completion import *
        sage: v = pAdicValuation(QQ, 5)
        sage: K = Completion(QQ, v)
        sage: K.valuation()
        5-adic valuation

    """
    def __init__(self, parent):
        r"""
        TESTS::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 5)
            sage: K = Completion(QQ, v)
            sage: v = K.valuation()
            sage: isinstance(v, CompletionValuation)
            True
            sage: TestSuite(v).run() # long time

        """
        super(CompletionValuation, self).__init__(parent)
        self._base_valuation = self.domain()._base_valuation

    def _repr_(self):
        r"""
        Return a printable representation of this valuation.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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

            sage: sys.path.append(os.getcwd()); from completion import *
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
        return super(CompletionValuation, self).extensions(ring)

    def _relative_size(self, x):
        r"""
        Return an estimate on the coefficient size of ``x``.

        The number returned is an estimate on the factor between the number of
        Bits used by ``x`` and the minimal number of bits used by an element
        congruent to ``x``.

        This is used by :meth:`simplify` to decide whether simplification of
        coefficients is going to lead to a significant shrinking of the
        coefficients of ``x``.

        EXAMPLES:: 

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation()._relative_size(1024)
            11

        """
        return self.domain().coerce(x)._relative_size()

    def simplify(self, x, error=None, force=False):
        r"""
        Return a simplified version of ``x``.

        Produce an element which differs from ``x`` by an element of
        valuation strictly greater than the valuation of ``x`` (or strictly
        greater than ``error`` if set.)

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation().simplify(1025, force=True)
            1

        """
        return self.domain().coerce(x).simplify(error, force=force)

    def upper_bound(self, x):
        r"""
        Return an upper bound of this valuation at ``x``.

        Use this method to get an approximation of the valuation of ``x``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation().upper_bound(1025)
            0

        """
        return self.domain().coerce(x)._upper_bound()

    def lower_bound(self, x):
        r"""
        Return an lower bound of this valuation at ``x``.

        Use this method to get an approximation of the valuation of ``x``
        when speed is more important than accuracy.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation().lower_bound(1025)
            0

        """
        return self.domain().coerce(x)._lower_bound()

    def restriction(self, ring):
        r"""
        Return the restriction of this valuation to ``ring``.

        EXAMPLES::

            sage: sys.path.append(os.getcwd()); from completion import *
            sage: v = pAdicValuation(QQ, 2)
            sage: K = Completion(QQ, v)
            sage: K.valuation().restriction(ZZ)
            2-adic valuation
        
        """
        if ring.is_subring(self.domain()._base):
            return self.domain()._base_valuation.restriction(ring)
        return super(CompletionValuation, self).restriction(ring)
