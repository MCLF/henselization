# -*- coding: utf-8 -*-

#*****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .util import AbstractMonkey

class Monkey(AbstractMonkey):
    _trac = "https://trac.sagemath.org/ticket/24934"

    def _test(self):
        from sage.all import QuadraticField, loads, dumps
        K = QuadraticField(-1, 'i')
        if loads(dumps(K.maximal_order())) is not K.maximal_order():
            raise Exception("#24934 has not been fixed")
    
    def _patch(self):
        import sage.rings.number_field.order
        sage.rings.number_field.order.AbsoluteOrderFactory = AbsoluteOrderFactory
        sage.rings.number_field.order.absolute_order_from_module_generators = AbsoluteOrderFactory("sage.rings.number_field.order.absolute_order_from_module_generators")
        sage.rings.number_field.order.RelativeOrderFactory = RelativeOrderFactory
        sage.rings.number_field.order.relative_order_from_ring_generators = RelativeOrderFactory("sage.rings.number_field.order.relative_order_from_ring_generators")
        del sage.rings.number_field.order.AbsoluteOrder.__reduce__
        del sage.rings.number_field.order.RelativeOrder.__reduce__
        

from sage.structure.factory import UniqueFactory
class AbsoluteOrderFactory(UniqueFactory):
    def create_key_and_extra_args(self, gens, check_integral=True, check_rank=True, check_is_ring=True, is_maximal=None, allow_subfield=False):
        if allow_subfield:
            raise NotImplementedError("the allow_subfield parameter is not supported yet")
        if len(gens) == 0:
            raise ValueError("gens must span an order over ZZ")

        from sage.all import Sequence
        gens = Sequence(gens)

        K = gens.universe()
        from sage.rings.number_field.order import is_NumberFieldOrder
        if is_NumberFieldOrder(K):
            K = K.number_field()
        gens = frozenset([K(g) for g in gens])
        return (K, gens), {"check_integral": check_integral,
                           "check_rank": check_rank,
                           "check_is_ring": check_is_ring,
                           "is_maximal": is_maximal}

    def create_object(self, version, key, check_integral, check_rank, check_is_ring, is_maximal):
        K, gens = key

        from sage.rings.number_field.order import each_is_integral
        if check_integral and not each_is_integral(gens):
            raise ValueError("each generator must be integral")

        K = iter(gens).next().parent()
        V, from_V, to_V = K.vector_space()
        mod_gens = [to_V(x) for x in gens]

        from sage.all import ZZ
        ambient = ZZ**V.dimension()
        W = ambient.span(mod_gens)

        if check_rank:
            if W.rank() != K.degree():
                raise ValueError("the rank of the span of gens is wrong")

        if check_is_ring:
            # Is there a faster way?
            from sage.rings.monomials import monomials
            alg = [to_V(x) for x in monomials(gens, [f.absolute_minpoly().degree() for f in gens])]
            if ambient.span(alg) != W:
                raise ValueError("the module span of the gens is not closed under multiplication.")

        from sage.rings.number_field.order import AbsoluteOrder
        return AbsoluteOrder(K, W, check=False, is_maximal=is_maximal)  # we have already checked everything

class RelativeOrderFactory(UniqueFactory):
    def create_key_and_extra_args(self, gens, check_is_integral=True, check_rank=True, is_maximal = None, allow_subfield=False):
        if allow_subfield:
            raise NotImplementedError("the allow_subfield parameter is not supported yet")

        from sage.all import Sequence
        gens = Sequence(gens)

        K = gens.universe()
        from sage.rings.number_field.order import is_NumberFieldOrder
        if is_NumberFieldOrder(K):
            K = K.number_field()
        gens = frozenset([K(g) for g in gens])

        return (K, gens), {"check_is_integral": check_is_integral,
                           "check_rank": check_rank,
                           "is_maximal": is_maximal}

    def create_object(self, version, key, check_is_integral, check_rank, is_maximal):
        K, gens = key

        from sage.rings.number_field.order import each_is_integral
        if check_is_integral and not each_is_integral(gens):
            raise ValueError("each generator must be integral")

        # The top number field that contains the order.
        K = iter(gens).next().parent()
    
        # The absolute version of that field.
        Kabs = K.absolute_field('z')
        from_Kabs, to_Kabs = Kabs.structure()
    
        module_gens = [to_Kabs(a) for a in gens]
        n = [a.absolute_minpoly().degree() for a in gens]
        from sage.rings.monomials import monomials
        absolute_order_module_gens = monomials(module_gens, n)
    
        from sage.rings.number_field.order import absolute_order_from_module_generators
        abs_order =  absolute_order_from_module_generators(absolute_order_module_gens,
                                                           check_integral=False, check_is_ring=False,
                                                           check_rank=check_rank)
    
        from sage.rings.number_field.order import RelativeOrder
        return RelativeOrder(K, abs_order, check=False, is_maximal=is_maximal)


Monkey().patch()
