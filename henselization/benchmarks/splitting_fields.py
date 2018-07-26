# -*- coding: utf-8 -*-
r"""
Benchmarks for the constrtuction of splitting fields
====================================================
"""
#*****************************************************************************
#       Copyright (C) 2018 Julian Rüth <julian.rueth@fsfe.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import
from sage.all import QQ, ZZ, PolynomialRing, GF

def splitting_field(polynomial, ramified_variable_name = 'a', unramified_variable_name = 'u'):
    r"""
    TODO: Eventually we should monkey patch this into polynomials over Henselizations.
    """
    ret = polynomial.base_ring()

    while True:
        F = polynomial.change_ring(ret)
        absolute_degree = ret.base().degree()
        ramified_degree = ret.valuation().value_group().index(polynomial.base_ring().valuation().value_group())
        unramified_degree = ZZ(absolute_degree/ramified_degree)

        print("Factoring %s over a field of degree %s * %s…"%(polynomial, unramified_degree, ramified_degree))
        F = list(F.factor())
        F = [f for f,e in F]
        F = sorted(F, key=lambda f:-f.degree())
        print("…factors with degrees %s"%([f.degree() for f in F],))

        for f in F:
            if f.degree() == 1:
                continue

            from sage.rings.padics.henselization.mac_lane_element import MacLaneElement_base
            if isinstance(f[0], MacLaneElement_base):
                valuation = f[0]._limit_valuation
            else:
                valuation = ret.valuation().mac_lane_approximants(f)[0]

            ramified_part = valuation.value_group().index(ret.valuation().value_group())
            unramified_part = ZZ(f.degree()/ ramified_part)

            if unramified_part != 1:
                print("Found unramified part of degree %s"%unramified_part)
                ret = polynomial.base_ring().extension(GF(polynomial.base_ring().valuation().residue_field().characteristic() ** (unramified_degree * unramified_part)).polynomial().change_ring(QQ).change_ring(polynomial.base_ring()))
                break
        else:
            for f in F:
                if f.degree() == 1:
                    continue
                print("Found totally ramified part of degree %s"%f.degree())
                ret = ret.extension(f, ramified_variable_name + str(ramified_degree * f.degree()))
                break
            else:
                break
    return ret

class SplittingField:
    timeout = 1800

    def setup(self):
        import henselization
        henselization; # silence pyflakes warning about an unused import

    def time_6_8(self):
        r"""
        TESTS::
 
            sage: import henselization
            sage: from henselization.benchmarks.splitting_fields import SplittingField
            sage: SplittingField().time_6_8()
            Factoring T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 14 over a field of degree 1 * 1…
            …factors with degrees [4, 1, 1]
            Found totally ramified part of degree 4
            Factoring T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 14 over a field of degree 1 * 4…
            …factors with degrees [2, 1, 1, 1, 1]
            Found totally ramified part of degree 2
            Factoring T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 14 over a field of degree 1 * 8…
            …factors with degrees [1, 1, 1, 1, 1, 1]

        """
        K = QQ.henselization(2)
        R = PolynomialRing(K, 'T')
        T = R.gen()
        f = T**6 + 168*T**5 - 209*T**4 + 52*T**3 + 26*T**2 + 8*T - 14
        splitting_field(f)

    def time_10_(self):
        r"""
        TESTS::
 
            sage: import henselization
            sage: from henselization.benchmarks.splitting_fields import SplittingField
            sage: SplittingField().time_6_8()
            Factoring T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 14 over a field of degree 1 * 1…
            …factors with degrees [4, 1, 1]
            Found totally ramified part of degree 4
            Factoring T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 14 over a field of degree 1 * 4…
            …factors with degrees [2, 1, 1, 1, 1]
            Found totally ramified part of degree 2
            Factoring T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 14 over a field of degree 1 * 8…
            …factors with degrees [1, 1, 1, 1, 1, 1]

        """
        K = QQ.henselization(2)
        R = PolynomialRing(K, 'T')
        T = R.gen()
        f = T**10 + T**6 + 168*T**5 - 209*T**4 + 52*T**3 + 26*T**2 + 8*T - 14
        splitting_field(f)

    def time_12_384(self):
        r"""
        TESTS::
 
            sage: import henselization
            sage: from henselization.benchmarks.splitting_fields import SplittingField
            sage: SplittingField().time_12_384() # long time
            Factoring T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13 over a field of degree 1 * 1…
            …factors with degrees [12]
            Found totally ramified part of degree 12
            Factoring T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13 over a field of degree 1 * 12…
            …factors with degrees [8, 2, 1, 1]
            Found unramified part of degree 2
            Factoring T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13 over a field of degree 2 * 1…
            …factors with degrees [12]
            Found totally ramified part of degree 12
            Factoring T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13 over a field of degree 2 * 12…
            …factors with degrees [4, 4, 1, 1, 1, 1]
            Found totally ramified part of degree 4
            Factoring T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13 over a field of degree 2 * 48…
            …factors with degrees [4, 1, 1, 1, 1, 1, 1, 1, 1]
            Found totally ramified part of degree 4
            Factoring T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13 over a field of degree 2 * 192…
            …factors with degrees [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
 
        """
        K = QQ.henselization(2)
        R = PolynomialRing(K, 'T')
        T = R.gen()
        f = T**12 - 4*T**11 + 2*T**10 + 13*T**8 - 16*T**7 - 36*T**6 + 168*T**5 - 209*T**4 + 52*T**3 + 26*T**2 + 8*T - 13
        splitting_field(f)

