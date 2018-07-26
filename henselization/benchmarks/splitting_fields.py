# -*- coding: utf-8 -*-
r"""
Benchmarks for the constrtuction of splitting fields.
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
        K._splitting_field_univariate_polynomial(f)

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
        K._splitting_field_univariate_polynomial(f)

