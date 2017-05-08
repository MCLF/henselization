*This project is far from stable. If you want to use it, please contact the author. I would be very happy to assist you at trying this out.*

*This requires some changes to the Sage library to work. The necessary changes are tracked here: https://trac.sagemath.org/ticket/22956*

All p-adic rings currently available in Sage are backed by elements that consist of an approximation and a precision. It can sometimes be tedious to work with such elements since error propagation needs to be analyzed and many algorithms in Sage don't play nice with such inexact elements.

This package implements an exact alternative where p-adic rings are backed by absolute number fields. The idea is to use exact elements in number fields where possible and describe algebraic elements as limits of Mac Lane valuations (see http://github.com/saraedum/mac_lane). Since computations in number fields (and in particular in relative extensions) are slow in Sage, extensions are always rewritten as isomorphic rings defined by an absolute number field with defining (Eisenstein) polynomials with small coefficients.

Splitting Fields
================

One can of course not expect arithmetic to be as fast as in the inexact p-adic rings but the approach seems to have its merits. With a few tweaks in Sage, this implementation can compute the degree 384 splitting field (unramified of degree 2) of a degree 12 polynomial over Q2 in a little more than four minutes:

```
import sys, os
sys.path.append(os.getcwd())
sys.path.append(os.path.dirname(os.getcwd()))
from mac_lane import *
from completion import *

def splitting_field(polynomial, ramified_variable_name = 'a', unramified_variable_name = 'u'):
    ret = polynomial.base_ring()

    while True:
        F = polynomial.change_ring(ret)
        absolute_degree = ret.base().degree()
        ramified_degree = ret.valuation().value_group().index(polynomial.base_ring().valuation().value_group())
        unramified_degree = ZZ(absolute_degree/ramified_degree)

        print "Factoring %s over a field of degree %s * %s…"%(polynomial, unramified_degree, ramified_degree)
        F = list(F.factor())
        F = [f for f,e in F]
        F = sorted(F, key=lambda f:-f.degree())
        print "…factors with degrees %s"%([f.degree() for f in F],)

        for f in F:
            if f.degree() == 1:
                continue

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
                ret = ret.extension(f, ramified_variable_name + str(ramified_degree + f.degree()))
                break
            else:
                break
    return ret

v = pAdicValuation(QQ, 2)
C = Completion(QQ, v)
R.<T> = C[]
f = T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13
print(splitting_field(f))
```
