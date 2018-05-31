<!-- [![Documentation Status](https://readthedocs.org/projects/henselization/badge/)](http://henselization.readthedocs.io/?badge=latest) -->
[![CircleCI](https://circleci.com/gh/MCLF/henselization/tree/master.svg?style=svg)](https://circleci.com/gh/MCLF/henselization/tree/master)
[![Coverage Status](https://coveralls.io/repos/github/MCLF/henselization/badge.svg?branch=master)](https://coveralls.io/github/MCLF/henselization?branch=master)
[![asv](http://img.shields.io/badge/benchmarked%20by-asv-green.svg?style=flat)](https://mclf.github.io/henselization-asv)

### Henselizations in Sage

This project is in an early alpha stage. It could already be useful but there
are likely quite some issue (that we'd like to hear about 🙂)

All p-adic rings currently available in Sage are backed by elements that
consist of an approximation and a precision. It can sometimes be tedious to
work with such elements since error propagation needs to be analyzed and many
algorithms in Sage don't play nice with such inexact elements.

This package implements an exact alternative where p-adic rings are backed by
absolute number fields. The idea is to use exact elements in number fields
where possible and describe algebraic elements as
[limits of Mac Lane valuations](https://doc.sagemath.org/html/en/reference/valuations/sage/rings/valuation/limit_valuation.html).
Since computations in number fields (and in particular in relative extensions)
are slow in Sage, extensions are always rewritten as isomorphic rings defined
by an absolute number field with defining (Eisenstein) polynomials with small
coefficients.

You need at least [Sage 8.2](https://www.sagemath.org) for the following examples to work.

If you can not install Sage on your local machine, you can also click
[![Launch on mybinder.org](https://camo.githubusercontent.com/d57df63fab21897847014ebaec3e7f5f48951ad2/68747470733a2f2f626574612e6d7962696e6465722e6f72672f62616467652e737667)](https://mybinder.org/v2/gh/mclf/henselization/master?filepath=example.ipynb)
to try this in an interactive Jupyter notebook.

The package can be loaded with
```
sage: from henselization import *
```

#### Example: Splitting Fields

One can of course not expect arithmetic to be as fast as in the inexact p-adic
rings but the approach seems to have its merits. With a few
[tweaks](https://github.com/MCLF/henselization/issues/17) in Sage, this
implementation can compute the degree 384 splitting field (unramified of degree
2) of a degree 12 polynomial over ℚ<sub>2</sub> in less than three minutes (having
`PYTHONOPTIMIZE=yes` set):

```
sage: from henselization import *
sage: from henselization.benchmarks.splitting_fields import splitting_field

sage: K = QQ.henselization(2)
sage: R.<T> = K[]
sage: f = T^12 - 4*T^11 + 2*T^10 + 13*T^8 - 16*T^7 - 36*T^6 + 168*T^5 - 209*T^4 + 52*T^3 + 26*T^2 + 8*T - 13
sage: splitting_field(f)
Extension defined by a192^4 + … of Extension defined by a48^4 + … of Extension defined by a12^12 - 4*a12^11 + 2*a12^10 + 13*a12^8 - 16*a12^7 - 36*a12^6 + 168*a12^5 - 209*a12^4 + 52*a12^3 + 26*a12^2 + 8*a12 - 13 of Extension defined by z2^2 + z2 + 1 of Henselization of Rational Field with respect to 2-adic valuation
```

#### Known bugs and issues

See our [issues list](https://github.com/MCLF/henselization/issues), and tell us of any bugs or ommissions that are not covered there.
