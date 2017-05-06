*This project is far from stable. If you want to use it, please contact the author. I would be very happy to assist you at trying this out.*

*This requires some changes to the Sage library to work. The necessary changes are tracked here: https://trac.sagemath.org/ticket/22956*

All p-adic rings currently available in Sage are backed by elements that consist of an approximation and a precision. It can sometimes be tedious to work with such elements since error propagation needs to be analyzed and many algorithms in Sage don't play nice with such inexact elements.

This package implements an exact alternative where p-adic rings are backed by absolute number fields. The idea is to use exact elements in number fields where possible and describe algebraic elements as limits of Mac Lane valuations (see http://github.com/saraedum/mac_lane). Since computations in number fields (and in particular in relative extensions) are slow in Sage, extensions are always rewritten as isomorphic rings defined by an absolute number field with defining (Eisenstein) polynomials with small coefficients.

One can of course not expect arithmetic to be as fast as in the inexact p-adic rings but the approach seems to have its merits. With a few tweaks in Sage, this implementation can compute the degree 384 splitting field (unramified of degree 2) of a degree 12 polynomial over Q2 in less than three minutes. 
