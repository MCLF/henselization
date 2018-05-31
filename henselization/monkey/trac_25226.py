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
    _trac = "https://trac.sagemath.org/ticket/25226"

    def _test(self):
        # not testable, see #25226
        pass
    
    def _patch(self):
        import patchy
        import sage.rings.valuation.inductive_valuation
        patchy.patch(sage.rings.valuation.inductive_valuation.NonFinalInductiveValuation.equivalence_decomposition, r"""
@@ -1209,7 +1209,7 @@ class NonFinalInductiveValuation(FiniteInductiveValuation, DiscreteValuation):
         v = self.extension(domain)
         ret = v.equivalence_decomposition(v.domain()(f))
         return Factorization([(self._eliminate_denominators(g), e)
-                              for (g,e) in ret], unit=self._eliminate_denominators(ret.unit()))
+                              for (g,e) in ret], unit=self._eliminate_denominators(ret.unit()), sort=False)
 
        """)


Monkey().patch()
