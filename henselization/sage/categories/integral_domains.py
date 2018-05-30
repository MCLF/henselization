class IntegralDomains:
    class ParentMethods:
        def henselization(self, prime):
            from sage.all import Henselization
            return Henselization(self, self.valuation(prime))
