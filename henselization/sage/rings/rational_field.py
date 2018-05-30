class RationalField:
    def henselization(self, prime):
        from sage.rings.padics.henselization.henselization import Henselization
        return Henselization(self, self.valuation(prime))
