class RationalField:
    def henselization(self, prime):
        r"""
        Return the Henselization of the rationals at the prime number or ideal
        ``prime``.

        EXAMPLES::

            sage: from henselization import *
            sage: H = QQ.henselization(5)

        ::

            sage: R.<x> = QQ[]
            sage: f = x^2 + 1
            sage: f.is_irreducible()
            True
            sage: f.change_ring(H).is_irreducible()
            False

        """
        from sage.rings.padics.henselization.henselization import Henselization
        return Henselization(self, self.valuation(prime))
