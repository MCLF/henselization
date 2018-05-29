from sage.all import QQ
from henselization.henselization import Henselization

class Rational:
    r"""
    Creation of relatively trivial Henselizations.

    TESTS::

        sage: from henselization.benchmarks.henselizations import Rational
        sage: Rational().time_create()

    """
    def time_create(self):
        Henselization(QQ, QQ.valuation(2))
