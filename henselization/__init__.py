from __future__ import absolute_import
from . import sage as monkey_sage
import sage.all

from recursive_monkey_patch import monkey_patch
monkey_patch(monkey_sage, sage)

from sage.rings.padics.henselization.henselization import Henselization
sage.all.Henselization = Henselization
