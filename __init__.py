# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

# Fix doctests so they work in standalone mode (when invoked with sage -t, they run within the completion/ directory)
import sys, os
if hasattr(sys.modules['__main__'], 'DC') and 'standalone' in sys.modules['__main__'].DC.options.optional:
    sys.path.append(os.path.dirname(os.getcwd()))

from mac_lane import *

from completion import Completion, CompleteRing_base, CompleteDomain, CompleteField, CompleteExtension_base, CompleteExtensionField, CompleteExtensionDomain, Extension
from base_element import BaseElementRing, BaseElementField, BaseElement_base
from mac_lane_element import MacLaneElement
import valuation
from valuation import CompletionValuation
from completions import CompleteDiscreteValuationRings, CompleteDiscreteValuationFields

from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("Completion", Completion)
register_factory_unpickle("valuation.Valuation", valuation.Valuation)
register_factory_unpickle("Extension", Extension)
