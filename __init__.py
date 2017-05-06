# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2016 Julian Rüth <julian.rueth@fsfe.org>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from mac_lane import *

from completion import Completion, Completion_base, Completion_Ring, Completion_Field, CompletionExtension, CompletionExtension_Field, CompletionExtension_Ring, CompletionExtensionAbsolute, CompletionExtensionAbsolute_Field, CompletionExtensionIteratedAbsolute, CompletionExtensionIteratedAbsolute_Field, CompletionExtensionIteratedAbsolute_Ring, CompletionExtensionIteratedQuotient, CompletionExtensionIteratedQuotient_Field, CompletionExtensionIteratedQuotient_Ring, CompletionExtensionSimple, CompletionExtensionSimple_Field, CompletionExtensionSimple_Ring, Extension
from completion_element import CompletionElement_base, CompletionElement_Field, CompletionElement_Ring
from base_element import BaseElement_Ring, BaseElement_Field, BaseElement_base
from mac_lane_element import MacLaneElement_base, MacLaneElement_Field, MacLaneElement_Ring
import valuation
from valuation import CompletionValuation
from completions import CompleteDiscreteValuationRings, CompleteDiscreteValuationFields
from maps import ConvertMap_generic, CompletionToVectorSpace, VectorSpaceToCompletion, VectorSpaceCompletionIsomorphism, RelativeExtensionCoercion_generic
from generator_element import GeneratorElement

from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("Completion", Completion)
register_factory_unpickle("valuation.Valuation", valuation.Valuation)
register_factory_unpickle("Extension", Extension)
