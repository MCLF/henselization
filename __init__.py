# -*- coding: utf-8 -*-
from mac_lane import *

from .completion import Completion, CompleteRing_base, CompleteDomain, CompleteField
from .base_element import BaseElement
from .mac_lane_element import MacLaneElement
import valuation
from .valuation import CompletionValuation


from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("Completion", Completion)
register_factory_unpickle("valuation.Valuation", valuation.Valuation)
