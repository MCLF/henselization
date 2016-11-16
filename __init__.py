# -*- coding: utf-8 -*-
from mac_lane import *

from .completion import Completion, CompleteRing_base, CompleteDomain, CompleteField
from .base_element import BaseElement

from sage.structure.factory import register_factory_unpickle
register_factory_unpickle("Completion", Completion)
