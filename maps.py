from sage.categories.morphism import Morphism

class BaseExtensionCoercion(Morphism):
    def __init__(self, parent):
        Morphism.__init__(self, parent)

    def is_injective(self):
        return True

    def _call_(self, x):
        from base_element import BaseElement_base
        if isinstance(x, BaseElement_base):
            return self.codomain()(x._x)
        raise NotImplementedError
