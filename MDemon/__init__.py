import periodictable


class ddict(dict):
    """
    A 'disposable' :class:`dict`.
    """

    def __setitem__(self, key, value):
        if key not in self:
            super().__setitem__(key, value)


class dlist(list):
    """
    A 'disposable' :class:`list`.
    """

    def append(self, value):
        for s in self:
            if s.__name__ == value.__name__:
                break
        else:
            super().append(value)


# Registry of Readers, Writers and Structures known to MDemon.
# Metaclass magic fills these as classes are declared.

_READERS = ddict()
_READER_HINTS = ddict()
_STRUCTURE_NAMES = ddict()
_STRUCTURES = dlist()
_STRUCTURE_ATTRS = dlist()
_PARTICLES = dlist()
_TOPOLOGIES = dlist()
_UNIVERSE_ATTRS = dlist()

_ELEMENTS = [
    (element.mass, element)
    for element in periodictable.elements
    if element.mass is not None
]
_ELEMENTS.sort(key=lambda x: x[0])  # 按质量排序

from . import analysis, constants
from .core import *
from .reader import *
from .utils import *
