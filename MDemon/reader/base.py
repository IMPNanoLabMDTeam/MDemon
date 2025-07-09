import numpy as np

from .. import _READER_HINTS, _READERS
from ..lib import util


class _ReaderMeta(type):
    def __init__(cls, name, bases, classdict):
        type.__init__(type, name, bases, classdict)
        try:
            fmt = util.asiterable(classdict["format"])
        except KeyError:
            pass
        else:
            for fmt_name in fmt:
                fmt_name = fmt_name.upper()
                _READERS[fmt_name] = cls

                if "_format_hint" in classdict:
                    _READER_HINTS[fmt_name] = classdict["_format_hint"].__func__


class ReaderBase(metaclass=_ReaderMeta):
    def __init__(self, filename):
        self.filename = filename

    def parse(self, **kwargs):  # pragma: no cover
        raise NotImplementedError("Override this in each subclass")

    def close(self):
        pass  # pylint: disable=unnecessary-pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # see http://docs.python.org/2/library/stdtypes.html#typecontextmanager
        self.close()
        return False  # do not suppress exceptions


class DynamicReaderBase(ReaderBase):
    pass


def squash_by(child_parent_ids, *attributes):
    """Squash a child-parent relationship

    Arguments
    ---------
    child_parent_ids - array of ids (unique values that identify the parent)
    *attributes - other arrays that need to follow the sorting of ids

    Returns
    -------
    child_parents_idx - an array of len(child) which points to the index of
                        parent
    parent_ids - len(parent) of the ids
    *parent_attrs - len(parent) of the other attributes
    """
    unique_resids, sort_mask, atom_idx = np.unique(
        child_parent_ids, return_index=True, return_inverse=True
    )

    return atom_idx, unique_resids, [attr[sort_mask] for attr in attributes]
