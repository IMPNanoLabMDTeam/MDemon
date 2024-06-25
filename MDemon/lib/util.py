import os

import numpy as np


def unique_int_1d_unsorted(array):
    values, indices = np.unique(array, return_index=True)
    return array[np.sort(indices)]


def asiterable(obj):
    """Returns `obj` so that it can be iterated over.

    A string is *not* detected as and iterable and is wrapped into a :class:`list`
    with a single element.

    See Also
    --------
    iterable

    """
    if not iterable(obj):
        obj = [obj]
    return obj


def astuple(obj):
    if not isinstance(obj, tuple):
        obj = asiterable(obj)
        obj = tuple(obj)
    return obj


def iterable(obj):
    """Returns ``True`` if `obj` can be iterated over and is *not* a  string
    nor a :class:`NamedStream`"""
    if isinstance(obj, (str,)):
        return False  # avoid iterating over characters of a string

    if hasattr(obj, "next"):
        return True  # any iterator will do
    try:
        len(obj)  # anything else that might work
    except (TypeError, AttributeError):
        return False
    return True


def guess_format(filename):
    """Return the format of `filename`

    The current heuristic simply looks at the filename extension and can work
    around compressed format extensions.

    Parameters
    ----------
    filename : str or stream
        path to the file or a stream, in which case ``filename.name`` is looked
        at for a hint to the format

    Returns
    -------
    format : str
        format specifier (upper case)

    Raises
    ------
    ValueError
        if the heuristics are insufficient to guess a supported format


    .. versionadded:: 0.11.0
       Moved into lib.util

    """
    format = (
        format_from_filename_extension(filename) if not iterable(filename) else "CHAIN"
    )

    return format.upper()


def format_from_filename_extension(filename):
    """Guess file format from the file extension.

    Parameters
    ----------
    filename : str

    Returns
    -------
    format : str

    Raises
    ------
    TypeError
        if the file format cannot be determined
    """
    try:
        root, ext = get_ext(filename)
    except Exception:
        errmsg = (
            f"Cannot determine file format for file '{filename}'.\n"
            f"           You can set the format explicitly with "
            f"'Universe(..., format=FORMAT)'."
        )
        raise TypeError(errmsg) from None
    format = ext.upper()

    return format


def get_ext(filename):
    """Return the lower-cased extension of `filename` without a leading dot.

    Parameters
    ----------
    filename : str

    Returns
    -------
    root : str
    ext : str
    """
    root, ext = os.path.splitext(filename)

    if ext.startswith(os.extsep):
        ext = ext[1:]

    return root, ext.lower()


def wishnotiterable(obj):
    if iterable(obj) and len(obj) == 1 and not isinstance(obj, dict):
        obj = obj[0]

    return obj
