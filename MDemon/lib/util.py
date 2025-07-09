import os
from functools import wraps

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


def check_coords(*coord_names, **options):
    """
    Decorator for automated coordinate array checking.

    This decorator is derived from :func:`MDAnalysis.lib.util.check_coords`
    with minimal neccessary modifications, and is intended for use especially
    in :mod:`MDemon.lib.distances`.
    It takes an arbitrary number of positional arguments which must correspond
    to names of positional arguments of the decorated function.
    It then checks if the corresponding values are valid coordinate arrays or
    an :class:`~MDemon.core.structure.Structure`.
    If the input is an array and all these arrays are single coordinates
    (i.e., their shape is ``(3,)``), the decorated function can optionally
    return a single coordinate (or angle) instead of an array of coordinates
    (or angles). This can be used to enable computations of single observables
    using functions originally designed to accept only 2-d coordinate arrays.

    If the input is an :class:`~MDemon.core.structure.Structure` it is
    converted into its corresponding position array via a call to
    `Structure.coordinate`.

    The checks performed on each individual coordinate array are:

    * Check that coordinate arrays are of type :class:`numpy.ndarray`.
    * Check that coordinate arrays have a shape of ``(n, 3)`` (or ``(3,)`` if
      single coordinates are allowed; see keyword argument `allow_single`).
    * Automatic dtype conversion to ``numpy.float32``.
    * Optional replacement by a copy; see keyword argument `enforce_copy` .
    * If coordinate arrays aren't C-contiguous, they will be automatically
      replaced by a C-contiguous copy.
    * Optional check for equal length of all coordinate arrays; see optional
      keyword argument `check_lengths_match`.

    Parameters
    ----------
    *coord_names : tuple
        Arbitrary number of strings corresponding to names of positional
        arguments of the decorated function.
    **options : dict, optional
        * **enforce_copy** (:class:`bool`, optional) -- Enforce working on a
          copy of the coordinate arrays. This is useful to ensure that the input
          arrays are left unchanged. Default: ``True``
        * **enforce_dtype** (:class:`bool`, optional) -- Enforce a conversion
          to float32.  Default: ``True``
        * **allow_single** (:class:`bool`, optional) -- Allow the input
          coordinate array to be a single coordinate with shape ``(3,)``.
        * **convert_single** (:class:`bool`, optional) -- If ``True``, single
          coordinate arrays will be converted to have a shape of ``(1, 3)``.
          Only has an effect if `allow_single` is ``True``. Default: ``True``
        * **reduce_result_if_single** (:class:`bool`, optional) -- If ``True``
          and *all* input coordinates are single, a decorated function ``func``
          will return ``func()[0]`` instead of ``func()``. Only has an effect if
          `allow_single` is ``True``. Default: ``True``
        * **check_lengths_match** (:class:`bool`, optional) -- If ``True``, a
          :class:`ValueError` is raised if not all coordinate arrays contain the
          same number of coordinates. Default: ``True``
        * **allow_structure** (:class:`bool`, optional) -- If ``False``, a
          :class:`TypeError` is raised if an :class:`Structure` is supplied
          Default: ``False``

    Raises
    ------
    ValueError
        If the decorator is used without positional arguments (for development
        purposes only).

        If any of the positional arguments supplied to the decorator doesn't
        correspond to a name of any of the decorated function's positional
        arguments.

        If any of the coordinate arrays has a wrong shape.
    TypeError
        If any of the coordinate arrays is not a :class:`numpy.ndarray` or an
        :class:`~MDemon.core.structure.Structure`.

        If the dtype of any of the coordinate arrays is not convertible to
          ``numpy.float32``.
    """
    enforce_copy = options.get("enforce_copy", True)
    enforce_dtype = options.get("enforce_dtype", True)
    allow_single = options.get("allow_single", True)
    convert_single = options.get("convert_single", True)
    reduce_result_if_single = options.get("reduce_result_if_single", True)
    check_lengths_match = options.get("check_lengths_match", len(coord_names) > 1)
    allow_structure = options.get("allow_structure", False)
    if not coord_names:
        raise ValueError(
            "Decorator check_coords() cannot be used without " "positional arguments."
        )

    def check_coords_decorator(func):
        fname = func.__name__
        code = func.__code__
        varnames = code.co_varnames
        nargs = code.co_argcount
        argnames = varnames[:nargs]
        ndefaults = len(func.__defaults__) if func.__defaults__ else 0
        # Create a tuple of positional argument names:
        nposargs = nargs - ndefaults
        posargnames = varnames[:nposargs]
        # The check_coords() decorator is designed to work only for
        # positional arguments:
        for name in coord_names:
            if name not in posargnames:
                raise ValueError(
                    f"In decorator check_coords(): Name '{name}' "
                    "doesn't correspond to any positional "
                    f"argument of the decorated function {func.__name__}()."
                )

        def _check_coords(coords, argname):
            is_single = False
            if isinstance(coords, np.ndarray):
                if allow_single:
                    if (coords.ndim not in (1, 2)) or (coords.shape[-1] != 3):
                        errmsg = (
                            f"{fname}(): {argname}.shape must be (3,) or "
                            f"(n, 3), got {coords.shape}"
                        )
                        raise ValueError(errmsg)
                    if coords.ndim == 1:
                        is_single = True
                        if convert_single:
                            coords = coords[None, :]
                else:
                    if (coords.ndim != 2) or (coords.shape[1] != 3):
                        errmsg = (
                            f"{fname}(): {argname}.shape must be (n, 3) "
                            f"got {coords.shape}"
                        )
                        raise ValueError(errmsg)
                if enforce_dtype:
                    try:
                        coords = coords.astype(np.float32, order="C", copy=enforce_copy)
                    except ValueError:
                        errmsg = (
                            f"{fname}(): {argname}.dtype must be"
                            f"convertible to float32, got"
                            f" {coords.dtype}."
                        )
                        raise TypeError(errmsg) from None
                # coordinates should now be the right shape
                ncoord = coords.shape[0]
            else:
                try:
                    coords = coords.coordinate  # homogenise to a numpy array
                    if coords.ndim == 1:
                        is_single = True
                        coords = coords[None, :]
                    ncoord = coords.shape[0]
                    if not allow_structure:
                        err = TypeError(
                            "Structure or other class with a"
                            "`.positions` method supplied as an"
                            "argument, but allow_structure is"
                            " False"
                        )
                        raise err
                except AttributeError:
                    raise TypeError(
                        f"{fname}(): Parameter '{argname}' must be"
                        f" a numpy.ndarray or an Structure,"
                        f" got {type(coords)}."
                    )

            return coords, is_single, ncoord

        @wraps(func)
        def wrapper(*args, **kwargs):
            # Check for invalid function call:
            if len(args) != nposargs:
                # set marker for testing purposes:
                wrapper._invalid_call = True
                if len(args) > nargs:
                    # too many arguments, invoke call:
                    return func(*args, **kwargs)
                for name in posargnames[: len(args)]:
                    if name in kwargs:
                        # duplicate argument, invoke call:
                        return func(*args, **kwargs)
                for name in posargnames[len(args) :]:
                    if name not in kwargs:
                        # missing argument, invoke call:
                        return func(*args, **kwargs)
                for name in kwargs:
                    if name not in argnames:
                        # unexpected kwarg, invoke call:
                        return func(*args, **kwargs)
                # call is valid, unset test marker:
                wrapper._invalid_call = False
            args = list(args)
            ncoords = []
            all_single = allow_single
            for name in coord_names:
                idx = posargnames.index(name)
                if idx < len(args):
                    args[idx], is_single, ncoord = _check_coords(args[idx], name)
                    all_single &= is_single
                    ncoords.append(ncoord)
                else:
                    kwargs[name], is_single, ncoord = _check_coords(kwargs[name], name)
                    all_single &= is_single
                    ncoords.append(ncoord)
            if check_lengths_match and ncoords:
                if ncoords.count(ncoords[0]) != len(ncoords):
                    raise ValueError(
                        "{}(): {} must contain the same number of "
                        "coordinates, got {}."
                        "".format(fname, ", ".join(coord_names), ncoords)
                    )
            # If all input coordinate arrays were 1-d, so should be the output:

            if all_single and reduce_result_if_single:
                return func(*args, **kwargs)[0]
            return func(*args, **kwargs)

        return wrapper

    return check_coords_decorator


def check_box(box):
    """
    Take a box input and deduce what type of system it represents based on
    the shape of the array and whether all angles are 90 degrees.
    This function is derived from :func:`MDAnalysis.lib.util.check_box`
    with minimal neccessary modifications.

    Parameters
    ----------
    box : array_like
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.

    Returns
    -------
    boxtype : {``'ortho'``, ``'tri_vecs'``}
        String indicating the box type (orthogonal or triclinic).
    checked_box : numpy.ndarray
        Array of dtype ``numpy.float32`` containing box information:
          * If `boxtype` is ``'ortho'``, `cecked_box` will have the shape ``(3,)``
            containing the x-, y-, and z-dimensions of the orthogonal box.
          * If  `boxtype` is ``'tri_vecs'``, `cecked_box` will have the shape
            ``(3, 3)`` containing the triclinic box vectors in a lower triangular
            matrix as returned by
            :meth:`~MDemon.lib.mdamath.triclinic_vectors`.

    Raises
    ------
    ValueError
        If `box` is not of the form ``[lx, ly, lz, alpha, beta, gamma]``
        or contains data that is not convertible to ``numpy.float32``.

    See Also
    --------
    MDemon.lib.mdamath.triclinic_vectors
    """
    if box is None:
        raise ValueError("Box is None")
    from .mdamath import triclinic_vectors  # avoid circular import

    box = np.asarray(box, dtype=np.float32, order="C")
    if box.shape != (6,):
        raise ValueError(
            "Invalid box information. Must be of the form "
            "[lx, ly, lz, alpha, beta, gamma]."
        )
    if np.all(box[3:] == 90.0):
        return "ortho", box[:3]
    return "tri_vecs", triclinic_vectors(box)


def flat(nums):
    """
    Flat the multi-layer nesting list.
    """
    res = []
    nums = asiterable(nums)
    for i in nums:
        if isinstance(i, dict):
            i = list(i.keys())

        if isinstance(i, list):
            res.extend(flat(i))
        else:
            res.append(i)
    return res
