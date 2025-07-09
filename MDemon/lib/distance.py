from typing import TYPE_CHECKING, Callable, Optional, Union

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:  # pragma: no cover
    from ..core.structure import Structure

# hack to select backend with backend=<backend> kwarg. Note that
# the cython parallel code (prange) in parallel.distances is
# independent from the OpenMP code
import importlib

from .util import check_box, check_coords

_distances = {}
_distances["serial"] = importlib.import_module(".c_distances", package="MDemon.lib")
try:
    _distances["openmp"] = importlib.import_module(
        ".c_distances_openmp", package="MDemon.lib"
    )
except ImportError:
    pass

del importlib


def _run(
    funcname: str,
    args: Optional[tuple] = None,
    kwargs: Optional[dict] = None,
    backend: str = "serial",
) -> Callable:
    """Helper function to select a backend function `funcname`."""
    args = args if args is not None else tuple()
    kwargs = kwargs if kwargs is not None else dict()
    backend = backend.lower()
    try:
        func = getattr(_distances[backend], funcname)
    except KeyError:
        errmsg = (
            f"Function {funcname} not available with backend {backend} "
            f"try one of: {_distances.keys()}"
        )
        raise ValueError(errmsg) from None
    return func(*args, **kwargs)


# serial versions are always available (and are typically used within
# the core and topology modules)
from .c_distances import _UINT64_MAX


def _check_result_array(result: Optional[npt.NDArray], shape: tuple) -> npt.NDArray:
    """Check if the result array is ok to use.

    The `result` array must meet the following requirements:
      * Must have a shape equal to `shape`.
      * Its dtype must be ``numpy.float64``.

    Paramaters
    ----------
    result : numpy.ndarray or None
        The result array to check. If `result` is `None``, a newly created
        array of correct shape and dtype ``numpy.float64`` will be returned.
    shape : tuple
        The shape expected for the `result` array.

    Returns
    -------
    result : numpy.ndarray (``dtype=numpy.float64``, ``shape=shape``)
        The input array or a newly created array if the input was ``None``.

    Raises
    ------
    ValueError
        If `result` is of incorrect shape.
    TypeError
        If the dtype of `result` is not ``numpy.float64``.
    """
    if result is None:
        return np.zeros(shape, dtype=np.float64)
    if result.shape != shape:
        raise ValueError(
            f"Result array has incorrect shape, should be {shape}, got "
            f"{result.shape}."
        )
    if result.dtype != np.float64:
        raise TypeError(
            f"Result array must be of type numpy.float64, got {result.dtype}."
        )
    # The following two lines would break a lot of tests. WHY?!
    #    if not coords.flags['C_CONTIGUOUS']:
    #        raise ValueError("{0} is not C-contiguous.".format(desc))
    return result


@check_coords(
    "reference",
    "configuration",
    reduce_result_if_single=False,
    check_lengths_match=False,
    allow_structure=True,
)
def distance_array(
    reference: Union[npt.NDArray, "Structure"],
    configuration: Union[npt.NDArray, "Structure"],
    box: Optional[npt.NDArray] = None,
    result: Optional[npt.NDArray] = None,
    backend: str = "serial",
) -> npt.NDArray:
    """
    Calculate all possible distances between a reference set and another
    configuration.

    This function is derived from :func:`MDAnalysis.lib.distance.distance_array`
    with minimal neccessary modifications.

    If there are ``n`` positions in `reference` and ``m`` positions in
    `configuration`, a distance array of shape ``(n, m)`` will be computed.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    If a 2D numpy array of dtype ``numpy.float64`` with the shape ``(n, m)``
    is provided in `result`, then this preallocated array is filled. This can
    speed up calculations.

    Parameters
    ----------
    reference :numpy.ndarray or :class:`~MDAnalysis.core.groups.AtomGroup`
        Reference coordinate array of shape ``(3,)`` or ``(n, 3)`` (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally). Also
        accepts an :class:`~MDAnalysis.core.groups.AtomGroup`.
    configuration : numpy.ndarray or :class:`~MDAnalysis.core.groups.AtomGroup`
        Configuration coordinate array of shape ``(3,)`` or ``(m, 3)`` (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally). Also
        accepts an :class:`~MDAnalysis.core.groups.AtomGroup`.
    box : array_like, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array which must have the shape ``(n, m)`` and dtype
        ``numpy.float64``.
        Avoids creating the array which saves time when the function
        is called repeatedly.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    d : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n, m)``)
        Array containing the distances ``d[i,j]`` between reference coordinates
        ``i`` and configuration coordinates ``j``.


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
       Now also accepts single coordinates as input.
    .. versionchanged:: 2.3.0
       Can now accept an :class:`~MDAnalysis.core.groups.AtomGroup` as an
       argument in any position and checks inputs using type hinting.
    """
    confnum = configuration.shape[0]
    refnum = reference.shape[0]

    # check resulting array will not overflow UINT64_MAX
    if refnum * confnum > _UINT64_MAX:
        raise ValueError(
            f"Size of resulting array {refnum * confnum} elements"
            " larger than size of maximum integer"
        )

    distances = _check_result_array(result, (refnum, confnum))
    if len(distances) == 0:
        return distances
    if box is not None:
        boxtype, box = check_box(box)
        if boxtype == "ortho":
            _run(
                "calc_distance_array_ortho",
                args=(reference, configuration, box, distances),
                backend=backend,
            )
        else:
            _run(
                "calc_distance_array_triclinic",
                args=(reference, configuration, box, distances),
                backend=backend,
            )
    else:
        _run(
            "calc_distance_array",
            args=(reference, configuration, distances),
            backend=backend,
        )

    return distances


@check_coords("reference", reduce_result_if_single=False, allow_structure=True)
def self_distance_array(
    reference: Union[npt.NDArray, "Structure"],
    box: Optional[npt.NDArray] = None,
    result: Optional[npt.NDArray] = None,
    backend: str = "serial",
) -> npt.NDArray:
    """
    Calculate all possible distances within a configuration `reference`.

    This function is derived from :func:`MDAnalysis.lib.distance.self.distance_array`
    with minimal neccessary modifications.

    If the optional argument `box` is supplied, the minimum image convention is
    applied when calculating distances. Either orthogonal or triclinic boxes are
    supported.

    If a 1D numpy array of dtype ``numpy.float64`` with the shape
    ``(n*(n-1)/2,)`` is provided in `result`, then this preallocated array is
    filled. This can speed up calculations.

    Parameters
    ----------
    reference : numpy.ndarray or :class:`~MDAnalysis.core.groups.AtomGroup`
        Reference coordinate array of shape ``(3,)`` or ``(n, 3)`` (dtype is
        arbitrary, will be converted to ``numpy.float32`` internally). Also
        accepts an :class:`~MDAnalysis.core.groups.AtomGroup`.
    box : array_like, optional
        The unitcell dimensions of the system, which can be orthogonal or
        triclinic and must be provided in the same format as returned by
        :attr:`MDAnalysis.coordinates.timestep.Timestep.dimensions`:
        ``[lx, ly, lz, alpha, beta, gamma]``.
    result : numpy.ndarray, optional
        Preallocated result array which must have the shape ``(n*(n-1)/2,)`` and
        dtype ``numpy.float64``. Avoids creating the array which saves time when
        the function is called repeatedly.
    backend : {'serial', 'OpenMP'}, optional
        Keyword selecting the type of acceleration.

    Returns
    -------
    d : numpy.ndarray (``dtype=numpy.float64``, ``shape=(n*(n-1)/2,)``)
        Array containing the distances ``dist[i,j]`` between reference
        coordinates ``i`` and ``j`` at position ``d[k]``. Loop through ``d``:

        .. code-block:: python

            for i in range(n):
                for j in range(i + 1, n):
                    k += 1
                    dist[i, j] = d[k]


    .. versionchanged:: 0.13.0
       Added *backend* keyword.
    .. versionchanged:: 0.19.0
       Internal dtype conversion of input coordinates to ``numpy.float32``.
    .. versionchanged:: 2.3.0
       Can now accept an :class:`~MDAnalysis.core.groups.AtomGroup` as an
       argument in any position and checks inputs using type hinting.
    """
    refnum = reference.shape[0]
    distnum = refnum * (refnum - 1) // 2
    # check resulting array will not overflow UINT64_MAX
    if distnum > _UINT64_MAX:
        raise ValueError(
            f"Size of resulting array {distnum} elements larger"
            " than size of maximum integer"
        )

    distances = _check_result_array(result, (distnum,))
    if len(distances) == 0:
        return distances
    if box is not None:
        boxtype, box = check_box(box)
        if boxtype == "ortho":
            _run(
                "calc_self_distance_array_ortho",
                args=(reference, box, distances),
                backend=backend,
            )
        else:
            _run(
                "calc_self_distance_array_triclinic",
                args=(reference, box, distances),
                backend=backend,
            )
    else:
        _run("calc_self_distance_array", args=(reference, distances), backend=backend)

    return distances
