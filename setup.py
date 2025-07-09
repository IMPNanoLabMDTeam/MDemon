import os
import sys

from Cython.Build import cythonize
from setuptools import Extension, setup

use_cython = True
cython_linetrace = True
annotate_cython = True


class MDExtension(Extension):
    def __init__(self, name, sources, *args, **kwargs):
        self._mda_include_dirs = []
        # don't abspath sources else packaging fails on Windows (issue #3129)
        super(MDExtension, self).__init__(name, sources, *args, **kwargs)


def get_numpy_include():
    # Obtain the numpy include directory. This logic works across numpy
    # versions.
    # setuptools forgets to unset numpy's setup flag and we get a crippled
    # version of it unless we do it ourselves.
    import builtins

    builtins.__NUMPY_SETUP__ = False
    try:
        import numpy as np
    except ImportError:
        print('*** package "numpy" not found ***')
        sys.exit(-1)
    return np.get_include()


include_dirs = [get_numpy_include()]
# Windows automatically handles math library linking
# and will not build MDemon if we try to specify one
if os.name == "nt":
    mathlib = []
else:
    mathlib = ["m"]

source_suffix = ".pyx" if use_cython else ".c"
cpp_source_suffix = ".pyx" if use_cython else ".cpp"

extra_compile_args = []
cpp_extra_compile_args = []
define_macros = []

if cython_linetrace:
    extra_compile_args.append("-DCYTHON_TRACE_NOGIL")
    cpp_extra_compile_args.append("-DCYTHON_TRACE_NOGIL")


distances = MDExtension(
    "MDemon.lib.c_distances",
    ["MDemon/lib/c_distances" + source_suffix],
    include_dirs=include_dirs + ["MDemon/lib/include"],
    libraries=mathlib,
    define_macros=define_macros,
    extra_compile_args=extra_compile_args,
)

source = MDExtension(
    "MDemon.core.source",
    ["MDemon/core/source" + source_suffix],
    include_dirs=include_dirs,
    libraries=mathlib,
    define_macros=define_macros,
    extra_compile_args=extra_compile_args,
)

pre_exts = [distances, source]

cython_generated = []
if use_cython:
    extensions = cythonize(
        pre_exts,
        annotate=annotate_cython,
        compiler_directives={
            "linetrace": cython_linetrace,
            "embedsignature": False,
            "language_level": "3",
        },
        force=True,
    )
    if cython_linetrace:
        print("Cython coverage will be enabled")
    for pre_ext, post_ext in zip(pre_exts, extensions):
        for source in post_ext.sources:
            if source not in pre_ext.sources:
                cython_generated.append(source)

# 只保留ext_modules配置，其他配置已移至pyproject.toml
setup(
    ext_modules=extensions,
)
