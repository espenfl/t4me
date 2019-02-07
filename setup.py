# -*- coding: utf-8 -*-
"""
Install the T4ME python package.

usage: pip install -e .[graphs]
"""

import os
import json
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from os.path import expanduser
import warnings
from distutils.errors import CCompilerError
from distutils.errors import DistutilsExecError
from distutils.errors import DistutilsPlatformError

SETUP_JSON_PATH = os.path.join(os.path.dirname(__file__), 'setup.json')
README_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.rst')
with open(README_PATH, 'r') as readme:
    LONG_DESCRIPTION = readme.read()

# The tailored check for compilers and presence of extensions dependencies is
# heavily motivated from
# http://charlesleifer.com/blog/misadventures-in-python-packaging-optional-c-extensions/
# Thanks Charles!

# Let us try to build the extensions using Cython (do not bother to supply the C files)
build_extensions = True

# Check if Cython is installed
try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    build_extensions = False

if build_extensions:
    import shutil
    import tempfile
    from distutils.ccompiler import new_compiler
    from distutils.sysconfig import customize_compiler

# set home of current user
home = expanduser("~")

# set standard local locations of lib and include directories
locallib = home + "/lib"
localinclude = home + "/include"

class BuildFailure(Exception):
    """Custom exception to indicate a failure to build C extensions."""

    pass

def _have_extension_support(name, libraries, library_dirs, include_dirs):
    # Make a simple code that checks for the presence of the header file
    c_code = ('#include <{name}.h>\n\n'
              'int main(int argc, char **argv) {{ return 0; }}'.format(name=name))
    tmp_dir = tempfile.mkdtemp(prefix='tmp_{name}_'.format(name=name))
    bin_file = os.path.join(tmp_dir, 'test_{name}'.format(name=name))
    src_file = bin_file + '.cpp'
    with open(src_file, 'w') as fh:
        fh.write(c_code)

    compiler = new_compiler()
    for directory in include_dirs:
        compiler.add_include_dir(directory)
    customize_compiler(compiler)
    success = False
    try:
        compiler.link_executable(
            compiler.compile([src_file], output_dir=tmp_dir),
            bin_file,
            libraries=libraries, library_dirs=library_dirs)
    except CCompilerError:
        print('unable to compile {name} C extensions - missing headers?'.format(name=name))
    except DistutilsExecError:
        print('unable to compile {name} C extensions - no C compiler?'.format(name=name))
    except DistutilsPlatformError:
        print('unable to compile {name} C extensions - platform error'.format(name=name))
    else:
        success = True
    #shutil.rmtree(tmp_dir)

    return success

# Now follows the definitions that sets up the locations of the (optional) dependencies
extensions = []
# Spglib (only needed if you want to use the tetrahedron method for integration, otherwise
# we utilize the PyPI installed Spglib)
spglib_lib = locallib
spglib_include = localinclude
spglib_libraries = ["tetrahedron", "symspg", "tetrahedron", "stdc++"]
spglib_include_dirs = [spglib_include]
spglib_library_dirs = [spglib_lib]
spglib_extension_support = False
if build_extensions:
    if _have_extension_support(name="spglib", libraries=spglib_libraries,
                               include_dirs=spglib_include_dirs, library_dirs=spglib_library_dirs):
        spglib_extension_support = True
        import numpy as np
        extensions.append(Extension("spglib_interface", ["spglib_interface/spglib.pyx"],
                                    include_dirs=spglib_include_dirs.extend(np.get_include()),
                                    library_dirs=spglib_library_dirs,
                                    libraries=spglib_libraries,
                                    extra_compile_args=[
                                        "-std=c++11", "-g", "-w", "-fno-omit-frame-pointer",
                                        "-fno-builtin-malloc -fno-builtin-calloc "
                                        "-fno-builtin-realloc -fno-builtin-free"],
                                    extra_link_args=["-g"]))

class _T4MEBuildExt(build_ext):
    """Custom build_ext class that handles the checks for the presence of compilers and headers."""
    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            raise BuildFailure()

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except (CCompilerError, DistutilsExecError, DistutilsPlatformError):
            raise BuildFailure()

def _do_setup(build_extensions, extensions):
    if not build_extensions:
        ext_modules = None
    else:
        ext_modules = extensions

    with open(SETUP_JSON_PATH, 'r') as info:
        SETUP_KWARGS = json.load(info)

    setup(
        packages=find_packages(),
        keywords='vasp, materials, transport, boltzmann',
        long_description=LONG_DESCRIPTION,
        cmdclass={'build_ext': _T4MEBuildExt},
        ext_modules=ext_modules,
        **SETUP_KWARGS)

if build_extensions:
    try:
        _do_setup(build_extensions, extensions)
    except BuildFailure:
        print('#' * 75)
        print('Error compiling C extensions, C extensions will not be built.')
        print('#' * 75)
        _do_setup(False, False)
else:
    _do_setup(False, False)
    
# # Intel MKL (only needed for the SKW routines)
# mkllib = "/opt/intel/mkl/lib/intel64"
# mklinclude = "/opt/intel/mkl/include"
# fftwlib = mkllib
# fftwinclude = mklinclude + "/fftw"
# # GNU GSL (for the closed Fermi integrals)
# gsllib = locallib
# gslinclude = localinclude
# # Einspline
# einsplinelib = locallib
# einsplineinclude = localinclude
# # Cubature
# cubaturelib = locallib
# cubatureinclude = localinclude
# # Wildmagic/GeometricTools
# wildmagiclib = locallib
# wildmagicinclude = localinclude
# # SKW interpolation
# skwlib = "skw"
# skwinclude = "skw"


# And then the extension definitions
# ext = [
#     # Extension("gsl", ["gsl_interface/gsl.pyx"],
#     #           include_dirs=[gslinclude, np.get_include()],
#     #           library_dirs=[gsllib],
#     #           libraries=["gsl", "gslcblas"]),
#     # Extension("einspline", ["einspline_interface/einspline.pyx"],
#     #           include_dirs=[einsplineinclude, np.get_include()],
#     #           library_dirs=[einsplinelib],
#     #           libraries=["einspline"],
#     #           extra_compile_args=["-std=c++11"],
#     #           language="c++"),
#     # Extension("wildmagic", ["wildmagic_interface/wildmagic.pyx"],
#     #           include_dirs=[wildmagicinclude, np.get_include()],
#     #           library_dirs=[wildmagiclib],
#     #           libraries=["Wm5Core", "Wm5Mathematics"],
#     #           language="c++"),
#     # Extension("cubature_wildmagic", ["cubature_wildmagic_interface/cubature_wildmagic.pyx"],
#     #           include_dirs=[cubatureinclude,
#     #                         wildmagicinclude, np.get_include()],
#     #           library_dirs=[cubaturelib, wildmagiclib],
#     #           libraries=["cubature", "Wm5Core", "Wm5Mathematics"],
#     #           extra_compile_args=["-O3", "-w",
#     #                               "-std=c++11"],
#     #           language="c++"),
#     # Extension("skw_interface", ["skw_interface/skw.pyx"],
#     #           include_dirs=[spglibinclude, mklinclude,
#     #                         skwinclude, fftwinclude, np.get_include()],
#     #           library_dirs=[spgliblib, fftwlib, mkllib, skwlib],
#     #           libraries=["stdc++", "mkl_rt", "pthread",
#     #                      "m", "dl", "skw", "symspg", "fftw3xc_intel"],
#     #           extra_compile_args=[
#     #               "-std=c++11"],
#     #           language="c++"),
#     Extension("spglib_interface", ["spglib_interface/spglib.pyx"],
#               include_dirs=[spglibinclude,
#                             "spglib/src", np.get_include()],
#               library_dirs=[spgliblib],
#               libraries=["tetrahedron", "symspg", "tetrahedron", "stdc++"],
#               extra_compile_args=[
#                   "-std=c++11", "-g", "-w", "-fno-omit-frame-pointer",
#                   "-fno-builtin-malloc -fno-builtin-calloc "
#                   "-fno-builtin-realloc -fno-builtin-free"],
#               extra_link_args=["-g"],
#               language="c++"),
# ]
