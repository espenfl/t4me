# -*- coding: utf-8 -*-
"""
Install the T4ME python package.

usage: pip install -e .[graphs]
"""

import os
import json
import tempfile
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from os.path import expanduser
import numpy as np
from distutils.errors import CCompilerError
from distutils.errors import DistutilsExecError
from distutils.errors import DistutilsPlatformError
from distutils.command.sdist import sdist as _sdist
from distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler

SETUP_JSON_PATH = os.path.join(os.path.dirname(__file__), 'setup.json')
README_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), 'README.rst')
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
    print(
        "No extensions will be built with Cython as this is not present. Using the supplied Cython compiled files instead."
    )
    build_extensions = False

# set home of current user
home = expanduser("~")

# set standard local locations of lib and include directories
locallib = home + "/lib"
localinclude = home + "/include"

# inspect MKL location
mklroot = os.environ.get('MKLROOT')
if mklroot is not None:
    mkl_include = mklroot + "/include"
    mkl_lib = mklroot + "/lib/intel64"
    fftw_lib = mkl_lib
    fftw_include = mklinclude + "/fftw"


class BuildFailure(Exception):
    """Custom exception to indicate a failure to build C extensions."""

    pass


def _have_extension_support(name, libraries, library_dirs, include_dirs):
    # Make a simple code that checks for the presence of the header file
    c_code = ('#include <{name}.h>\n\n'
              'int main(int argc, char **argv) {{ return 0; }}'.format(
                  name=name))
    tmp_dir = tempfile.mkdtemp(prefix='tmp_{name}_'.format(
        name=name.replace('/', '')))
    bin_file = os.path.join(tmp_dir,
                            'test_{name}'.format(name=name.replace('/', '')))
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
            libraries=libraries,
            library_dirs=library_dirs)
    except CCompilerError:
        print('Unable to compile {name} C extensions - missing headers? '
              'This extension will not be installed.'.format(name=name))
    except DistutilsExecError:
        print('Unable to compile {name} C extensions - no C compiler? '
              'This extension will not be installed.'.format(name=name))
    except DistutilsPlatformError:
        print('Unable to compile {name} C extensions - platform error. '
              'This extension will not be installed'.format(name=name))
    else:
        success = True
    #shutil.rmtree(tmp_dir)

    return success


# Now follows the definitions that sets up the locations of the (optional) dependencies
extensions = []

# Spglib (only needed if you want to use the tetrahedron method for integration, otherwise
# we utilize the PyPI installed Spglib)
spglib_libraries = ["tetrahedron", "symspg", "tetrahedron", "stdc++"]
spglib_include_dirs = [localinclude, "spglib_interface"]
spglib_library_dirs = [locallib]
spglib_sources = ["src/t4me/spglib_interface/spglib_interface.cpp"]

if build_extensions:
    spglib_sources.append("src/t4me/spglib_interface/spglib.pyx")
else:
    spglib_sources.append("src/t4me/spglib_interface/spglib.cpp")

if _have_extension_support(
        name="spglib",
        libraries=spglib_libraries,
        include_dirs=spglib_include_dirs,
        library_dirs=spglib_library_dirs):
    spglib_include_dirs.append(np.get_include())
    extensions.append(
        Extension(
            "t4me.spglib_interface",
            include_dirs=spglib_include_dirs,
            library_dirs=spglib_library_dirs,
            libraries=spglib_libraries,
            sources=spglib_sources,
            extra_compile_args=[
                "-std=c++11", "-w", "-fno-omit-frame-pointer",
                "-fno-builtin-malloc -fno-builtin-calloc "
                "-fno-builtin-realloc -fno-builtin-free"
            ],
            language="c++"))

# GSL (only if you need access to the closed Fermi integrals)
gsl_libraries = ["gsl", "gslcblas"]
gsl_include_dirs = [localinclude]
gsl_library_dirs = [locallib]
gsl_sources = ["src/t4me/gsl_interface/gsl.c"]

if build_extensions:
    gsl_sources = ["src/t4me/gsl_interface/gsl.pyx"]

if _have_extension_support(
        name="gsl/gsl_sf_fermi_dirac",
        libraries=gsl_libraries,
        include_dirs=gsl_include_dirs,
        library_dirs=gsl_library_dirs):
    gsl_include_dirs.append(np.get_include())
    extensions.append(
        Extension(
            "t4me.gsl",
            include_dirs=gsl_include_dirs,
            library_dirs=gsl_library_dirs,
            libraries=gsl_libraries,
            sources=gsl_sources))

# Einspline (only if you are going to use the interpolation routines present here)
einspline_libraries = ["einspline"]
einspline_include_dirs = [localinclude]
einspline_library_dirs = [locallib]
einspline_sources = ["src/t4me/einspline_interface/einspline_interface.cpp"]

if build_extensions:
    einspline_sources.append("src/t4me/einspline_interface/einspline.pyx")
else:
    einspline_sources.append("src/t4me/einspline_interface/einspline.cpp")

if _have_extension_support(
        name="einspline/bspline",
        libraries=einspline_libraries,
        include_dirs=einspline_include_dirs,
        library_dirs=einspline_library_dirs):
    einspline_include_dirs.append(np.get_include())
    extensions.append(
        Extension(
            "t4me.einspline",
            include_dirs=einspline_include_dirs,
            library_dirs=einspline_library_dirs,
            libraries=einspline_libraries,
            sources=einspline_sources,
            extra_compile_args=["-std=c++11"],
            language="c++"))

# GeometricTools (we still use the old package called WildMagic5 until the developer is
# done with the transition)
geometrictools_libraries = ["Wm5Core", "Wm5Mathematics", "pthread", "stdc++"]
geometrictools_include_dirs = [localinclude]
geometrictools_library_dirs = [locallib]
geometrictools_sources = [
    "src/t4me/wildmagic_interface/wildmagic_interface.cpp"
]
if build_extensions:
    geometrictools_sources.append("src/t4me/wildmagic_interface/wildmagic.pyx")
else:
    geometrictools_sources.append("src/t4me/wildmagic_interface/wildmagic.cpp")

if _have_extension_support(
        name="WildMagic5/Wm5Math",
        libraries=geometrictools_libraries,
        include_dirs=geometrictools_include_dirs,
        library_dirs=geometrictools_library_dirs):
    geometrictools_include_dirs.append(np.get_include())
    extensions.append(
        Extension(
            "t4me.wildmagic",
            include_dirs=geometrictools_include_dirs,
            library_dirs=geometrictools_library_dirs,
            libraries=geometrictools_libraries,
            sources=geometrictools_sources,
            extra_compile_args=["-std=c++11"],
            language="c++"))

# GeometricTools and Cubature (offers also cubature integration routines)
cubature_libraries = ["cubature", "m"]
cubature_include_dirs = [localinclude]
cubature_library_dirs = [locallib]
geometrictools_sources = [
    "src/t4me/cubature_wildmagic_interface/cubature_wildmagic_interface.cpp"
]

if build_extensions:
    geometrictools_sources.append(
        "src/t4me/cubature_wildmagic_interface/cubature_wildmagic.pyx")
else:
    geometrictools_sources.append(
        "src/t4me/cubature_wildmagic_interface/cubature_wildmagic.cpp")

if _have_extension_support(name="WildMagic5/Wm5Math", libraries=geometrictools_libraries,
                           include_dirs=geometrictools_include_dirs, library_dirs=geometrictools_library_dirs) and \
    _have_extension_support(name="cubature", libraries=cubature_libraries,
                           include_dirs=cubature_include_dirs, library_dirs=cubature_library_dirs):
    geometrictools_include_dirs = cubature_include_dirs + geometrictools_include_dirs
    geometrictools_library_dirs = geometrictools_library_dirs + cubature_library_dirs
    geometrictools_libraries = geometrictools_libraries + cubature_libraries
    extensions.append(
        Extension(
            "t4me.cubature_wildmagic",
            include_dirs=geometrictools_include_dirs,
            library_dirs=geometrictools_library_dirs,
            libraries=geometrictools_libraries,
            sources=geometrictools_sources,
            extra_compile_args=["-std=c++11"],
            language="c++"))

# SKW interpolation
if mklroot is not None:
    skw_libraries = [
        "stdc++", "mkl_rt", "pthread", "m", "dl", "skw", "symspg",
        "fftw3xc_intel"
    ],
    skw_include_dirs = [localinclude, mklinclude, fftwinclude, "skw"]
    skw_library_dirs = [locallib, fftwlib, mkllib, "skw"]
    skw_sources = ["src/t4me/skw_interface/skw_interface.cpp"]

    if build_extensions:
        skw_sources.append("src/t4me/skw_interface/skw.pyx")
    else:
        skw_sources.append("src/t4me/skw_interface/skw.cpp")

    if _have_extension_support(
            name="skw",
            libraries=skw_libraries,
            include_dirs=skw_include_dirs,
            library_dirs=skw_library_dirs):
        skw_include_dirs.append(np.get_include())
        extensions.append(
            Extension(
                "t4me.skw_interface",
                include_dirs=skw_include_dirs,
                library_dirs=skw_library_dirs,
                libraries=skw_libraries,
                sources=skw_sources,
                language="c++"))
else:
    print("No Intel MKL installed. SKW extension will not be installed.")


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


def _do_setup(extensions):
    if not extensions:
        ext_modules = None
    else:
        ext_modules = extensions

    with open(SETUP_JSON_PATH, 'r') as info:
        SETUP_KWARGS = json.load(info)

    setup(
        packages=find_packages(where="src"),
        long_description=LONG_DESCRIPTION,
        cmdclass={
            'build_ext': _T4MEBuildExt,
            'sdist': sdist_override
        },
        ext_modules=ext_modules,
        **SETUP_KWARGS)


class sdist_override(_sdist):
    # We override the sdist and force Cython to run in order to make sure the supplied
    # files are up to date.
    def run(self):
        try:
            from Cython.Build import cythonize
            cythonize([
                'src/t4me/spglib_interface/spglib.pyx',
                'src/t4me/gsl_interface/gsl.pyx',
                'src/t4me/einspline_interface/einspline.pyx',
                'src/t4me/cubature_wildmagic_interface/cubature_wildmagic.pyx',
                'src/t4me/wildmagic_interface/wildmagic.pyx',
                'src/t4me/skw_interface/skw.pyx'
            ])
            _sdist.run(self)
        except ImportError:
            print(
                "You do not have Cython installed and are thus not allowed to issue sdist."
            )


try:
    _do_setup(extensions)
except BuildFailure:
    print('#' * 75)
    print('Error compiling C extensions, C extensions will not be built.')
    print('#' * 75)
