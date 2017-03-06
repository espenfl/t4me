from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy

# THE ONLY DEPENDENCIES NEEDED ARE THE GSL AND
# SPGLIB 1.7.4. THE OTHERS ARE OPTIONAL AND IF NOT
# INSTALLED WILL LIMIT THE CHOICE OF INTEGRATION AND
# INTERPOLATION TECHNIQUES. UNCOMMENT THE OPTIONAL
# IF YOU FEEL YOU WANT TO TEST THESE FEATURES.

home = "/home/flagre"
locallib = home + "/local/lib"
localinclude = home + "/local/include"

# ABSOLUTELY REQUIRED

# Spglib
# this is for a local version of spglib
#spgliblib = locallib
#spglibinclude = localinclude + "/spglib"
# this is for the submodule version (recommended)
spgliblib = "spglib/lib"
spglibinclude = "spglib/include/spglib"

# OPTIONAL
# GNU GSL
gsllib = locallib
# (the gsl .h files use gsl/somefile.h)
gslinclude = localinclude

# Intel MKL (only needed for SKW)
mkllib = "/opt/intel/mkl/lib/intel64"
mklinclude = "/opt/intel/mkl/include"
fftwlib = mkllib
fftwinclude = mklinclude + "/fftw"

# SKW interpolation
skwlib = "skw"
skwinclude = "skw"

# Einspline
einsplinelib = locallib
einsplineinclude = localinclude + "/einspline"

# Cubature
cubaturelib = locallib
cubatureinclude = localinclude + "/cubature"

# Wildmagic
wildmagiclib = "/usr/lib64"
wildmagicinclude = "/usr/include/WildMagic"

# For development purposes
#gptoolslib = locallib
#gptoolsinclude = localinclude + "/gptools"

ext = [
    #        Extension("gsl", ["gsl_interface/gsl.pyx"],
    #                 include_dirs=[gslinclude, numpy.get_include()],
    #                 library_dirs=[gsllib],
    #                 libraries=["gsl", "gslcblas"]),

    #       Extension("einspline", ["einspline_interface/einspline.pyx"],
    #                 include_dirs=[einsplineinclude, numpy.get_include()],
    #                 library_dirs=[einsplinelib],
    #                 libraries=["einspline"],
    #                 extra_compile_args=["-std=c++11"],
    #                 language="c++"),

    #       Extension("wildmagic", ["wildmagic_interface/wildmagic.pyx"],
    #                 include_dirs=[wildmagicinclude, numpy.get_include()],
    #                 library_dirs=[wildmagiclib],
    #                 libraries=["Wm5Core", "Wm5Mathematics"],
    #                 language="c++"),

    #       Extension("cubature_wildmagic", ["cubature_wildmagic_interface/cubature_wildmagic.pyx"],
    #                 include_dirs=[cubatureinclude,
    #                               wildmagicinclude, numpy.get_include()],
    #                 library_dirs=[cubaturelib, wildmagiclib],
    #                 libraries=["cubature", "Wm5Core", "Wm5Mathematics"],
    #                 extra_compile_args=["-O3", "-w",
    #                                     "-std=c++11"],
    #                 language="c++"),

    #       Extension("skw_interface", ["skw_interface/skw.pyx"],
    #                 include_dirs=[spglibinclude, mklinclude,
    #                               skwinclude, fftwinclude, numpy.get_include()],
    #                 library_dirs=[spgliblib, fftwlib, mkllib, skwlib],
    #                 libraries=["stdc++", "mkl_rt", "pthread",
    #                            "m", "dl", "skw", "symspg", "fftw3xc_intel"],
    #                 extra_compile_args=[
    #                     "-std=c++11"],
    #                 language="c++"),
    # special include for tetrahedron_method.c
    # (to be fixed in the future) when this is fully separted
    # in spglib, tetrahedron_method is compiled and linked manually
    # by the compile script in the base directory
    Extension("spglib_interface", ["spglib_interface/spglib.pyx"],
              include_dirs=[spglibinclude,
                            "spglib/src", numpy.get_include()],
              library_dirs=[spgliblib],
              libraries=["tetrahedron", "symspg", "tetrahedron", "stdc++"],
              extra_compile_args=[
        "-std=c++11", "-g", "-w", "-fno-omit-frame-pointer", "-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"],
        extra_link_args=["-g"],
        language="c++"),
    #       Extension("prof", ["prof/prof.pyx"],
    #                 include_dirs=[gptoolsinclude],
    #                 library_dirs=[gptoolslib],
    #                 libraries=["profiler", "tcmalloc"]),
]

setup(name='T4M',
      version='1.0',
      description='',
      author='Espen Flage-Larsen',
      author_email='espen.flage-larsen@sintef.no',
      url='',
      cmdclass={'build_ext': build_ext},
      ext_modules=ext)
