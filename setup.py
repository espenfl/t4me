from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import numpy

# THE ONLY DEPENDENCIES NEEDED ARE THE GSL AND
# SPGLIB 1.7.4. THE OTHERS ARE OPTIONAL AND IF NOT
# INSTALLED WILL LIMIT THE CHOICE OF INTEGRATION AND
# INTERPOLATION TECHNIQUES. UNCOMMENT THE OPTIONAL
# IF YOU FEEL YOU WANT TO TEST THESE FEATURES.

# required (global version)
# USE THIS AND TAILOR IF YOU HAVE USED SU INSTALLB
# globallib = "/usr/lib64"
# globalinclude = "/usr/include"

# spgliblib = globallib
# spglibinclude = globalinclude + "/spglib"
# gsllib = globallib
# gslinclude = globalinclude + "/gsl"

# required (local version)
# USE THIS AND TAILOR IF YOU INSTALLED AS LOCAL USER
home = "/home/flagre"
locallib = home + "/local/lib"
localinclude = home + "/local/include"

spgliblib = locallib
spglibinclude = localinclude + "/spglib"
gsllib = locallib
# (the gsl .h files use gsl/somefile.h)
gslinclude = localinclude

# optional
mkllib = "/opt/intel/mkl/lib/intel64"
mklinclude = "/opt/intel/mkl/include"
fftwlib = mkllib
fftwinclude = mklinclude + "/fftw"

skwlib = "skw"
skwinclude = "skw"

einsplinelib = locallib
einsplineinclude = localinclude + "/einspline"

cubaturelib = locallib
cubatureinclude = localinclude + "/cubature"

wildmagiclib = "/usr/lib64"
wildmagicinclude = "/usr/include/WildMagic"

#gptoolslib = locallib
#gptoolsinclude = localinclude + "/gptools"

ext = [Extension("gsl", ["gsl_interface/gsl.pyx"],
                 include_dirs=[gslinclude, numpy.get_include()],
                 library_dirs=[gsllib],
                 libraries=["gsl", "gslcblas"]),

       Extension("einspline", ["einspline_interface/einspline.pyx"],
                 include_dirs=[einsplineinclude, numpy.get_include()],
                 library_dirs=[einsplinelib],
                 libraries=["einspline"],
                 extra_compile_args=["-std=c++11"],
                 language="c++"),

       Extension("wildmagic", ["wildmagic_interface/wildmagic.pyx"],
                 include_dirs=[wildmagicinclude, numpy.get_include()],
                 library_dirs=[wildmagiclib],
                 libraries=["Wm5Core", "Wm5Mathematics"],
                 language="c++"),

       Extension("cubature_wildmagic", ["cubature_wildmagic_interface/cubature_wildmagic.pyx"],
                 include_dirs=[cubatureinclude,
                               wildmagicinclude, numpy.get_include()],
                 library_dirs=[cubaturelib, wildmagiclib],
                 libraries=["cubature", "Wm5Core", "Wm5Mathematics"],
                 extra_compile_args=["-O3", "-w",
                                     "-std=c++11"],
                 language="c++"),

       Extension("skw_interface", ["skw_interface/skw.pyx"],
                 include_dirs=[spglibinclude, mklinclude,
                               skwinclude, fftwinclude, numpy.get_include()],
                 library_dirs=[spgliblib, fftwlib, mkllib, skwlib],
                 libraries=["stdc++", "mkl_rt", "pthread",
                            "m", "dl", "skw", "symspg", "fftw3xc_intel"],
                 extra_compile_args=[
                     "-std=c++11"],
                 language="c++"),

       Extension("spglib_interface", ["spglib_interface/spglib.pyx"],
                 include_dirs=[spglibinclude, numpy.get_include()],
                 library_dirs=[spgliblib],
                 libraries=["symspg", "stdc++"],
                 extra_compile_args=[
                     "-std=c++11", "-g", "-w", "-fno-omit-frame-pointer", "-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free"],
                 extra_link_args=["-g"],
                 language="c++"),
       #       Extension("prof", ["prof/prof.pyx"],
       #                 include_dirs=[gptoolsinclude],
       #                 library_dirs=[gptoolslib],
       #                 libraries=["profiler", "tcmalloc"]),
       ]

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=ext)
