############################
T4ME - Transport 4 MatErials
############################

Routines to calculate the transport properties of materials
using the linearized Boltzmann Transport Equations (BTE)
in the Relaxtion-Time-Approximation (RTA).

Getting started
***************

Testing. How to get going is here briefly explained.

Prerequisites
=============

In its basic form T4ME only need the following
dependencies: `Spglib <https://atztogo.github.io/spglib/>`_
and `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_.
Additional optional dependencies include:
`Einspline <http://einspline.sourceforge.net/>`_,
`Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_,
`GeometricTools <https://www.geometrictools.com/>`_
(use Wildmagic 5.14).

Installing
==========

Please consult the [setup]{setup.py} file for details on
how to configure the location of the include and libraries.
Then execute::

       python setup.py build_ext --inplace

And you should be able to now execute t4me.py which is
the main driver routine::

    python t4me.py

Structure
=========

The structure of the program is simple: the main routines
are written in Python utlizing NumPy and SciPy where
necessary. In addition there are calls to external
routines through Cython. It should support Python > 2.7 and
Python 3.

The main driver for the program is the `t4me.py` file
located in the main directory. This is used to execute
the program. The parameters and input files are located,
and should be located (for user supplied input files) in
the `input` folder. Here three parameters files should
also be present:

* `param.yml` - contains the main parameters, which
  routines to execute, what type of integration,
  interpolation etc. to perform and what input
  files to use.

* `cellparam.yml` - contains details of the unit cell, its atoms
  and the reciprocal sampling density.

* `bandparam.yml` - contains parameters for band generation, 
  the scattering parameters for each band
  and the parameters to
  use in the tight binding generation
  (if that is needed).

Running tests
*************

Some very simple tests are included located inside the
`tests` folder. They contain the necessary input files
to run the tests and should not be modified. The tests
can be executed by setting `run_tests` to True in the
main parameter file located by default in `input/param.yml`.
Then all tests are executed by issuing::

     python t4me.py

Built With
**********
* `Spglib <https://atztogo.github.io/spglib/>`_
* `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
* `Einspline <http://einspline.sourceforge.net/>`_
* `Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_
* `GeometricTools <https://www.geometrictools.com/>`_

Versioning
**********

Standard Git versioning is utilized.

Author
******
Espen Flage-Larsen with finances from the Norwegian
Research Council, Thelma project (228854).

License
*******

This project is licensed under the GNU GPLv3. Please see
:file:`LICENSE.md` for additional details.
