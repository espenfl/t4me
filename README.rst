############################
T4ME - Transport 4 MatErials
############################

Routines to calculate the transport properties of materials
using the linearized Boltzmann Transport Equations (BTE)
in the Relaxtion-Time-Approximation (RTA).

Features
********

Here comes a condensed list of features.

Getting started
***************

How to get going is here briefly explained. Consult the
documentation for additional details.

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

Please consult the :file:`setup.py` file for details on
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
routines through Cython, particularly the optional libraries.
It should support Python > 2.7 and Python 3.

The main driver for the program is the :file:`t4me.py` file
located in the main directory. This is used to execute
the program.

Input
=====

The program relies on three parameter files written in YAML.
They should all reside inside the `input` directory.

* :file:`input/param.yml` - contains the main parameters, which
  routines to execute, what type of integration,
  interpolation etc. to perform and what input
  files to use.

* :file:`input/cellparam.yml` - contains details of the unit cell, its atoms
  and the reciprocal sampling density.

* :file:`input/bandparam.yml` - contains parameters for band generation, 
  the scattering parameters for each band
  and the parameters to
  use in the tight binding generation
  (if that is needed).
  
These files are documented, please consult respective files
for detailed information of how to set the relevant
parameters.

In addition if external input files be used they should also
be placed in the `input` folder.

Output
======

The transport coefficients are written to the `output`
directory. The log file `info.log` is also written in
this directory and can be monitored during a run to
check the process.

Documentation
*************

Download and install
====================

Here comes details regarding download and installation.

Tutorial
========

Here comes a few tutorials.

Examples
========

Here comes a few specific examples.

Interfaces
==========

Here comes a short description of how to write a custom interface.

Running tests
=============

Some very simple tests are included located inside the
`tests` folder. They contain the necessary input files
to run the tests and should not be modified. The tests
can be executed by setting `run_tests` to True in the
main parameter file located by default in :file:`input/param.yml`.
One can also specify "fast", "slow" and "basic". The former
is the same as setting `run_tests` to True, while the middle
runs more elaborate, but also more time consuming tests. The former
is reserved for usage where one would want to test only the
basic functionality (i.e. `Spglib <https://atztogo.github.io/spglib/>`_
and `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_
support).

After `run_tests` is specified tests are executed by issuing::

     python t4me.py


Contributing and versioning
**********

Standard Git versioning is utilized. Contributions are welcome and
encouraged.

Author
******
Espen Flage-Larsen with finances from the Norwegian
Research Council, Thelma project (228854).

License
*******

This project is licensed under the GNU GPLv3. Please see
:file:`LICENSE.md` for additional details.
