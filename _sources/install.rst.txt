Installing
==========

Basic install
-------------
First make sure Spglib is installed

::

   pip install spglib

Then install T4ME by executing the command

::

   pip install T4ME

This will give the user the posibility to calculate the transport coefficients
using integration routines in SciPy. For other integration and interpolation routines
the user needs to follow the following recipe.

Advanced install
----------------
For more advanced functionality (interpolation and other integration routines) the
user should determine which external libraries are needed and install them based
on their respective documentation. Please also fetch the repository from github and work from its
base directory when executing the following commands.

The :file:`setup.py` file assume in its supplied form that the user installs the libraries
in the standard folders, e.g. `$HOME/include` and `$HOME/lib` for the include and library
files, respectively. If other locations are needed, please adapt the :file:`setup.py`
file.

As an example, we want to enable the tetrahedron integration. A Spglib interface needs to be compiled.
This can be build with the included :file:`build_spglib` file.

::

   ./build_spglib

If that was successfull, T4ME can then be built by issuing the following command

::

   pip install .

or

::

   pip install -e .[dev]

Another example. We want to enable SKW interpolation. The SKW routines can be built (assuming Intel MKL is installed) by issuing

::

   ./build_skw

If successfull, T4ME can then be installed by issuing one of the two commands listed above. If other FFT routines are to be used, please modify :file:`skw/Makefile`.
   
All other libraries need to be built externally and linked in.

Upon successfull completion of the installation, T4ME is executed with the command

::

   t4me

An `input` directory is needed which should contain the input files.
