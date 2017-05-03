Installing
==========

Please consult the :file:`setup.py` file for details on
how to configure the location of the include and libraries.

Only Spglib (in principle the only requirement for minimal
functionality against calculations based on first-principle
data) is included as a submodule. This can be activated by

::

   git submodule init
   git submodule update

And Spglib can be build with the included :file:`build_spglib`

::

   ./build_spglib

It is possible to build Spglib elsewhere, but consult
:file:`build_spglib` and make sure that the
:file:`tetrahedron_method.h` is also inside the include
directory.

All other libraries need to be built externally and linked in.

Then execute

::

       python setup.py build_ext --inplace
       
To build the necessary parts of interfaces used in T4ME. You should
now be able to now execute t4me.py which is
the main driver routine

::

    python t4me.py

This should execute cleanly. However, it might be usefull
to run T4ME from other locations. In order to for this to
work, the main directory of T4ME needs to be added to
your ``PYTHONPATH`` environmental variable. Usually, this can
be done (in bash) by adding the following to your
:file:`.bashrc` or :file:`.bash_profile`

::
   
   export PYTHONPATH=$PYTHONPATH:/home/username/somepath/t4me
