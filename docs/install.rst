Installing
==========

Please consult the :file:`setup.py` file for details on
how to configure the location of the include and libraries.
Then execute::

       python setup.py build_ext --inplace

And you should be able to now execute t4me.py which is
the main driver routine::

    python t4me.py

This should execute cleanly. However, it might be usefull
to run T4ME from other locations. In order to for this to
work, the main directory of T4ME needs to be added to
your ``PYTHONPATH`` environmental variable. Usually, this can
be done (in bash) by adding the following to your
:file:`.bashrc` or :file:`.bash_profile`

::
   
   export PYTHONPATH=$PYTHONPATH:/home/username/somepath/t4me
