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
