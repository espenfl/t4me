Prerequisites
=============

In its basic form T4ME only need the following
dependency:

- `Spglib <https://atztogo.github.io/spglib/>`_.
  This can be installed with

::

   pip install spglib

Additional optional dependencies include:

- `Spglib <https://atztogo.github.io/spglib/>`_,
  A custom interface to Spglib to enable tetrahedron integration.
- `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_,
  Used to solve the analytical Fermi-Dirac integrals.
- `Einspline <http://einspline.sourceforge.net/>`_,
  Used for interpolation of band structure data.
- `Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_,
  Used to integrate the transport and density of states
  using a predetermined accuracy. Used on the fly interpolation
  from the GeometricTools/WildMagic package.
- `ALGLIB <http://www.alglib.net/>`_,
  Used to interpolate the band structure data. Offers a less
  memory intensive RBF method than the one included in SciPy.
- `GeometricTools <https://www.geometrictools.com/>`_
  (use Wildmagic 5.14).
  Used to interpolate band structure data.
