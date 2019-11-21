############################
T4ME - Transport 4 MatErials
############################

.. parsed-literal::
       _________________  ____  ___    _______________

      /            /   / /   / /   \__/   /          /

     /____    ____/   / /   / /          /   _______/

         /   /   /   /_/   /_/          /   /___

        /   /   /           /   /\_/   /   ____/

       /   /   /_____    __/   /  /   /   /_______

      /   /         /   / /   /  /   /           /

     /___/         /___/ /___/  /___/___________/

Routines to calculate the transport properties of materials
using the linearized Boltzmann Transport Equations (BTE)
in the Relaxtion-Time-Approximation (RTA).

Please go to the `T4ME documentation <https://t4me.readthedocs.io/en/latest/>`_
for more extensive documentation and
information regarding usage (the API documentation is currently not operational).

Features
********
- Modular, easily extendable by users

- Band structures:

  - Generate the band structure from analytic function

    - Parabolic bands
    - Parabolic bands pluss a quartic correction
    - Kane type of bands

  - Read from first-principle codes

    - Interface to VASP is included
    - Interface to read Wannier90 input and output files
      and use these to construct tight binding orbitals using
      PythTB is included.

  - Read from NumPy datafiles

- Scattering properties:

  - Parabolic energy dispersion models:

    - Acoustic phonon scattering from deformations
    - Non-polar optical phonon scattering (not fully tested)
    - Piezoelectric acoustic phonon scattering (not fully tested)
    - Polar optical phonon scattering (not fully tested)
    - Intervalley phonon scattering (not fully tested)
    - Ionized impurity scattering

  - Density of states models:
    - Acoustic phonon scattering from deformations
    - Non-polar optical phonon scattering (not fully tested)
    - Polar optical phonon scattering (not fully tested)
    - Intervalley phonon scattering (not fully tested)

  - Alloy scattering

- Solution of the transport and density of states integrals:

  - Trapezoidal, Simpson and Romberg integration of a static
    input grid
  - Linear tetrahedron method (Spglib needed)
  - Weighed sum method

- Interpolation of the band structure and scattering properties:

  - All routines available in SciPy
  - GeometricTools/WildMagic regular grid routines


Structure
*********

The structure of the program is simple: the main routines
are written in Python utlizing NumPy and SciPy where
necessary. In addition there are calls to external
routines through Cython, particularly the optional libraries.
Only support for Python3 is confirmed.

Contributing and versioning
***************************

Standard Git versioning is utilized. Contributions are welcome,
encouraged and (greatly) appreciated. Please go here:
`T4ME@GitHub <https://github.com/espenfl/t4me>`_

Author
******

Espen Flage-Larsen with funding from the Norwegian
Research Council, Thelma project (228854).

License
*******

This project is licensed under the BSD 3-clause license. Please see
``LICENSE.md`` included in the root folder of T4ME for additional details.
