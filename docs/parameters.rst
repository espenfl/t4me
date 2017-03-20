.. parameters:

General input parameters
================

.. contents::
   depth: 2

Dispersion
-----------------------------

The following parameters are related to the energy and velocity
dispersions.

``dispersion_interpolate``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If set to ``True`` the bandstructure is interpolated on a
k-point grid.

``dispersion_interpolate_sampling``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The target k-point sampling when performing interpolation.

::

   dispersion_interpolate_sampling: [45,45,45]

Interpolates the input band structure to a grid density of
45, 45 and 45 k-points along the unit axis of the supplied
k-point grid.

``dispersion_interpolate_method``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Choses which interpolative method to use

::
   
   dispersion_interpolate_method: "wildmagic"

Will for instance use the Wildmagic library.

``dispersion_interpolate_type``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Additional selective layer for the method chosen by
:ref'`dispersion_interpolate_method`.

::
   dispersion_interpolate_type: "akima"

Uses the Akima interpolation in the WildMagic library.

``dispersion_velocities_numdiff``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Use numerical differentiation to calculate the
velocities if they are not present on entry, or/and
use numerical differentiation to extract the
velocities after the dispersions have been
interpolated (used by default for the interpolat
routines that do not support velocity extraction)

::
   dispersion_velocities_numdiff: False

Turns for instance of the numerical difference calculation
of the velocities. In this case please make sure that
the velocities are present on input or that they are
genrated by other means.


``dispersion_write_preinter``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Selects if a line extraction of the band structure is written to
the file :file:`bands` before interpolation. If velocities are present
this is also written to the :file:`velocities`

::

   dispersion_write_preinter: False

Writes the extracted band structure values along a line to file(s).

``dispersion_write_postinter``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Selects if a line extraction of the band structure is written to
the file :file:`bands_inter` after interpolation. If velocities
are present this is also written to the :file:`velocities_inter`

::
   dispersion_write_postinter: False

Does not write the extracted band structure values along a line
to file(s).

``dispersion_write_start``
~~~~~~~~~~~~~~~~~~~~~~~~~~
The start point (in direct coordinates) for the line extraction.

::
   
   dispersion_write_start: [0.0, 0.0, 0.0]

An example start point, here the Gamma point.

``dispersion_write_end``
~~~~~~~~~~~~~~~~~~~~~~~~
The end point (in direct coordinates) for the line extraction.

::
   
   dispersion_write_end: [0.5, 0.0, 0.0]

``num_kpoints_along_line``
~~~~~~~~~~~~~~~~~~~~~~~~~~
How many samples to use along the line to be extracted.

::

   num_kpoints_along_line: 20

Here 20 points is used along the line.

``dispersion_w90_tb_zero_energy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sets the zero energy in the band structure. This parameter is
passed to `zero_energy` in the :func:`model` function in the :class:`w90`
class in PythTB and is used if the Wannier90 interface of PythTB is to be
used to set up the input. This needs to be enabled in the :ref:`param`
parameter. Please consult the
:ref:`PythTB manual <http://physics.rutgers.edu/pythtb/usage.html>`_
for additional details. In units of eV. Usually set to the Fermi level or
the top of the valence band.

::
   
   dispersion_w90_tb_zero_energy:  5.0

Sets it to 5.0 eV and this value is then subtracted from the energies.

``dispersion_w90_tb_min_hopping_norm``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hopping terms with a complex norm less than this value will not be included
in the tight binding model. This parameter is
passed to `min_hopping_norm` in the :func:`model` function in
the :class:`w90` class in PythTB. Please consult the
:ref:`PythTB manual <http://physics.rutgers.edu/pythtb/usage.html>`_
for additional details. In units of eV.

::
   
   dispersion_w90_tb_min_hopping_norm: 0.01

Tight binding hopping parameters with a norm less than 0.01 eV is not included
in the reconstruction of the tight binding model in PythTB.

``dispersion_w90_tb_max_distance``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hopping terms between two sites will be ignored if the distance is larger than
max_distance.
This parameter is passed to `max_distance` in the :func:`model` function in
the :class:`w90` class in PythTB. Please consult the
:ref:`PythTB manual <http://physics.rutgers.edu/pythtb/usage.html>`_
for additional details. In units of AA.

::
   dispersion_w90_tb_max_distance: 4.0

Hopping terms with a distance larger than 4 AA is not included in the
reconstruction of the tight binding model in PythTB.

Electron transport
------------------

The following parameters determines how the transport of electrons
is to be determined.

``transport_calc``
~~~~~~~~~~~~~~~~~~
Determines if the transport calculations are to executed.

::

   transport_calc: True

Calculate the transport properties.

``transport_method``
~~~~~~~~~~~~~~~~~~~~
Selects which mode to use to calculate the transport properties.
Currently three different modes are accepted;

- `closed` The integrals are solved using the closed Fermi-Dirac
  integrals. Only available if the band structure is generated by
  means of analytic models. Only one scattering mechnism can be used
  for each band in this approach.

- `numeric` A numerical integration of the Fermi-Dirac integrals,
  which allows to concatenate different scattering mechanisms for each
  band.

- `numerick` The integrals are solved by integrating over the k-point
  grid or by utilizing the spectral function.

::
   
   transport_method: "numerick"

In this example the transport integrals are solved using the closed
analytical expressions for the Fermi-Dirac integrals.

``transport_integration_method``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Selects which method to use for solving the integral over the k-points.
Only applicable if ``transport_method`` is set to `numerick`.

- `trapz` Use the trapezoidal integration scheme implemented in SciPy
- `simps` Use the Simpson integration scheme implemented in SciPy
- `romberg` Use the Romberg integration scheme implemented in SciPy
- `tetra` Use the linear tetrahedron method
- `smeared` Use the weighted sum approach with a smearing factor
- `cubature` Use the
  `Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_
  integration library together with one of the interpolation routines
  available in the
  `GeometricTools/WildMagic <https://www.geometrictools.com/>`_
  library. Yields the posibility to specify a target accuracy. This
  approach currently only works for cubic, tetragonal and orthorhombic
  unit cell.
  
``transport_integration_spectral_smearing``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Gaussian smearing factor for the weighted sum approach.
In units of eV. Only relevant if ``transport_integration_method``
is set to `smeared`.

::
   
   transport_integration_spectral_smearing: 0.1

Would set it to 0.1 eV.
   
``transport_integration_spectral_density``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The sampling density of the spectral function. Only relevant if
``transport_integration_method`` is set to `tetra` or `smeared`.

::
   
   transport_integration_spectral_density: 1000

An example requesting 1000 samples.
   
``transport_integration_spectral_energy_cutoff``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines the extra padding that is used for the spectral function on
both sides of the requested chemical potential. If multiple chemical 
potentials are requested, the lowest and the highest value is checked and
the range of the energy interval on which the spectral function is
calculated is padded with the specified value. Only relevant if
``transport_integration_method`` is set to `tetra` or `smeared`. In
units of eV.


::

   transport_integration_spectral_energy_cutoff: 1.0

Here, 1.0 eV is subtracted (added) to the smallest (largest) requested
chemical potential.


``transport_interpolate_method``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines which on-the-fly interpolation method is to be used while
performing the Cubature integration. Only relevant if
``transport_integration_method`` is set to `cubature`. Currently
the only option is `wildmagic` which uses the
`GeometricTools/WildMagic <https://www.geometrictools.com/>`_  library.
Which particular interpolation type to use is set with
``transport_interpolate_type``.

::
   
   transport_integration_method: "wildmagic"

Selects the only available method of interpolation during the
Cubature integration.


``transport_interpolate_type``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines which on-the-fly interpolation type to be used while
performing the Cubature integration. Only relevant if
``transport_integration_method`` is set to `cubature`. Currently
the following options are available:

- `trilinear` Using trilinear interpolation
- `tricubic_exact` Using exact tricubic interpolation
- `tricubic_bspline` Using b-splines
- `akima` Using Akima interpolation

Consult the documentation at
`GeometricTools/WildMagic <https://www.geometrictools.com/>`_ for
additional details. Akima is particularly usefull since it is a
special spline interpolation with local character.

::

   transport_interpolate_type: "akima"

Perform on-the-fly Akima interpolation during Cubature integration.

``transport_chempot_min``
~~~~~~~~~~~~~~~~~~~~~~~~~
The minimum chemical potential requested for which the transport
coefficients are calculated. In units of eV.

::
   
   transport_chempot_min: -1.0

Starts the calculation of the transport properties at -1.0 eV.

``transport_chempot_max``
~~~~~~~~~~~~~~~~~~~~~~~~~
The maximum chemical potential requested for which the transport
coefficients are calculated. In units of eV.

::
   
   transport_chempot_max: 1.0

Ends the calculation of the transport properties at 1.0 eV.

``transport_chempot_samples``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The number of chemical potential samples to use between
``transport_chempot_min`` and ``transport_chempot_max``.

::

   transport_chempot_samples: 100

Extract the transport coefficients at 100 points between
``transport_chempot_min`` and ``transport_chempot_max``.

``transport_energycutband``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Bands that reside ``transport_energycutband`` outside
the chemical potential is dropped from the calculation of the
transport coefficients. All k-points
are currently analyzed in order to determine which bands fall inside
the energy range
[``transport_chempot_min``-``transport_energycutband``,``transport_chempot_max``+``transport_energycutband``]
. Units in eV.

::

   transport_energycutband: 1.0

Substract and add 1.0 eV to ``transport_chempot_min`` and
``transport_chempot_max``, respectively. Bands that does not have
any k-point with energy in the range [-2.0 eV, 2.0 eV] is not included
in the calculation of the transport coefficients.
   
``transport_include_bands``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
A list containing specific bands on which to calculate the transport
coefficients. If the list is empty, use all bands within the range set by
:ref:``transport_energycutband``. Band index starts at 1.

::
   
   transport_include_bands: [3, 4, 10]

Calculate the transport coefficients for band 3, 4 and 10. 

``transport_use_analytic_scattering``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines if the analytic spherical scattering models should be used.
They can be applied also to dispersions which are not spherical, but
such an application have to be physically justified.

::

   transport_use_analytic_scattering: False

Use the density-of-states to set up the scattering mechanisms.
   
``transport_use_scattering_ontfly``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines if the scattering values should also be integrated on-the-fly
when performing Cubature integration. Only relevant if
``transport_integration_method`` is set to `cubature`.

::
   
   transport_use_scattering_ontfly: False

Do not use on-the-fly interpolation of the scattering values.

``transport_drop_valence``
~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines if all valence band should be dropped while reading
e.g. external data. Currently only works for the VASP interface.

::

   transport_drop_valence: False

Do not exclude the valence bands during read-in.

``transport_drop_conduction``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Determines if all conduction bands should be dropped while reading
e.g. external data. Currently only works for the VASP interface.

::

   transport_drop_conduction: False

Do not exclude the conduction bands during read-in.

``transport_isotropic``
~~~~~~~~~~~~~~~~~~~~~~~
Only calculate the first element of the transport tensors during
Cubature integration. Only relevant if ``transport_integration_method``
is set to `cubature`

::

   transport_isotropic: False

``temperature_min``
~~~~~~~~~~~~~~~~~~~
The minimum temperature in K.

::

   temperature_min: 100

The minimum temperature is set at 100 K.

``temperature_max``
~~~~~~~~~~~~~~~~~~~
The maximum temperature in K.

::
   
   temperature_max: 700

The maximum temperature is set at 700 K.

``temperature_steps``
~~~~~~~~~~~~~~~~~~~~~
The number of temperature steps from ``temperature_min``
to ``temperature_max``.

::
   
   temperature_steps: 7

In total 7 temperature steps, resulting in temperature
samplings at 100, 200, 300, 400, 500, 600 and 700 K.

``gamma_center``
~~~~~~~~~~~~~~~~
:math:`\\Gamma` centered k-point grids? Anything else is currently
not supported (or tested).

::

   gamma_center: True

Notifies that the k-point grids are :math:`\\Gamma` centered.

``maxeint``
~~~~~~~~~~~
The limites of the dimensionless carrier energy :math:`\\eta`
used for the numerical solution of the Fermi-Dirac integrals.
Only relevant if ``transport_method`` is set to `numerick`.

::
   
   maxeint: 100

Sets the limits of the Fermi-Dirac integrals to 100 :math:`\\eta`.

``occ_cutoff``
~~~~~~~~~~~~~~
The cutoff to use when detecting occupancies. Used for detecting
the valence band maximum, conduction band minimum and then also for
the band gap.

::
   
   occ_cutoff: 1.0e-4

The occupancy cutoff is set at 1.0e-4, which means that states with
an occupancy less than this will be assumed not occupied and vice
versa.

``e_fermi_in_gap``
~~~~~~~~~~~~~~~~~~
Determines if the Fermi level is to be placed in the middle of
the gap.

::

   e_fermi_in_gap: False

Do not place the Fermi level in the middle of the gap.

``e_fermi``
~~~~~~~~~~~
Determine if one should shift the energies to the supplied
Fermi level (usually read in the interface).

::
   e_fermi: True

Shift the energies such that zero is placed at the supplied
Fermi level.


``e_vbm``
~~~~~~~~~
Determines if to set the Fermi level at the valence band
maximum.

::
   
   e_vbm: False

Do not set the Fermi level at the top valence band.

``e_shift``
~~~~~~~~~~~
After all alignments have been performed, perform
this additional shift. Units in eV.

::

   e_shift: 0.0

Sets the additional energy shift to 0 eV.

``cubature_h``
~~~~~~~~~~~~~~
Determines if to use p- or h-cubature for the Cubature integration.
Consult the manual at
`Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_
Only relevant if ``transport_integration_method`` is set to `cubature`.

::
   
   cubature_h: False

Use p-cubature.


``cubature_max_it``
~~~~~~~~~~~~~~~~~~~
The maximum number of iterations while performing Cubature
integration.
Consult the manual at
`Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_
Only relevant if ``transport_integration_method`` is set to `cubature`.

::
   
   cubature_max_it: 0

No maximum limit to the number of iterations (integration stops
when ``cubature_abs_err`` or ``cubature_rel_err`` is reached)

``cubature_abs_err``
~~~~~~~~~~~~~~~~~~~~
The absolute error when the Cubature integration is truncated.
Consult the manual at
`Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_
Only relevant if ``transport_integration_method`` is set to `cubature`.

::
   
   cubature_abs_err: 0.0

The relative error is set at 0.0. If ``cubature_rel_err`` is set
larger than zero, it takes precense.

``cubature_rel_err``
~~~~~~~~~~~~~~~~~~~~
The relative error when the Cubature integration is truncated.
Consult the manual at
`Cubature <http://ab-initio.mit.edu/wiki/index.php/Cubature>`_
Only relevant if ``transport_integration_method`` is set to `cubature`.

::

   cubature_rel_err: 0.01

Truncate the Cubature integration after a relative error of 0.01
is reached. Notice that sometimes, if the transport coefficients are
small (think off-diagonal elements in a isotropic system) it can be
difficult to obtain the requested relative error and the one
enters in practice an infinite loop. Carefully setting
``cubature_max_it`` can alleviate this.

``skw_expansion_factor``
~~~~~~~~~~~~~~~~~~~~~~~~
The expansion factor used in the SKW routine. It is basically
tells how many unit cells that can be used. Only relevant if
``dispersion_interpolate_method`` is set to `skw`.

::
   
   skw_expansion_factor: 5

Use 5 unit cells in each direction. In a second step a sphere is cut
from this volume, thus removing the points in the far corners of
this volume in the interpolation procedure.

``carrier_valence_energy``
~~~~~~~~~~~~~~~~~~~~~~~~~~
The cutoff in which where to interpret the carriers as p-type.
Used in the calculation of the carrier concentration. Units in
eV.

::
   
   carrier_valence_energy: 0.0

Would make sure all carriers at negative energies are interpreted
as p-type.

``carrier_conduction_energy``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The cutoff in which where to interpret the carriers as n-type.
Used in the calculation of the carrier concentration. Units in
eV.

::
   
   carrier_valence_energy: 0.0

Would make sure all carriers at positive energies are interpreted
as n-type.


``carrier_dos_analytick``
~~~~~~~~~~~~~~~~~~~~~~~~~
Determines if the carrier concentration should be recaculated after
being set up with analytical models. Only relevant if the band structure
is generated from analytical models.

::
   
   carrier_dos_analytick: True

Do not recalculate and use the analytical expressions for the carrier
concentration.

``defect_ionization``
~~~~~~~~~~~~~~~~~~~~~

# Use defect ionization?
defect_ionization: False # n/a
# How many donors?
donor_number: 0.0 # 10^21 cm^-3
# Donor degeneration factor:
donor_degen_fact: 0.75 # #
# Donor energy:
donor_energy: 0.0 # eV
# How many acceptors?
acceptor_number: 0.0 # 10^21 cm^-3
# Acceptor degeneration factor:
acceptor_degen_fact: 0.25 # #
# Acceptor energy:
acceptor_energy: 0.0 # eV
# Input data and parameter read in. Consult manual.
read: vasp # n/a
# Filename for the input data etc. Consult manual.
readfile: "" # n/a
# Do you want to use scissor operator?
scissor: False # eV or n/a
# What symprec to use?
symprec: 1.0e-6
# Dump tight binding construction data to std out?
displaytb: False # n/a
# Print output from libraries bundled with T4ME
libinfo: True # n/a
# Only store total scattering rate
# (saves memory), but need to be enabled to see
# the values of each mechanism
onlytotalrate: True
# Paralel mode?
parallel: False
# Run tests? If so, nothing else is done
run_tests: False


