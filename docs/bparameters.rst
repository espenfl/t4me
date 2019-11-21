.. gparameters:

Band structure input parameters
===============================

.. contents::
   :depth: 2
   :local:

Notes about format
------------------
The input files follow normal YAML conventions.
Please inspect the sample file :file:`input/bandparam.yml`.
Even though many parameters have default values if not
specified the user should always run the calculations with
fully specified input files for consistency and reproducibility.

There is one entry per band. If many bands are used one can
specify a range, e.g. Band X-Y: to set the same parameters
for bands X to Y. Or if one would want to set the same
parameters for all bands one should use Band 1-:
This is quite usefull when reading data from a
full-band calculation of some sort.

Remember to use two spaces indent after each Band
entry (before a new Band entry) in order to comply with
the YAML formatting standard.

Also, the parameters should be indented with two spaces
from the Band entry:

| Band 1-2:
|   aparameter: somevalue
|   anotherparameter: someothervalue
|
| Band 3-5:
|   aparameter: somevalue
|   anotherparameter: someothervalue

Make sure all bands are specified. In this example, five
bands was included.

General parameters
------------------

The following parameters are general and does not relate
directly to a specific scattering mechanism etc.

``type``
~~~~~~~~
Determines how to generate the bands if not read. Relevant only
if ``read`` is set to `param`. The following options are
possible:

- `0` - parabolic bands according to the relation

  .. :math:`E=\\frac{\\hbar^2k^2}{2m}`

  where the effective mass :math:`m` is set by ``effmass``.

- `1` - parabolic bands pluss a quartic correction according
  to the relation:

  .. :math:`E\\frac{\\hbar^2k^2}{2m}+ak^4`

  where the effective mass :math:`m` is set by ``effmass`` and
  the correction factor :math:`a` is set by ``a``.

- `2` - Kane types (alpha correction) according to the relation:

  .. :math:`E(1+\\alpha)=\\frac{\\hbar^2k^2}{2m}`

  where the effective mass :math:`m` is set by ``effmass`` and
  the correction factor :math:`\\alpha` is set by ``a``.

No band folding is performed, except for the tight-binding case.
It is important thus to scale the unit cell such that there is
enough band coverage within the requested region of the chemical
potentials pluss the excess needed for the thermal broadening.

``effmass``
~~~~~~~~~~~
The effective mass in units of the free electron mass along
each configured unit vector in reciprocal cell. Use negative
values to generate bands that curve down and vice versa.

Example:
::

   effmass: [-1.0,-1.0,-1.0]

Generates band that for the parabolic case curves down
with an effective mass along each unit vector of the
configured recirprocal cell equal to the free electron mass.

``a``
~~~~~
The correction factor to be applied. See ``type`` for
additional description. Is given along each unit vector
in the configured reciprocal cell similar to the effective mass.

Example:
::

   a: [-100.0,-100.0,-100.0]

Applies a correction factor of -100.0 along each unit vector
direction in the currently configured reciprocal cell.

``e0``
~~~~~~
An energy shift in units of eV. Applies to the current band.

Example:
::

   e0: 0.0

Shift the band with 0.0 eV.

``status``
~~~~~~~~~~
Determines if this is a valence or a conduction band.
The following options are available:

- `v` - valence band
- `c` - conduction band

Example:
::

  status: v

This band is a valence band.

``kshift``
~~~~~~~~~~
Shift the band by a reciprocal vector, otherwise it
is centered at Gamma. Have to be specified in cartesian
coordinates.

Example:
::

   kshift: [0.0,0.0,0.0]

Do not apply any shift to the current band.

``spin_degen``
~~~~~~~~~~~~~~
The spin degeneracy of the current band. The following options
are available:

- `1` - not spin degenerated
- `2` - spin degenerated

Example:
::

  spin_degen: 2

The current band is spin degenerated.

General scattering related parameters
-------------------------------------

In the following the parameters related to the setup of
the scattering mechanisms are given.


``select_scattering``
~~~~~~~~~~~~~~~~~~~~~
Determines which scattering mechnisms to apply for the current
band. Set element to 1 to include
scattering, 0 otherwise.
Currently the following scattering mechanisms have been
implemented (the number indicate array index, starting at 1):

- 1 elastic acoustic phonon scattering from def. pot.
- 2 non-polar optical phonon scattering
- 3 intervalley phonon scattering
- 4 polar optical phonon scattering
- 5 piezoelectric phonon scattering
- 6 ionized impurity scattering (Brooks-Herring)
- 7 ionized impority scattering (Conwell-Weiskopf)
- 8 alloy scattering
- 9-11 empty slots
- 12 constant scattering

If one does not use the analytic (parabolic)
scattering models and instead use the density of
states to generate the scattering rate, then only
the first four and the last have been implemented
(currently only the first and last have been properly
tested)

Example:
::

  select_scattering: [1,0,0,0,0,0,0,0,0,0,0,0]

Apply acoustic-phonon scattering by deformation potential to the
following band.

``explicit_prefact``
~~~~~~~~~~~~~~~~~~~~
Set an explicit prefactor for the relaxation time instead of using
the prefactor from the density of states or parabolic band models.
This behavior is enabled by setting the relevant element to `1`
for the mechanism where one would like to
specify an explicit prefact (constant tau0 is not
included and is set below) for. Make sure that the total units that
come out should be in fs. This is not always so easy to do due to
temperature variations etc. Thus if the user also perform calculations
at different temperatures, please consider that the prefactor usually
change. This option should only be used by experts. If all elements
in the array is `0`, the scattering models based on density of states
or parabolic bands is used.

Example:
::

   explicit_prefact: [0,0,0,0,0,0,0,0,0,0,0]

Disable the use of explicit prefactors.

``explicit_prefact_values``
~~~~~~~~~~~~~~~~~~~~~~~~~~~
The values of the explicit prefactors. Only relevant for the
entries in ``explicit_prefact`` with a value of `1`.
Remember that the units of the relaxation time come out as fs,
including the energy dependency (density of states or
parabolic band). Depending on the model, the prefactor thus
have different units. Also consider that the prefactor usually
has a temperature and effective mass dependence.

Example:
::

   explicit_prefact_values: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                            0.0, 0.0, 0.0, 0.0, 0.0]

All explicit prefactors of the relaxtion time is set to zero.


Acoustic phonon scattering parameters
-------------------------------------
This model uses the linear Debye model.


``d_a``
~~~~~~~
Acoustic deformation potential in units of eV. Remember to
rescale this is the overlap matrix is not one.

Example:
::

   d_a: 10

Use a deformation potential of 10 eV.

``speed_sound``
~~~~~~~~~~~~~~~
The speed of sound. This is the group velocity of the
low energy acoustic branch that is in the Debye model assumed
to be linear. In units of m/s.

Example:
::

   speed_sound: 10000

Use a group velocity of 10000 m/s.

Piezoelectric phonon scattering parameters
------------------------------------------
This model uses the polarization that is set up
due to strain effects to describes acoustic
phonon scattering. Typically important for polar materials.

``p``
~~~~~
The piezoelectric constant in units of
:math:`\mathrm{C}/\mathrm{m}^2`

Example:
::

   p: 0.0

The piezoelectric constant is set to zero.

``isl``
~~~~~~~
The inverse screening length in the Debye formulation in
units of inverse AA.

Example:
::

   isl: 0.0

The inverse screening length is set to zero.

Non-polar optical phonon scattering
-----------------------------------
This model uses the Einstein model of a optical
phonon mode (dispersion assumed to be flat so a
constant value is used for the frequency).

``d_o``
~~~~~~~
The optical deformation potential in units of eV/AA.

Example:
::

   d_o: 35.0

The optical deformation potential is set to 35.0 eV/AA.

``n_o``
~~~~~~~
The occupation number of the optical phonon.

Example:
::

   n_o: 0.0

The occupation number of the optical phonon is set to zero.

``omega_o``
~~~~~~~~~~~
The optical phonon frequency to use from the Einstein model. In
units of THz.

Example:
::

   omega_o: 0.0

The optical phonon frequency is set to zero.

Polar optical phonon scattering
-------------------------------
After the Froelich model. Should be replaced for a more
explicit model in the future.

``epsi``
~~~~~~~~
The permitivity of the electron in units of the vacuum
permitivity.

Example:
::

   epsi: 0.0

The permitivity is set to zero.

``f``
~~~~~
The Froehlich term.

Example:
::

   f: 0.0

The Froechlich term is set to zero.

Intervalley acoustic phonon scattering
--------------------------------------
A model where the electron scatters both of acoustic and
optical phonon modes. E.g. phonons connect two valleys.

``n_vv``
~~~~~~~~
The intervalley phonon occupation number.

Example:
::

   n_vv: 0.0

The intervalley phonon occupation number is set to zero.


``omega_vv``
~~~~~~~~~~~~
The transition frequency in units of THz.

Example:
::

   omega_vv: 0.0

The transition frequency is set to zero.

``etrans``
~~~~~~~~~~
The transition energy between the bottom of the two values. In
units of eV.

Example:
::

   etrans: 0.0

The transition energy is set to zero.

``zf``
~~~~~~
The number of possible final states (final state degeneracy).

Example:
::

   zf: 0.0

The number of final states is set to zero.

``q_energy_trans``
~~~~~~~~~~~~~~~~~~
The scattering vector connecting the two valleys in direct
reciprocal coordinates.

Example:
::

   q_energy_trans: [[0,0,0],[0.5,0.5,0.5]]

The scattering vector is set along the diagonal reciprocal
cell.

Ionized impurity scattering parameters
--------------------------------------
Parameters using either the  Conwell and Weisskopf (CW) or
the Broks and Herring (BH) model to describe ionized
impurity scattering.


``n_i``
~~~~~~~
The density of ionized impurities in units of
:math:`10^{21} \mathrm{cm}^{-3}`. Used for both the
CW and BH model.

Example:
::

   n_i: 0.01

The density of ionized impurities is set to
:math:`10^{19} \mathrm{cm}^{-3}`.

``isl_i``
~~~~~~~~~
The inverse screening length in units of inverse AA. Only
used for the BH model.

Example:
::

   isl_i: 0.3

The inverse screening length is set to 0.3 inverse AA.

``z``
~~~~~
The number of charge units of the impurity. In units of the
electron charge.

Example:
::

   z: 1.0

The charge of the impurity is set to one electron charge.

Alloy scattering parameters
---------------------------
A scattering model for the alloy
:math:`\mathrm{A}_x\mathrm{B}_{1-x}\mathrm{C}`.


``vdiff``
~~~~~~~~~
The atomic potential difference between the species A and
B in eV.

Example:
::

   vdiff: 1.0

The potential difference is set to 1.0 eV.

``alloyconc``
~~~~~~~~~~~~~
The concentration, :math:`x` of the alloy.

Example:

::

   alloyconc: 0.5

The concentration is set to 50%, i.e. 50% of A and 50% of B.

Common scattering parameters
----------------------------

Here follows scattering parameters that are shared between
the different scattering mechnisms.

``eps``
~~~~~~~
The dielectric constant in units of the vacuum value.

Example:
::

   eps: 12.0

The dielectric constant is set to 12.0 times the vacuum value.

``rho``
~~~~~~~
The mass density of the material in
:math:`\mathrm{g}/\mathrm{cm}^3`.

Example:
::

   rho: 2.4

The mass density of the material is set to 2.4
:math:`\mathrm{g}/\mathrm{cm}^3`.

``tau0_c``
~~~~~~~~~~
The value of the constant relaxation time in units of fs.

Example:
::

   tau0_c: 100.0

The constant relaxation time is set at 100.0 fs.

``emmission``
~~~~~~~~~~~~~
Determines if the considered scattering mechnism is by
emmision or absorption. Acoustic phonon scattering includes both
so this is only relevant where scattering of optical phonons
is encountered.

Example:
::

   emission: False

Use absorption, i.e. a phonon is absorbed in the scattering event.
