Output
======

The transport coefficients are written to the `output` directory.
The log file :file:`info.log` is also written in this directory
and can be monitored during a run to check the process.

The following files are of importance:

Files
-----

Electrical conductivity
~~~~~~~~~~~~~~~~~~~~~~~
The electrical conductivity can be found in the :file:`sigma`.
In units of :math:`\mathrm{S/m}`. Consult the documentation in
the header of the output file for layout.

Seebeck coefficients
~~~~~~~~~~~~~~~~~~~~
The Seebeck coefficient can be found in the :file:`seebeck`.
In units of :math:`\mu \mathrm{V/K}`. Consult the documentation in
the header of the output file for layout.

Lorenz coefficient
~~~~~~~~~~~~~~~~~~
The Lorenz coefficient can be found in the :file:`lorenz`.
In units of :math:`10^{-8} \mathrm{V^2/K^2}`. Consult the
documentation in the header of the output file for layout.

The electrothermal conductivity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The electrical part of the thermal conductivity can be found
in the :file:`kappae`. In units of :math:`\mathrm{W}/\mathrm{mK}`.

The charge carrier concentration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The charge carrier concentration can be located in the :file:`cc`.
In units of :math:`10^{21} \mathrm{cm}^{-3}`. Consult the documentation
in the header of the output file for layout.

The Hall coefficient
~~~~~~~~~~~~~~~~~~~~
The Hall coefficient (big R) can be located in the :file:`hall`.
In units of :math:`\mathrm{cm}^3/\mathrm{C}`. Consult the documentation
in the header of the output file for layout.

.. warning: NOT YET IMPLEMENTED WHEN FIRST-PRINCIPLE INPUT IS UTILIZED
   (ONLY FILLED WITH BOGUS DATA). ONLY WORKS FOR SPHERICAL
   BANDS AT THE MOMENT.

The relaxation times
~~~~~~~~~~~~~~~~~~~~
The total relaxation time for each band can be found in
the files :file:`scattering_band_n`, for band number `n`. In units of fs.

Visualization
-------------
Data can easily be visualized with Gnuplot.
Currently no automatic visualization is performed.
Data is blocked on temperature.
