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
