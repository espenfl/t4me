.. cparameters:

Unit cell input parameters
==========================

.. contents::
   :depth: 2
   :local: 

Notes about format
------------------
The input files follow normal YAML conventions.
Please inspect the sample file :file:`input/cellparam.yml`.
Even though many parameters have default values if not
specified the user should always run the calculations with
fully specified input files for consistency and reproducibility.
      
The unit cell
-------------

``a``
~~~~~
The first lattice vector describing the unit cell in AA.

``b``
~~~~~
The second lattice vector describing the unit cell in AA.

``c``
~~~~~
The third lattice vector describing the unit cell in AA.

Example:
::

   a: [5.0,0.0,0.0]
   b: [0.0,5.0,0.0]
   c: [0.0,0.0,5.0]

Generates a unit cell that is 5.0 by 5.0 by 5.0 AA.

Atomic positions
----------------

``direct``
~~~~~~~~~~
Determines if the atomic positions are given in direct
coordinates.

Example:
::
   
   direct: True

The atomic positions are given in direct coordinates

``pos``
~~~~~~~
A list of the atomic positions. In direct or cartersian
coordinates depending on the parameter ``direct``.

Example:
::
   
   pos: [[0.0,0.0,0.0]]

One atom centered at origo.

``atomtypes``
~~~~~~~~~~~~~
The type of atoms as a list in the same order as ``pos``.
Use abbreviations that are standard to the periodic table.
An addition element X is added for unknown types.

Example:
::

   atomtypes: [X]

K-point grid density
--------------------

``ksampling``
~~~~~~~~~~~~~
The k-point sampling along each axis of the configured
reciprocal unit cell.

Example:
::

   ksampling: [15,15,15]

Use 15 by 15 by 15 samples in the full reciprocal unit cell.
