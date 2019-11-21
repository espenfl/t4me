/* Copyright 2016 Espen Flage-Larsen

   This file is part of T4ME and covered by the BSD 3-clause license.

   You should have received a copy of the BSD 3-clause license
   along with T4ME.  If not, see <https://opensource.org/licenses/BSD-3-Clause/>.

*/

#include "skw.h"

int interpolate_skw_interface(double *energies, int num_bands, double *kpoints, int num_kpoints, int *ksampling, double *lattice, double *positions, int *species, int num_atoms, double star_radius_factor, double *ienergies, double *ivelocities, double *ikpoints, int num_ikpoints, int *iksampling);
