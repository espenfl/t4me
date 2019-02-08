/* Copyright 2016 Espen Flage-Larsen

    This file is part of T4ME.

    T4ME is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    T4ME is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with T4ME.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "skw.h"

int interpolate_skw_interface(double *energies, int num_bands, double *kpoints, int num_kpoints, int *ksampling, double *lattice, double *positions, int *species, int num_atoms, double star_radius_factor, double *ienergies, double *ivelocities, double *ikpoints, int num_ikpoints, int *iksampling);
