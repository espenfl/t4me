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

#include <iostream>
#include <string.h>
#include <einspline/bspline.h>
#include <einspline/nubspline.h>

void einspline_execute_uniform(int *num_points, double *boundaryx, double *boundaryy, double *boundaryz, double *data, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata, double *igradx, double *igrady, double *igradz, int grad, char *bc);

void einspline_execute_nonuniform(int *num_points, double *gridx, double *gridy, double *gridz, double *data, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata, double *igradx, double *igrady, double *igradz, int grad, char *bc);
