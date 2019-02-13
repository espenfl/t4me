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

#include <WildMagic5/Wm5IntpTrilinear3.h>
#include <WildMagic5/Wm5IntpTricubic3.h>
#include <WildMagic5/Wm5IntpAkimaUniform3.h>
#include <WildMagic5/Wm5Math.h>

void wildmagic_execute_interpolation(int *num_points,double *domainx, double *domainy, double *domainz, double *data_temp, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata,int itype);

void wildmagic_gradient_execute_interpolation(int *num_points,double *domainx, double *domainy, double *domainz, double *data_temp, double *ix, double *iy, double *iz, int ip, int num_bands,double *idata, double *igradx, double *igrady, double *igradz, int itype);
