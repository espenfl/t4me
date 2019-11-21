/* Copyright 2016 Espen Flage-Larsen

   This file is part of T4ME and covered by the BSD 3-clause license.

   You should have received a copy of the BSD 3-clause license
   along with T4ME.  If not, see <https://opensource.org/licenses/BSD-3-Clause/>.

*/

#include <WildMagic5/Wm5IntpTrilinear3.h>
#include <WildMagic5/Wm5IntpTricubic3.h>
#include <WildMagic5/Wm5IntpAkimaUniform3.h>
#include <WildMagic5/Wm5Math.h>

void wildmagic_execute_interpolation(int *num_points,double *domainx, double *domainy, double *domainz, double *data_temp, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata,int itype);

void wildmagic_gradient_execute_interpolation(int *num_points,double *domainx, double *domainy, double *domainz, double *data_temp, double *ix, double *iy, double *iz, int ip, int num_bands,double *idata, double *igradx, double *igrady, double *igradz, int itype);
