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

#include "einspline_interface.hpp"

void einspline_execute_uniform(int *num_points, double *boundaryx, double *boundaryy, double *boundaryz, double *data, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata, double *igradx, double *igrady, double *igradz, int grad, char *bc) {

  // grid for computation (single uniform)
  Ugrid x_grid;
  Ugrid y_grid;
  Ugrid z_grid;
  x_grid.start=boundaryx[0];
  x_grid.end=boundaryx[1];
  x_grid.num=num_points[0];
  y_grid.start=boundaryy[0];
  y_grid.end=boundaryy[1];
  y_grid.num=num_points[1];
  z_grid.start=boundaryz[0];
  z_grid.end=boundaryz[1];
  z_grid.num=num_points[2];

  bc_code boundary=NATURAL;

  if (strcmp(bc,"PERIODIC")==0) {
    boundary=PERIODIC;
  }
  else if (strcmp(bc,"DERIV1")==0) {
    boundary=DERIV1;
  }
  else if (strcmp(bc,"DERIV2")==0) {
    boundary=DERIV2;
  }
  else if (strcmp(bc,"FLAT")==0) {
    boundary=FLAT;
  }
  else if (strcmp(bc,"NATURAL")==0) {
    boundary=NATURAL;
  }
  else if (strcmp(bc,"ANTIPERIODIC")==0) {
    boundary=ANTIPERIODIC;
  }
  else {
    std::cout << "The supplied boundary condition " << bc << " for the Einspline interpolation is not recognized. Setting the condition to NATURAL. Continuing." << std::endl;
  }

  BCtype_d xBC = {boundary,boundary,0.0,0.0};
  BCtype_d yBC = {boundary,boundary,0.0,0.0};
  BCtype_d zBC = {boundary,boundary,0.0,0.0};

  int num_kpoints = num_points[0] * num_points[1] * num_points[2];

  // create the bspline object
  if (grad==0) {
    for (int band=0;band<num_bands;band++) {
      UBspline_3d_d *spline_3d = create_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,
						      (double*)data + band*num_kpoints);

      for (int i=0;i<ip;i++) {
	double value=0.0;
	//eval_UBspline_3d_d(spline_3d, ix[i], iy[i], iz[i], &idata[i+bands*num_bands]);
	eval_UBspline_3d_d(spline_3d, ix[i], iy[i], iz[i], &value);
      	// print this line if we encounter NaN values
	if (std::isnan(value)) {
	  std::cout << "element " << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	idata[band*ip+i]=value;
      }
      // destroy the bspline object
      destroy_Bspline(spline_3d);
    }
  }
  else if (grad==1) {
    double *grad = new double[3];
    for (int band=0;band<num_bands;band++) {
      UBspline_3d_d *spline_3d = create_UBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,
						      (double*)data + band*num_kpoints);
      for (int i=0;i<ip;i++) {
	double value=0.0;
	grad[0]=0.0;
	grad[1]=0.0;
	grad[2]=0.0;
	eval_UBspline_3d_d_vg(spline_3d, ix[i], iy[i], iz[i], &value, grad);
	// print this line if we encounter NaN values
	if (std::isnan(value)) {
	  std::cout << "element " << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	idata[band*ip+i]=value;
	if (std::isnan(grad[0])) {
	  std::cout << "element grad[0]" << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	igradx[band*ip+i]=grad[0];
	if (std::isnan(grad[1])) {
	  std::cout << "element grad[1]" << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	igrady[band*ip+i]=grad[1];
	if (std::isnan(grad[2])) {
	  std::cout << "element grad[2]" << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	igradz[band*ip+i]=grad[2];
      }
      // destroy the bspline object
      destroy_Bspline(spline_3d);
    }
    delete [] grad;
  }
}

void einspline_execute_nonuniform(int *num_points, double *gridx, double *gridy, double *gridz, double *data, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata, double *igradx, double *igrady, double *igradz, int grad, char *bc) {

  // grid for computation (single uniform)
  NUgrid *x_grid = create_general_grid (gridx,
                                        num_points[0]);
  NUgrid *y_grid = create_general_grid (gridy,
                                        num_points[1]);
  NUgrid *z_grid = create_general_grid (gridz,
                                        num_points[2]);

  bc_code boundary=NATURAL;

  if (strcmp(bc,"PERIODIC")==0) {
    boundary=PERIODIC;
  }
  else if (strcmp(bc,"DERIV1")==0) {
    boundary=DERIV1;
  }
  else if (strcmp(bc,"DERIV2")==0) {
    boundary=DERIV2;
  }
  else if (strcmp(bc,"FLAT")==0) {
    boundary=FLAT;
  }
  else if (strcmp(bc,"NATURAL")==0) {
    boundary=NATURAL;
  }
  else if (strcmp(bc,"ANTIPERIODIC")==0) {
    boundary=ANTIPERIODIC;
  }
  else {
    std::cout << "The supplied boundary condition " << bc << " for the Einspline interpolation is not recognized. Setting the condition to NATURAL. Continuing." << std::endl;
  }

  BCtype_d xBC = {boundary,boundary,0.0,0.0};
  BCtype_d yBC = {boundary,boundary,0.0,0.0};
  BCtype_d zBC = {boundary,boundary,0.0,0.0};

  int num_kpoints = num_points[0];

  // create the bspline object
  if (grad==0) {
    for (int band=0;band<num_bands;band++) {
      NUBspline_3d_d *spline_3d = create_NUBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,
                                                        (double*)data + band*num_kpoints);

      for (int i=0;i<ip;i++) {
	double value=0.0;
	//eval_UBspline_3d_d(spline_3d, ix[i], iy[i], iz[i], &idata[i+bands*num_bands]);
	eval_NUBspline_3d_d(spline_3d, ix[i], iy[i], iz[i], &value);
      	// print this line if we encounter NaN values
	if (std::isnan(value)) {
	  std::cout << "element " << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	idata[band*ip+i]=value;
      }
      // destroy the bspline object
      destroy_Bspline(spline_3d);
    }
  }
  else if (grad==1) {
    double *grad = new double[3];
    for (int band=0;band<num_bands;band++) {
      NUBspline_3d_d *spline_3d = create_NUBspline_3d_d(x_grid,y_grid,z_grid,xBC,yBC,zBC,
                                                        (double*)data + band*num_kpoints);
      for (int i=0;i<ip;i++) {
	double value=0.0;
	grad[0]=0.0;
	grad[1]=0.0;
	grad[2]=0.0;
	eval_NUBspline_3d_d_vg(spline_3d, ix[i], iy[i], iz[i], &value, grad);
	// print this line if we encounter NaN values
	if (std::isnan(value)) {
	  std::cout << "element " << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	idata[band*ip+i]=value;
	if (std::isnan(grad[0])) {
	  std::cout << "element grad[0]" << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	igradx[band*ip+i]=grad[0];
	if (std::isnan(grad[1])) {
	  std::cout << "element grad[1]" << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	igrady[band*ip+i]=grad[1];
	if (std::isnan(grad[2])) {
	  std::cout << "element grad[2]" << i << " contains a NaN, coordinate: " << ix[i] << ' ' <<
	    iy[i] << ' ' << iz[i] << std::endl;
	}
	igradz[band*ip+i]=grad[2];
      }
      // destroy the bspline object
      destroy_Bspline(spline_3d);
    }
    delete [] grad;
  }
  // destroy grid
  destroy_grid(x_grid);
  destroy_grid(y_grid);
  destroy_grid(z_grid);
}
