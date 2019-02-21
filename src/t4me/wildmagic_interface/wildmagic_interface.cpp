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

#include "wildmagic_interface.hpp"
using namespace Wm5;

void wildmagic_execute_interpolation(int *num_points,double *domainx, double *domainy, double *domainz, double *data_temp, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata,int itype) {

  // initialize domain
  int xBound=num_points[0];
  double xMin=domainx[0];
  double xSpacing=(domainx[1]-domainx[0])/(xBound-1);
  int yBound=num_points[1];
  double yMin=domainy[0];
  double ySpacing=(domainy[1]-domainy[0])/(yBound-1);
  int zBound=num_points[2];
  double zMin=domainz[0];
  double zSpacing=(domainz[1]-domainz[0])/(zBound-1);
  int num_kpoints = xBound*yBound*zBound;

  // now we do something that is probably _really_ uncessary, but need to get going...
  double ****data=new double ***[num_bands];
  for (int band=0;band<num_bands;band++) {
    data[band]=new double **[zBound];
    for (int k=0;k<zBound;k++) {
      data[band][k]=new double*[yBound];
      for (int j=0;j<yBound;j++) {
        data[band][k][j]=new double[xBound];
      }
    }
  }
  for (int band=0;band<num_bands;band++) {
    for (int k=0;k<zBound;k++) {
      for (int j=0;j<yBound;j++) {
        for (int i=0;i<xBound;i++) {
          data[band][k][j][i]=data_temp[band*num_kpoints+i*zBound*yBound+j*zBound+k];
        }
      }
    }
  }


  // generate interpolator
  if (itype==0) {
    for (int band=0;band<num_bands;band++) {
      IntpTrilinear3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,data[band]);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (value > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
      }
    }
  }
  else if (itype==1) {
    for (int band=0;band<num_bands;band++) {
      IntpTricubic3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,data[band],true);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (value > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
      }
    }
  }
  else if (itype==2) {
    for (int band=0;band<num_bands;band++) {
      IntpTricubic3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,data[band],false);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (idata[i] > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
      }
    }
  }
  else if (itype==3) {
    for (int band=0;band<num_bands;band++) {
      IntpAkimaUniform3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,data[band]);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (idata[i] > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
      }
    }
  }
  else {
    printf("Error in wildmagic_interface: integration type %d is not defined", itype);
    exit(1);
  }
  // kill data
  for (int band=0;band<num_bands;band++) {
    for (int k=0;k<zBound;k++) {
      for (int j=0;j<yBound;j++) {
	delete[] data[band][k][j];
      }
      delete[] data[band][k];
    }
    delete[] data[band];
  }
  delete[] data;
}

void wildmagic_gradient_execute_interpolation(int *num_points,double *domainx, double *domainy, double *domainz, double *data, double *ix, double *iy, double *iz, int ip, int num_bands, double *idata, double *igradx, double *igrady, double *igradz, int itype) {

  // initialize domain
  int xBound=num_points[0];
  double xMin=domainx[0];
  double xSpacing=(domainx[1]-domainx[0])/(xBound-1);
  int yBound=num_points[1];
  double yMin=domainy[0];
  double ySpacing=(domainy[1]-domainy[0])/(yBound-1);
  int zBound=num_points[2];
  double zMin=domainz[0];
  double zSpacing=(domainz[1]-domainz[0])/(zBound-1);
  int num_kpoints = xBound*yBound*zBound;

  // now we do something that is probably _really_ uncessary, but need to get going...
  double ****tempdata=new double ***[num_bands];
  for (int band=0;band<num_bands;band++) {
    tempdata[band]=new double **[zBound];
    for (int k=0;k<zBound;k++) {
      tempdata[band][k]=new double*[yBound];
      for (int j=0;j<yBound;j++) {
        tempdata[band][k][j]=new double[xBound];
      }
    }
  }
  for (int band=0;band<num_bands;band++) {
    for (int k=0;k<zBound;k++) {
      for (int j=0;j<yBound;j++) {
        for (int i=0;i<xBound;i++) {
          tempdata[band][k][j][i]=data[band*num_kpoints+i*zBound*yBound+j*zBound+k];
        }
      }
    }
  }

  // generate interpolator
  if (itype==0) {
    for (int band=0;band<num_bands;band++) {
      IntpTrilinear3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,tempdata[band]);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (value > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
	// also fetch gradients
	igradx[band*ip+i]=interpolator(1,0,0,ix[i], iy[i], iz[i]);
	igrady[band*ip+i]=interpolator(0,1,0,ix[i], iy[i], iz[i]);
	igradz[band*ip+i]=interpolator(0,0,1,ix[i], iy[i], iz[i]);
      }
    }
  }
  else if (itype==1) {
    for (int band=0;band<num_bands;band++) {
      IntpTricubic3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,tempdata[band],true);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (idata[i] > (std::numeric_limits<double>::max()/10.0)) {
	printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
	// also fetch gradients
	igradx[band*ip+i]=interpolator(1,0,0,ix[i], iy[i], iz[i]);
	igrady[band*ip+i]=interpolator(0,1,0,ix[i], iy[i], iz[i]);
	igradz[band*ip+i]=interpolator(0,0,1,ix[i], iy[i], iz[i]);
      }
    }
  }
  else if (itype==2) {
    for (int band=0;band<num_bands;band++) {
      IntpTricubic3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,tempdata[band],false);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (idata[i] > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
	// also fetch gradients
	igradx[band*ip+i]=interpolator(1,0,0,ix[i], iy[i], iz[i]);
	igrady[band*ip+i]=interpolator(0,1,0,ix[i], iy[i], iz[i]);
	igradz[band*ip+i]=interpolator(0,0,1,ix[i], iy[i], iz[i]);
      }
    }
  }
  else if (itype==3) {
    for (int band=0;band<num_bands;band++) {
      IntpAkimaUniform3<double> interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,ySpacing,zMin,zSpacing,tempdata[band]);
      // fetch points
      for (int i=0;i<ip;i++) {
	double value=interpolator(ix[i], iy[i], iz[i]);
	// print this line if we encounter values above the biggest machine machine represenable
	// number over 10 (WildMagic dumps the largest value if something funny happens)
	if (idata[i] > (std::numeric_limits<double>::max()/10.0)) {
	  printf("element %d did not interpolate correctly, coordinate=%f, %f, %f\n", i, ix[i], iy[i], iz[i]);
	}
	idata[band*ip+i]=value;
	// also fetch gradients
	igradx[band*ip+i]=interpolator(1,0,0,ix[i], iy[i], iz[i]);
	igrady[band*ip+i]=interpolator(0,1,0,ix[i], iy[i], iz[i]);
	igradz[band*ip+i]=interpolator(0,0,1,ix[i], iy[i], iz[i]);
      }
    }
  }
  else {
    printf("Error in wildmagic_interface: integration type %d is not defined", itype);
    exit(1);
  }

  // kill tempdata
  for (int band=0;band<num_bands;band++) {
    for (int k=0;k<zBound;k++) {
      for (int j=0;j<yBound;j++) {
	delete[] tempdata[band][k][j];
      }
      delete[] tempdata[band][k];
    }
    delete[] tempdata[band];
  }
  delete[] tempdata;
}
