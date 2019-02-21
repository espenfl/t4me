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

#include "cubature_wildmagic_interface.hpp"
using namespace Wm5;

double pi=3.1415926535897932;
double betafact=1e5/8.6173324;
double r_factor[12]={0.0,0.0,0.0,1.0,1.0,2.0,2.0,0.0,0.0,0.0,0.0,0.5};

double fetch_combined_scattering_parabolic(double energy,
					   double *tau0,
					   int num_scatterings,
					   double *energy_trans) {

  double tau=0.0;
  for (int i=0;i<num_scatterings;i++) {
    tau+=tau0[i]*pow(fabs(energy+energy_trans[i]),0.5-r_factor[i]);
  }
  tau = 1.0/tau;
  if ((std::isnan(tau)) || (std::isinf(tau))) {
    tau = std::numeric_limits<double>::max();
  }
  return tau;
}

int f_tdf_0_trilinear(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_trilinear_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_0_trilinear_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim, double
				*fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_trilinear(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_trilinear_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    (temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_trilinear_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim,
				double *fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_trilinear(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_trilinear_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_trilinear_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim,
				double *fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}


int f_tdf_0_trilinear_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_trilinear_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_trilinear_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_trilinear_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,
					      p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_trilinear_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    (temp_energy-chempot) / (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_1_trilinear_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_2_trilinear_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,
					      p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_trilinear_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_2_trilinear_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],p->derivl[1],p->derivl[2],x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTrilinear3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],p->derivm[1],p->derivm[2],x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_tricubic(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_tricubic_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_0_tricubic_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim, double
				*fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_tricubic(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_tricubic_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    (temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_tricubic_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim,
				double *fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_tricubic(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_tricubic_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_tricubic_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim,
				double *fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}


int f_tdf_0_tricubic_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_tricubic_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_tricubic_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_tricubic_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,
					      p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_tricubic_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    (temp_energy-chempot) / (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_1_tricubic_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_2_tricubic_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,
					      p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_tricubic_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_2_tricubic_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],p->derivl[1],p->derivl[2],x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpTricubic3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],p->derivm[1],p->derivm[2],x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_akima(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_akima_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_0_akima_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim, double
				*fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_1_akima(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));

  return 0;

}

int f_tdf_1_akima_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    (temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_akima_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim,
				double *fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_akima(unsigned ndim,const double *x, void *params,
		  unsigned fdim, double *fval) {

  struct f_tdf_params * p = (struct f_tdf_params *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_akima_num_scatter(unsigned ndim,const double *x,
			      void *params, unsigned fdim,
			      double *fval) {

  struct f_tdf_params_num_scatter * p =
    (struct f_tdf_params_num_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_akima_const_scatter(unsigned ndim,const double *x,
				void *params, unsigned fdim,
				double *fval) {

  struct f_tdf_params_const_scatter * p =
    (struct f_tdf_params_const_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_l)->operator()(x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->velocity_interpolator_m)->operator()(x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}


int f_tdf_0_akima_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_akima_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_0_akima_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_akima_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,
					      p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_1_akima_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    (temp_energy-chempot) / (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_1_akima_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*(temp_energy-chempot) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_2_akima_deriv(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_tdf_params_deriv * p =
    (struct f_tdf_params_deriv *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=fetch_combined_scattering_parabolic(fabs(temp_energy),
					      p->tau0,
					      p->num_scatterings,
					      p->tau_energy_trans)*
    temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_tdf_2_akima_deriv_num_scatter(unsigned ndim,const double *x,
				    void *params, unsigned fdim,
				    double *fval) {

  struct f_tdf_params_deriv_num_scatter * p =
    (struct f_tdf_params_deriv_num_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],
					 p->derivl[1],
					 p->derivl[2],
					 x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],
					 p->derivm[1],
					 p->derivm[2],
					 x[0],x[1],x[2]);
  double energy_trans=p->energy_trans;
  double temp_scattering=static_cast<IntpAkimaUniform1<double>*>
    (p->tau_interpolator)->operator()(temp_energy+energy_trans);
  double beta=p->beta;
  double chempot=p->chempot;
  if (std::isnan(temp_scattering)) {
    temp_scattering = std::numeric_limits<double>::max();
  }
  fval[0]=temp_scattering*temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;
}

int f_tdf_2_akima_deriv_const_scatter(unsigned ndim,const double *x,
				      void *params, unsigned fdim,
				      double *fval) {

  struct f_tdf_params_deriv_const_scatter * p =
    (struct f_tdf_params_deriv_const_scatter *)params;
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(x[0],x[1],x[2]);
  double temp_velocity_l=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivl[0],p->derivl[1],p->derivl[2],x[0],x[1],x[2]);
  double temp_velocity_m=static_cast<IntpAkimaUniform3<double>*>
    (p->energy_interpolator)->operator()(p->derivm[0],p->derivm[1],p->derivm[2],x[0],x[1],x[2]);
  double beta=p->beta;
  double chempot=p->chempot;
  fval[0]=temp_velocity_l*temp_velocity_m*
    pow(temp_energy-chempot,2.0) /
    (1+cosh((temp_energy-chempot)*beta));
  return 0;

}

int f_akima(unsigned ndim,const double *x, void *params,
	    unsigned fdim, double *fval) {

  fval[0]=static_cast<IntpAkimaUniform3<double>*>
    (params)->operator()(x[0],x[1],x[2]);
  return 0;

}

int f_akima_gaussian(unsigned ndim,const double *x,
		     void *params, unsigned fdim,
		     double *fval) {

  struct f_gaussian_params * p =
    (struct f_gaussian_params *)params;
  // fetch interpolated value
  double temp_energy=static_cast<IntpAkimaUniform3<double>*>
    (p->interpolator)->operator()(x[0],x[1],x[2]);
  double sigma=p->sigma;
  // apply gaussian smearing and return
  fval[0]=exp(-0.5*pow((p->energy_ref-temp_energy)/sigma,2.0)) /
    (sigma*sqrt(2*pi));
  return 0;

}

int f_tricubic(unsigned ndim,const double *x,
	       void *params, unsigned fdim,
	       double *fval) {

  fval[0]=static_cast<IntpTricubic3<double>*>
    (params)->operator()(x[0],x[1],x[2]);
  return 0;

}

int f_tricubic_gaussian(unsigned ndim,const double *x,
			void *params, unsigned fdim,
			double *fval) {

  struct f_gaussian_params * p =
    (struct f_gaussian_params *)params;
  // fetch interpolated value
  double temp_energy=static_cast<IntpTricubic3<double>*>
    (p->interpolator)->operator()(x[0],x[1],x[2]);
  double sigma=p->sigma;
  // apply gaussian smearing and return
  fval[0]=exp(-0.5*pow((p->energy_ref-temp_energy)/sigma,2.0)) /
    (sigma*sqrt(2*pi));
  return 0;

}

int f_trilinear(unsigned ndim,const double *x,
		void *params, unsigned fdim,
		double *fval) {

  fval[0]=static_cast<IntpTrilinear3<double>*>
    (params)->operator()(x[0],x[1],x[2]);
  return 0;

}

int f_trilinear_gaussian(unsigned ndim,const double *x,
			 void *params, unsigned fdim, double *fval) {

  struct f_gaussian_params * p =
    (struct f_gaussian_params *)params;
  // fetch interpolated value
  double temp_energy=static_cast<IntpTrilinear3<double>*>
    (p->interpolator)->operator()(x[0],x[1],x[2]);
  double sigma=p->sigma;
  // apply gaussian smearing and return
  fval[0]=exp(-0.5*pow((p->energy_ref-temp_energy)/sigma,2.0)) /
    (sigma*sqrt(2*pi));
  return 0;

}

int gaussian(int xBound, int yBound, int zBound,
	     double ***unsmeared, double ***smeared,
	     double energy_ref, double sigma) {

  for (int k=0;k<zBound;k++) {
    for (int j=0;j<yBound;j++) {
      for (int i=0;i<xBound;i++) {
	smeared[k][j][i]=exp(-0.5*pow((energy_ref-unsmeared[k][j][i]) /
				      sigma,2.0))/(sigma*sqrt(2*pi));
      }
    }
  }
  return 0;

}

void cubature_wildmagic_execute_integration(int *num_points,
					    double *domainx,
					    double *domainy,
					    double *domainz,
					    double *data,
					    double &value,
					    double &error,
					    int itype) {

  // initialize domain
  int xBound=num_points[0];
  double xMin=domainx[0];
  double xSpacing=(domainx[1]-domainx[0])/(xBound-1);
  double xMax=xMin+xSpacing*(xBound-1);
  int yBound=num_points[1];
  double yMin=domainy[0];
  double ySpacing=(domainy[1]-domainy[0])/(yBound-1);
  double yMax=yMin+ySpacing*(yBound-1);
  int zBound=num_points[2];
  double zMin=domainz[0];
  double zSpacing=(domainz[1]-domainz[0])/(zBound-1);
  double zMax=zMin+zSpacing*(zBound-1);

  // now we do something that is probably _really_ uncessary, but need to get going...
  double ***tempdata=new double **[zBound];
  for (int k=0;k<zBound;k++) {
    tempdata[k]=new double*[yBound];
    for (int j=0;j<yBound;j++) {
      tempdata[k][j]=new double[xBound];
    }
  }
  for (int k=0;k<zBound;k++) {
    for (int j=0;j<yBound;j++) {
      for (int i=0;i<xBound;i++) {
	tempdata[k][j][i]=data[i*zBound*yBound+j*zBound+k];
      }
    }
  }
  double vmin[3]={xMin,yMin,zMin};
  double vmax[3]={xMax,yMax,zMax};
  unsigned fdim=1;
  unsigned dim=3;
  unsigned max_it=0;
  double abs_err=0.0;
  double rel_err=1e-3;
  error_norm norm=ERROR_INDIVIDUAL;

  // generate interpolator and integrate with cubature
  if (itype==0) {
    IntpTrilinear3<double> interpolator(xBound,yBound,zBound,xMin,
					xSpacing,yMin,ySpacing,
					zMin,zSpacing,tempdata);
    hcubature(fdim,f_trilinear,&interpolator,dim,vmin,vmax,max_it,
	      abs_err,rel_err,norm,&value,&error);
  }
  else if (itype==1) {
    IntpTricubic3<double> interpolator(xBound,yBound,zBound,xMin,
				       xSpacing,yMin,ySpacing,
				       zMin,zSpacing,tempdata,true);
    hcubature(fdim,f_tricubic,&interpolator,dim,vmin,vmax,max_it,
	      abs_err,rel_err,norm,&value,&error);
  }
  else if (itype==2) {
    IntpTricubic3<double> interpolator(xBound,yBound,zBound,xMin,
				       xSpacing,yMin,ySpacing,
				       zMin,zSpacing,tempdata,false);
    hcubature(fdim,f_tricubic,&interpolator,dim,vmin,vmax,max_it,
	      abs_err,rel_err,norm,&value,&error);
  }
  else if (itype==3) {
    IntpAkimaUniform3<double> interpolator(xBound,yBound,zBound,xMin,
					   xSpacing,yMin,ySpacing,
					   zMin,zSpacing,tempdata);
    hcubature(fdim,f_akima,&interpolator,dim,vmin,vmax,max_it,
	      abs_err,rel_err,norm,&value,&error);
  }
  else {
    printf("Error in wildmagic_interface: integration type %d is \
not defined", itype);
    exit(1);
  }
  // kill tempdata
  for (int k=0;k<zBound;k++) {
    for (int j=0;j<yBound;j++) {
      delete[] tempdata[k][j];
    }
    delete[] tempdata[k];
  }
  delete[] tempdata;
}

void calc_transport_tensors_cubature(int *num_points,
				     double *domainx,
				     double *domainy,
				     double *domainz,
				     double *energy, int num_bands,
				     double *velocity,
				     double *tau,
				     double *energy_samples,
				     int *spin_degen,
				     int num_scatterings,
				     int num_energy_samples,
				     double *tau_energy_trans,
				     double *chempot, int num_chempot,
				     double *temperatures,
				     int num_temps,
				     int num_kpoints_ibz,
				     int *k_sampling,
				     double *rec_basis,
				     int tensor_type,
				     int interpol_type,
				     int integral_type,
				     int max_it, double abs_err,
				     double rel_err,
                                     int iso,
                                     int info,
				     double *cond,
				     double *seebeck,
				     double *lorenz,
				     int gen_velocities,
				     int scattering_type) {

  bool tetra=true;
  bool h=false;
  if (integral_type == 1) {
    h=true;
  }

  double g=7.7480917310;
  double hbar=6.582119514;
  double hbar_c=197.3269788;
  double kb=8.6173303;
  double sigmaunit=1e11*g/(16*pow(pi,2.0)*hbar*kb);
  double seebeckunit=1e6;
  double lorenzunit=1e8;

  // wildmagic interpolation
  if ((interpol_type >= 0) && (interpol_type < 4)) {
    // initialize domain
    int xBound=num_points[0];
    double xMin=domainx[0];
    double xSpacing=(domainx[1]-domainx[0])/(xBound-1);
    double xMax=xMin+xSpacing*(xBound-1);
    int yBound=num_points[1];
    double yMin=domainy[0];
    double ySpacing=(domainy[1]-domainy[0])/(yBound-1);
    double yMax=yMin+ySpacing*(yBound-1);
    int zBound=num_points[2];
    double zMin=domainz[0];
    double zSpacing=(domainz[1]-domainz[0])/(zBound-1);
    double zMax=zMin+zSpacing*(zBound-1);
    int tot_k_points=xBound*yBound*zBound;


    if (!gen_velocities) {
      if (info) {
        std::cout << "Cubature/Wildmagic calculations of the transport tensors with supplied velocities";
      }
    }
    else {
      if (info) {
        std::cout << "Cubature/Wildmagic calculations of the transport tensors, extracting the velocities";
      }
    }
    if (info) {
      if (h) {
        if (interpol_type == 0) {
          std::cout << ", performing trilinear interpolation and hcubature integration";
        }
        if ((interpol_type == 1) || (interpol_type == 2)) {
          std::cout << ", performing tricubic interpolation and hcubature integration";
        }
        if (interpol_type == 3) {
          std::cout << ", performing Akima interpolation and hcubature integration";
        }
      }
      else {
        if (interpol_type == 0) {
          std::cout << ", performing trilinear interpolation and pcubature integration";
        }
        if ((interpol_type == 1) || (interpol_type == 2)) {
          std::cout << ", performing tricubic interpolation and pcubature integration";
        }
        if (interpol_type == 3) {
          std::cout << ", performing Akima interpolation and pcubature integration";
        }
      }
    }

    // now we do something that is probably _really_ uncessary, but need to get going..

    // actually it might turn out to be necessary since there are so many different
    // standards for passing multidimensional data around...shame!
    // here we need to comply with Cython AND Wildmagic, so vectors are unfortunately
    // out of the question
    double ****energy_data=new double ***[num_bands];
    double *****velocity_data=new double ****[num_bands];
    for (int band=0;band<num_bands;band++) {
      energy_data[band]=new double **[zBound];
      velocity_data[band]=new double ***[3];
      for (int dir=0;dir<3;dir++) {
	velocity_data[band][dir]=new double **[zBound];
	for (int k=0;k<zBound;k++) {
	  if (dir==0) {
	    energy_data[band][k]=new double*[yBound];
	  }
	  velocity_data[band][dir][k]=new double*[yBound];
	  for (int j=0;j<yBound;j++) {
	    if (dir==0) {
	      energy_data[band][k][j]=new double[xBound];
	    }
	    velocity_data[band][dir][k][j]=new double[xBound];
	  }
	}
      }
    }
    for (int band=0;band<num_bands;band++) {
      for (int dir=0;dir<3;dir++) {
	for (int k=0;k<zBound;k++) {
	  for (int j=0;j<yBound;j++) {
	    for (int i=0;i<xBound;i++) {
	      if (dir==0) {
		energy_data[band][k][j][i]=
		  energy[band*tot_k_points+
			 i*zBound*yBound+j*zBound+k];
	      }
	      velocity_data[band][dir][k][j][i]=
		velocity[band*3*tot_k_points+dir*tot_k_points+
			 i*zBound*yBound+j*zBound+k];
	    }
	  }
	}
      }
    }

    // some more stuff on the stack
    double vmin[3]={xMin,yMin,zMin};
    double vmax[3]={xMax,yMax,zMax};
    unsigned fdim=1;
    unsigned dim=3;
    error_norm norm=ERROR_INDIVIDUAL;
    double value=0.0;
    double error=0.0;

    std::vector<std::vector<std::vector<std::vector<double> > > >
      sigma(num_temps, std::vector<std::vector<std::vector<double> > >
	    (num_chempot,std::vector<std::vector<double> >
	     (3, std::vector<double> (3))));
    std::vector<std::vector<std::vector<std::vector<double> > > >
      chi(num_temps, std::vector<std::vector<std::vector<double> > >
	  (num_chempot,std::vector<std::vector<double> >
	   (3, std::vector<double> (3))));
    std::vector<std::vector<std::vector<std::vector<double> > > >
      kappa(num_temps, std::vector<std::vector<std::vector<double> > >
	    (num_chempot,std::vector<std::vector<double> >
	     (3, std::vector<double> (3))));

    int dir2_limit;
    // non-parabolic scattering
    if (scattering_type==0) {
      if (info) {
        std::cout << ", numerical on-the-fly interpolation of the relaxation time" << std::endl;
      }
      double ***tau_data=new double **[num_temps];
      for (int temp=0;temp<num_temps;temp++) {
	tau_data[temp]=new double *[num_bands];
	for (int band=0;band<num_bands;band++) {
	  tau_data[temp][band]=new double[num_energy_samples];
	  for (int nrgy=0;nrgy<num_energy_samples;nrgy++) {
	    tau_data[temp][band][nrgy]=
	      tau[temp*num_energy_samples*num_bands+
		  band*num_energy_samples+nrgy];
	  }
	}
      }

      double energy_samples_min=energy_samples[0];
      double energy_samples_spacing=energy_samples[1]-energy_samples_min;

      double tau_energy_trans_data[num_bands];
      for (int band=0;band<num_bands;band++) {
	tau_energy_trans_data[band]=0.0;
        for (int scattering=0;scattering<num_scatterings;scattering++) {
          tau_energy_trans_data[band] += tau_energy_trans[band*num_scatterings+scattering];
        }
      }

      // loop temperature
      for (int temp=0;temp<num_temps;temp++) {
	double beta = betafact/temperatures[temp];
	// loop chemical potential
	for (int chemp=0;chemp<num_chempot;chemp++) {
	  // loop bands
	  for (int band=0;band<num_bands;band++) {
            if (info) {
              std::cout << "Cubature/Wildmagic: numerical tabulated scattering at " << temperatures[temp] << "K, chempot: " << chempot[chemp] << " eV" << std::endl;
            }
	    if (interpol_type == 0) {
	      IntpAkimaUniform1<double>
		tau_interpolator(num_energy_samples,energy_samples_min,
				 energy_samples_spacing,
				 tau_data[temp][band]);
	      IntpTrilinear3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,
				    xSpacing,yMin,ySpacing,
				    zMin,zSpacing,energy_data[band]);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTrilinear3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1]);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
                    IntpTrilinear3<double>
                      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2]);

		    struct f_tdf_params_num_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,&tau_interpolator,
		       beta,chempot[chemp],tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_trilinear_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_trilinear_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_trilinear_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_trilinear_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_trilinear_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_trilinear_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_num_scatter p=
		      {&energy_interpolator,&tau_interpolator,beta,
		       chempot[chemp],tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_trilinear_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_trilinear_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_trilinear_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_trilinear_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_trilinear_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_trilinear_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	    if (interpol_type == 1) {
	      IntpAkimaUniform1<double>
		tau_interpolator(num_energy_samples,energy_samples_min,
				 energy_samples_spacing,
				 tau_data[temp][band]);
	      IntpTricubic3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,
				    xSpacing,yMin,ySpacing,
				    zMin,zSpacing,energy_data[band], true);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTricubic3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1], true);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpTricubic3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2], true);

		    struct f_tdf_params_num_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,&tau_interpolator,
		       beta,chempot[chemp],tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
                    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_num_scatter p=
		      {&energy_interpolator,&tau_interpolator,beta,
		       chempot[chemp],tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	    if (interpol_type == 2) {
	      IntpAkimaUniform1<double>
		tau_interpolator(num_energy_samples,energy_samples_min,
				 energy_samples_spacing,
				 tau_data[temp][band]);
	      IntpTricubic3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,
				    xSpacing,yMin,ySpacing,
				    zMin,zSpacing,energy_data[band], false);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTricubic3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1], false);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpTricubic3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2], false);

		    struct f_tdf_params_num_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,&tau_interpolator,
		       beta,chempot[chemp],tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_num_scatter p=
		      {&energy_interpolator,&tau_interpolator,beta,
		       chempot[chemp],tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	    if (interpol_type == 3) {
	      IntpAkimaUniform1<double>
		tau_interpolator(num_energy_samples,energy_samples_min,
				 energy_samples_spacing,
				 tau_data[temp][band]);
	      IntpAkimaUniform3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,
				    xSpacing,yMin,ySpacing,
				    zMin,zSpacing,energy_data[band]);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpAkimaUniform3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1]);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpAkimaUniform3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2]);

		    struct f_tdf_params_num_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,&tau_interpolator,
		       beta,chempot[chemp],tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_akima_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
                    else {
		      pcubature(fdim,f_tdf_0_akima_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_akima_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_akima_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_akima_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_akima_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_num_scatter p=
		      {&energy_interpolator,&tau_interpolator,beta,
		       chempot[chemp],tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_akima_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_akima_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_akima_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_akima_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_akima_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_akima_deriv_num_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	  }
	  std::vector<std::vector<double> > sigma_inv(3,std::vector<double> (3));
          std::vector<std::vector<double> > seebeck_temp(3,std::vector<double> (3));
          std::vector<std::vector<double> > lorenz_highdeg(3,std::vector<double> (3));
          std::vector<std::vector<double> > lorenz_cor(3,std::vector<double> (3));
          std::vector<std::vector<double> > tmp(3,std::vector<double> (3));
          // invert sigma
	  invert_3x3_matrix(sigma[temp][chemp],sigma_inv);
          // now calculate the seebeck tensors and add units
          // also store the conductivity
          multiply_3x3_matrix(sigma_inv,chi[temp][chemp],seebeck_temp);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
	      cond[temp*num_chempot*9+9*chemp+3*dir1+dir2]=
		sigmaunit*sigma[temp][chemp][dir1][dir2]/temperatures[temp];
	      seebeck[temp*num_chempot*9+9*chemp+3*dir1+dir2]=-
		seebeckunit*seebeck_temp[dir1][dir2]/temperatures[temp];
            }
          }
          // now calculate the lorenz tensor and add units
          // first store the highly degenerate part (first term)
          multiply_3x3_matrix(kappa[temp][chemp],sigma_inv,lorenz_highdeg);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
	      lorenz[temp*num_chempot*9+9*chemp+3*dir1+dir2]=
		lorenzunit*lorenz_highdeg[dir1][dir2]/pow(temperatures[temp],2.0);
            }
          }
          // now calculate the correction term and subtract
          multiply_3x3_matrix(chi[temp][chemp],seebeck_temp,tmp);
          multiply_3x3_matrix(tmp,sigma_inv,lorenz_cor);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
              lorenz[temp*num_chempot*9+9*chemp+3*dir1+dir2]-=
                lorenzunit*lorenz_cor[dir1][dir2]/pow(temperatures[temp],2.0);
            }
          }
          // add symmetric values in the tensor 21=12 etc.
          cond[temp*num_chempot*9+9*chemp+3] = cond[temp*num_chempot*9+9*chemp+1];
          cond[temp*num_chempot*9+9*chemp+6] = cond[temp*num_chempot*9+9*chemp+2];
          cond[temp*num_chempot*9+9*chemp+7] = cond[temp*num_chempot*9+9*chemp+5];
          seebeck[temp*num_chempot*9+9*chemp+3] = seebeck[temp*num_chempot*9+9*chemp+1];
          seebeck[temp*num_chempot*9+9*chemp+6] = seebeck[temp*num_chempot*9+9*chemp+2];
          seebeck[temp*num_chempot*9+9*chemp+7] = seebeck[temp*num_chempot*9+9*chemp+5];
          lorenz[temp*num_chempot*9+9*chemp+3] = lorenz[temp*num_chempot*9+9*chemp+1];
          lorenz[temp*num_chempot*9+9*chemp+6] = lorenz[temp*num_chempot*9+9*chemp+2];
          lorenz[temp*num_chempot*9+9*chemp+7] = lorenz[temp*num_chempot*9+9*chemp+5];
	}
      }
      // clean up some mess
      for (int temp=0;temp<num_temps;temp++) {
        for (int band=0;band<num_bands;band++) {
          delete[] tau_data[temp][band];
        }
        delete[] tau_data[temp];
      }
      delete[] tau_data;
    }
    // analytic scattering models (parabolic bands)
    else if (scattering_type==1) {
      if (info) {
        std::cout << ", parabolic scattering models" << std::endl;
      }
      double tau0_data[num_scatterings];
      double **tau_energy_trans_data=new double *[num_bands];
      for (int band=0;band<num_bands;band++) {
	tau_energy_trans_data[band]=new double [num_scatterings];
        for (int scattering=0;scattering<num_scatterings;scattering++) {
          tau_energy_trans_data[band][scattering] = tau_energy_trans[band*num_scatterings+scattering];
        }
      }
      // loop temperature
      for (int temp=0;temp<num_temps;temp++) {
	double beta = betafact/temperatures[temp];
	// loop eta
	for (int chemp=0;chemp<num_chempot;chemp++) {
          if (info) {
            std::cout << "Cubature/Wildmagic: parabolic scattering at " << temperatures[temp] << " K, chempot: " << chempot[chemp] << " eV" << std::endl;
          }
	  // loop bands
	  for (int band=0;band<num_bands;band++) {
	    // fetch scattering
	    for (int scattering=0;
		 scattering<num_scatterings;scattering++) {
	      tau0_data[scattering]=
		tau[temp*num_scatterings*num_bands+
		    num_scatterings*band+scattering];
	    }
	    if (interpol_type == 0) {
	      IntpTrilinear3<double> energy_interpolator(xBound,yBound,zBound,xMin,
							 xSpacing,yMin,ySpacing,
							 zMin,zSpacing,energy_data[band]);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTrilinear3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1]);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpTrilinear3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2]);
		    struct f_tdf_params p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp],
		       tau0_data,num_scatterings,tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_trilinear,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_trilinear,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_trilinear,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_trilinear,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_trilinear,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_trilinear,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv p=
		      {&energy_interpolator,beta,chempot[chemp],tau0_data,
		       num_scatterings,tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_trilinear_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_trilinear_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_trilinear_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_trilinear_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_trilinear_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_trilinear_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	    if (interpol_type == 1) {
	      IntpTricubic3<double> energy_interpolator(xBound,yBound,zBound,xMin,
							xSpacing,yMin,ySpacing,
							zMin,zSpacing,energy_data[band], true);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTricubic3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1], true);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpTricubic3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2], true);
		    struct f_tdf_params p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp],
		       tau0_data,num_scatterings,tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv p=
		      {&energy_interpolator,beta,chempot[chemp],tau0_data,
		       num_scatterings,tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	    if (interpol_type == 2) {
	      IntpTricubic3<double> energy_interpolator(xBound,yBound,zBound,xMin,
							xSpacing,yMin,ySpacing,
							zMin,zSpacing,energy_data[band], false);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTricubic3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1], false);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
                    IntpTricubic3<double>
                      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
                                                 xSpacing,yMin,ySpacing,
                                                 zMin,zSpacing,
                                                 velocity_data[band][dir2], false);
                    struct f_tdf_params p=
                      {&energy_interpolator,&velocity_interpolator_dir1,
                       &velocity_interpolator_dir2,beta,chempot[chemp],
                       tau0_data,num_scatterings,tau_energy_trans_data[band]};
                    // sigma
                    if (h) {
                      hcubature(fdim,f_tdf_0_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv p=
		      {&energy_interpolator,beta,chempot[chemp],tau0_data,
		       num_scatterings,tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	    if (interpol_type == 3) {
	      IntpAkimaUniform3<double> energy_interpolator(xBound,yBound,zBound,xMin,
							    xSpacing,yMin,ySpacing,
							    zMin,zSpacing,energy_data[band]);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpAkimaUniform3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1]);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpAkimaUniform3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2]);
		    struct f_tdf_params p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp],
		       tau0_data,num_scatterings,tau_energy_trans_data[band]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_akima,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_akima,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_akima,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_akima,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_akima,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_akima,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv p=
		      {&energy_interpolator,beta,chempot[chemp],tau0_data,
		       num_scatterings,tau_energy_trans_data[band],deriv1,deriv2};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_akima_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_akima_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_akima_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_akima_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		    // kappa
		    if (h) {
		    hcubature(fdim,f_tdf_2_akima_deriv,&p,dim,vmin,vmax,
			      max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_akima_deriv,&p,dim,vmin,vmax,
				max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=spin_degen[band]*value;
		  }
		}
	      }
	    }
	  }
          std::vector<std::vector<double> > sigma_inv(3,std::vector<double> (3));
          std::vector<std::vector<double> > seebeck_temp(3,std::vector<double> (3));
          std::vector<std::vector<double> > lorenz_highdeg(3,std::vector<double> (3));
          std::vector<std::vector<double> > lorenz_cor(3,std::vector<double> (3));
          std::vector<std::vector<double> > tmp(3,std::vector<double> (3));
          // invert sigma
	  invert_3x3_matrix(sigma[temp][chemp],sigma_inv);
          // now calculate the seebeck tensors and add units
          // also store the conductivity
          multiply_3x3_matrix(sigma_inv,chi[temp][chemp],seebeck_temp);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
	      cond[temp*num_chempot*9+9*chemp+3*dir1+dir2]=
		sigmaunit*sigma[temp][chemp][dir1][dir2]/temperatures[temp];
	      seebeck[temp*num_chempot*9+9*chemp+3*dir1+dir2]=-
		seebeckunit*seebeck_temp[dir1][dir2]/temperatures[temp];
            }
          }
          // now calculate the lorenz tensor and add units
          // first store the highly degenerate part (first term)
          multiply_3x3_matrix(kappa[temp][chemp],sigma_inv,lorenz_highdeg);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
	      lorenz[temp*num_chempot*9+9*chemp+3*dir1+dir2]=
		lorenzunit*lorenz_highdeg[dir1][dir2]/pow(temperatures[temp],2.0);
            }
          }
          // now calculate the correction term and subtract
          multiply_3x3_matrix(chi[temp][chemp],seebeck_temp,tmp);
          multiply_3x3_matrix(tmp,sigma_inv,lorenz_cor);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
              lorenz[temp*num_chempot*9+9*chemp+3*dir1+dir2]-=
                lorenzunit*lorenz_cor[dir1][dir2]/pow(temperatures[temp],2.0);
            }
          }
          // add symmetric values in the tensor 21=12 etc.
          cond[temp*num_chempot*9+9*chemp+3] = cond[temp*num_chempot*9+9*chemp+1];
          cond[temp*num_chempot*9+9*chemp+6] = cond[temp*num_chempot*9+9*chemp+2];
          cond[temp*num_chempot*9+9*chemp+7] = cond[temp*num_chempot*9+9*chemp+5];
          seebeck[temp*num_chempot*9+9*chemp+3] = seebeck[temp*num_chempot*9+9*chemp+1];
          seebeck[temp*num_chempot*9+9*chemp+6] = seebeck[temp*num_chempot*9+9*chemp+2];
          seebeck[temp*num_chempot*9+9*chemp+7] = seebeck[temp*num_chempot*9+9*chemp+5];
          lorenz[temp*num_chempot*9+9*chemp+3] = lorenz[temp*num_chempot*9+9*chemp+1];
          lorenz[temp*num_chempot*9+9*chemp+6] = lorenz[temp*num_chempot*9+9*chemp+2];
          lorenz[temp*num_chempot*9+9*chemp+7] = lorenz[temp*num_chempot*9+9*chemp+5];
	}
      }
      // clean up some mess
      for (int band=0;band<num_bands;band++) {
        delete[] tau_energy_trans_data[band];
      }
      delete[] tau_energy_trans_data;
    }
    // constant scattering
    else {
      if (info) {
        std::cout << ", constant scattering" << std::endl;
      }
      // fetch constant scattering time for each band
      // constant scattering always at slot 11
      std::vector<double> tau0(num_bands);
      for (int band=0;band<num_bands;band++) {
	// assume similar for all temperatures (first index of tau)
	tau0[band]=tau[num_scatterings*band+11];
      }
      // loop temperature
      for (int temp=0;temp<num_temps;temp++) {
	double beta = betafact/temperatures[temp];
	// loop eta
	for (int chemp=0;chemp<num_chempot;chemp++) {
          if (info) {
            std::cout << "Cubature/Wildmagic: constant scattering at " << temperatures[temp] << "K, chempot: " << chempot[chemp] << " eV" << std::endl;
          }
	  // loop bands
	  for (int band=0;band<num_bands;band++) {
	    if (interpol_type == 0) {
	      IntpTrilinear3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,xSpacing,
				    yMin,ySpacing,zMin,zSpacing,
				    energy_data[band]);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTrilinear3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1]);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpTrilinear3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2]);
		    struct f_tdf_params_const_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_trilinear_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_trilinear_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_trilinear_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_trilinear_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_trilinear_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_trilinear_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3] = {0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
                    // set deriv arrays
                    int deriv2[3]={0,0,0};
                    deriv2[dir2]=1;
                    struct f_tdf_params_deriv_const_scatter p=
                      {&energy_interpolator,beta,chempot[chemp],deriv1,deriv2};
                    // sigma
                    if (h) {
                      hcubature(fdim,f_tdf_0_trilinear_deriv_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_0_trilinear_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // chi
                    if (h) {
                      hcubature(fdim,f_tdf_1_trilinear_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_1_trilinear_deriv_const_scatter,&p,dim,
                                vmin,vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // kappa
                    if (h) {
                      hcubature(fdim,f_tdf_2_trilinear_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_2_trilinear_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
		      }
                    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                  }
		}
	      }
	    }
	    if (interpol_type == 1) {
	      IntpTricubic3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,xSpacing,
				    yMin,ySpacing,zMin,zSpacing,
				    energy_data[band], true);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTricubic3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1], true);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    IntpTricubic3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2], true);
		    struct f_tdf_params_const_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_const_scatter p=
		      {&energy_interpolator,beta,chempot[chemp],deriv1,deriv2};
                    // sigma
                    if (h) {
                      hcubature(fdim,f_tdf_0_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_0_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // chi
                    if (h) {
			hcubature(fdim,f_tdf_1_tricubic_deriv_const_scatter,&p,dim,vmin,
				  vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_1_tricubic_deriv_const_scatter,&p,dim,
				vmin,vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // kappa
                    if (h) {
                      hcubature(fdim,f_tdf_2_tricubic_deriv_const_scatter,&p,dim,vmin,
				  vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_2_tricubic_deriv_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                  }
                }
              }
            }
	    if (interpol_type == 2) {
	      IntpTricubic3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,xSpacing,
				    yMin,ySpacing,zMin,zSpacing,
				    energy_data[band], false);
	      // loop tensor
              // supply velocity grid
              if (!gen_velocities) {
                for (int dir1=0;dir1<3;dir1++) {
                  IntpTricubic3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1], false);
                  for (int dir2=dir1;dir2<3;dir2++) {
		    IntpTricubic3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2], false);
		    struct f_tdf_params_const_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_0_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_2_tricubic_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		  }
                }
              }
              else {
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
                    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_const_scatter p=
		      {&energy_interpolator,beta,chempot[chemp],deriv1,deriv2};
                    // sigma
                    if (h) {
                      hcubature(fdim,f_tdf_0_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_0_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // chi
                    if (h) {
                      hcubature(fdim,f_tdf_1_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_1_tricubic_deriv_const_scatter,&p,dim,
                                vmin,vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // kappa
                    if (h) {
                      hcubature(fdim,f_tdf_2_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_2_tricubic_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                  }
		}
	      }
	    }
	    if (interpol_type == 3) {
	      IntpAkimaUniform3<double>
		energy_interpolator(xBound,yBound,zBound,xMin,xSpacing,
				    yMin,ySpacing,zMin,zSpacing,
				    energy_data[band]);
              // supply velocity grid
              if (!gen_velocities) {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  IntpAkimaUniform3<double>
                    velocity_interpolator_dir1(xBound,yBound,zBound,xMin,
                                               xSpacing,yMin,ySpacing,
                                               zMin,zSpacing,
                                               velocity_data[band][dir1]);
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
                    IntpAkimaUniform3<double>
		      velocity_interpolator_dir2(xBound,yBound,zBound,xMin,
						 xSpacing,yMin,ySpacing,
						 zMin,zSpacing,
						 velocity_data[band][dir2]);
		    struct f_tdf_params_const_scatter p=
		      {&energy_interpolator,&velocity_interpolator_dir1,
		       &velocity_interpolator_dir2,beta,chempot[chemp]};
		    // sigma
		    if (h) {
		      hcubature(fdim,f_tdf_0_akima_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }

		    else {
		      pcubature(fdim,f_tdf_0_akima_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;

		    // chi
		    if (h) {
		      hcubature(fdim,f_tdf_1_akima_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    else {
		      pcubature(fdim,f_tdf_1_akima_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		    // kappa
		    if (h) {
		      hcubature(fdim,f_tdf_2_akima_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
                    else {
		      pcubature(fdim,f_tdf_2_akima_const_scatter,&p,dim,vmin,
				vmax,max_it,abs_err,rel_err,norm,&value,&error);
		    }
		    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
		  }
                }
              }
              else {
                // loop tensor
                for (int dir1=0;dir1<3;dir1++) {
                  // set deriv arrays
                  int deriv1[3]={0,0,0};
                  deriv1[dir1]=1;
                  if (iso) {
                    dir2_limit = dir1 + 1;
                  }
                  else {
                    dir2_limit = 3;
                  }
                  for (int dir2=dir1;dir2<dir2_limit;dir2++) {
		    // set deriv arrays
		    int deriv2[3]={0,0,0};
		    deriv2[dir2]=1;
		    struct f_tdf_params_deriv_const_scatter p=
		      {&energy_interpolator,beta,chempot[chemp],deriv1,deriv2};
                    // sigma
                    if (h) {
                      hcubature(fdim,f_tdf_0_akima_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_0_akima_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    sigma[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // chi
                    if (h) {
                      hcubature(fdim,f_tdf_1_akima_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_1_akima_deriv_const_scatter,&p,dim,
                                vmin,vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    chi[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                    // kappa
                    if (h) {
                      hcubature(fdim,f_tdf_2_akima_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    else {
                      pcubature(fdim,f_tdf_2_akima_deriv_const_scatter,&p,dim,vmin,
                                vmax,max_it,abs_err,rel_err,norm,&value,&error);
                    }
                    kappa[temp][chemp][dir1][dir2]+=tau0[band]*spin_degen[band]*value;
                  }
		}
	      }
	    }
	  }
          std::vector<std::vector<double> > sigma_inv(3,std::vector<double> (3));
          std::vector<std::vector<double> > seebeck_temp(3,std::vector<double> (3));
          std::vector<std::vector<double> > lorenz_highdeg(3,std::vector<double> (3));
          std::vector<std::vector<double> > lorenz_cor(3,std::vector<double> (3));
          std::vector<std::vector<double> > tmp(3,std::vector<double> (3));
          // invert sigma
	  invert_3x3_matrix(sigma[temp][chemp],sigma_inv);
          // now calculate the seebeck tensors and add units
          // also store the conductivity
          multiply_3x3_matrix(sigma_inv,chi[temp][chemp],seebeck_temp);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
	      cond[temp*num_chempot*9+9*chemp+3*dir1+dir2]=
		sigmaunit*sigma[temp][chemp][dir1][dir2]/temperatures[temp];
	      seebeck[temp*num_chempot*9+9*chemp+3*dir1+dir2]=-
		seebeckunit*seebeck_temp[dir1][dir2]/temperatures[temp];
            }
          }
          // now calculate the lorenz tensor and add units
          // first store the highly degenerate part (first term)
          multiply_3x3_matrix(kappa[temp][chemp],sigma_inv,lorenz_highdeg);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
	      lorenz[temp*num_chempot*9+9*chemp+3*dir1+dir2]=
		lorenzunit*lorenz_highdeg[dir1][dir2]/pow(temperatures[temp],2.0);
            }
          }
          // now calculate the correction term and subtract
          multiply_3x3_matrix(chi[temp][chemp],seebeck_temp,tmp);
          multiply_3x3_matrix(tmp,sigma_inv,lorenz_cor);
          for (int dir1=0;dir1<3;dir1++) {
            for (int dir2=0;dir2<3;dir2++) {
              lorenz[temp*num_chempot*9+9*chemp+3*dir1+dir2]-=
                lorenzunit*lorenz_cor[dir1][dir2]/pow(temperatures[temp],2.0);
            }
          }
          // add symmetric values in the tensor 21=12 etc.
          cond[temp*num_chempot*9+9*chemp+3] = cond[temp*num_chempot*9+9*chemp+1];
          cond[temp*num_chempot*9+9*chemp+6] = cond[temp*num_chempot*9+9*chemp+2];
          cond[temp*num_chempot*9+9*chemp+7] = cond[temp*num_chempot*9+9*chemp+5];
          seebeck[temp*num_chempot*9+9*chemp+3] = seebeck[temp*num_chempot*9+9*chemp+1];
          seebeck[temp*num_chempot*9+9*chemp+6] = seebeck[temp*num_chempot*9+9*chemp+2];
          seebeck[temp*num_chempot*9+9*chemp+7] = seebeck[temp*num_chempot*9+9*chemp+5];
          lorenz[temp*num_chempot*9+9*chemp+3] = lorenz[temp*num_chempot*9+9*chemp+1];
          lorenz[temp*num_chempot*9+9*chemp+6] = lorenz[temp*num_chempot*9+9*chemp+2];
          lorenz[temp*num_chempot*9+9*chemp+7] = lorenz[temp*num_chempot*9+9*chemp+5];
	}
      }
    }
    // clean remaining mess
    for (int band=0;band<num_bands;band++) {
      for (int dir=0;dir<3;dir++) {
	for (int k=0;k<zBound;k++) {
	  for (int j=0;j<yBound;j++) {
	    if (dir==0) {
	      delete[] energy_data[band][k][j];
	    }
	    delete[] velocity_data[band][dir][k][j];
	  }
	  if (dir==0) {
	    delete[] energy_data[band][k];
	  }
	  delete[] velocity_data[band][dir][k];
	}
	if (dir==0) {
	  delete[] energy_data[band];
	}
	delete[] velocity_data[band][dir];
      }
      delete[] velocity_data[band];
    }
    delete[] energy_data;
    delete[] velocity_data;
    }
  else {
    std::cout << "This interpol type is not supported. Exiting." << std::endl;
    exit(1);
  }
}

void calc_total_dos(double *dos_decomp, int num_decomps,
		    int num_energy_samples,double *dos) {

  for (int energy=0;energy<num_energy_samples;energy++) {
    dos[energy]=0.0;
    for (int decomp=0;decomp<num_decomps;decomp++) {
      dos[energy]+=dos_decomp[decomp*num_energy_samples+energy];
    }
  }

}

void calc_dos_cubature(int *num_points,
		       double *domainx, double *domainy,
		       double *domainz, double *energies,
		       double *energy_samples, int *spin_degen,
		       int num_energy_samples, int num_bands,
		       double sigma, int interpol_method,
		       int smear_then_interpolate,
                       double volume, double volume_bz,
		       unsigned max_it,
		       double abs_err, double rel_err, int integral_type,
                       int info,
		       double *dos, double *int_dos) {
  bool h = false;
  // set type of cubature
  if (integral_type == 1) {
    h = true;
  }

  // initialize domain
  int xBound=num_points[0];
  double xMin=domainx[0];
  double xSpacing=(domainx[1]-domainx[0])/(xBound-1);
  double xMax=xMin+xSpacing*(xBound-1);
  int yBound=num_points[1];
  double yMin=domainy[0];
  double ySpacing=(domainy[1]-domainy[0])/(yBound-1);
  double yMax=yMin+ySpacing*(yBound-1);
  int zBound=num_points[2];
  double zMin=domainz[0];
  double zSpacing=(domainz[1]-domainz[0])/(zBound-1);
  double zMax=zMin+zSpacing*(zBound-1);
  int tot_num_points=xBound*yBound*zBound;

  double ****energies_temp=new double ***[num_bands];
  double ***delta_approx_energy_slit=new double **[zBound];
  for (int band=0;band<num_bands;band++) {
    energies_temp[band]=new double **[zBound];
    for (int k=0;k<zBound;k++) {
      energies_temp[band][k]=new double*[yBound];
      if (band==0) {
	delta_approx_energy_slit[k]=new double*[yBound];
      }
      for (int j=0;j<yBound;j++) {
	energies_temp[band][k][j]=new double[xBound];
	if (band==0) {
	  delta_approx_energy_slit[k][j]=new double[xBound];
	}
	for (int i=0;i<xBound;i++) {
	  energies_temp[band][k][j][i]=
	    energies[band*tot_num_points+
		     i*zBound*yBound+j*zBound+k];
	}
      }
    }
  }
  double vmin[3]={xMin,yMin,zMin};
  double vmax[3]={xMax,yMax,zMax};
  unsigned fdim=1;
  unsigned dim=3;
  error_norm norm=ERROR_INDIVIDUAL;
  double value=0.0;
  double error=0.0;

  // loop bands
  for (int band=0;band<num_bands;band++) {
    double idos = 0.0;
    for (int energy=0;energy<num_energy_samples;energy++) {
      if (info) {
        std::cout << "Cubature/Wildmagic: density of states at " << energy_samples[energy] << " eV for band: " << band << std::endl;
      }
      // now start the integration
      // smear first, then interpolate
      if (smear_then_interpolate) {
	gaussian(xBound,yBound,zBound,energies_temp[band],
		 delta_approx_energy_slit,energy_samples[energy],sigma);
	if (interpol_method==0) {
	  IntpTrilinear3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,
			 zMin,zSpacing,delta_approx_energy_slit);
	  if (h) {
	    hcubature(fdim,f_trilinear,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_trilinear,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
	else if (interpol_method==1) {
	  IntpTricubic3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,
			 zMin,zSpacing,delta_approx_energy_slit,true);
	  if (h) {
	    hcubature(fdim,f_tricubic,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_tricubic,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
	else if (interpol_method==2) {
	  IntpTricubic3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,
			 zMin,zSpacing,delta_approx_energy_slit,false);
	  if (h) {
	    hcubature(fdim,f_tricubic,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_tricubic,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
	else {
	  IntpAkimaUniform3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,
			 zMin,zSpacing,delta_approx_energy_slit);
	  if (h) {
	    hcubature(fdim,f_akima,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_akima,&interpolator,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
      }
      else {
	if (interpol_method==0) {
	  IntpTrilinear3<double>
	    interpolator(xBound,yBound,zBound,xMin,xSpacing,yMin,
			 ySpacing,zMin,zSpacing,energies_temp[band]);
	  struct f_gaussian_params p=
	    {&interpolator,sigma,energy_samples[energy]};
	  if (h) {
	    hcubature(fdim,f_trilinear_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_trilinear_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
	else if (interpol_method==1) {
	  IntpTricubic3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,zMin,
			 zSpacing,energies_temp[band],true);
	  struct f_gaussian_params p=
	    {&interpolator,sigma,energy_samples[energy]};
	  if (h) {
	    hcubature(fdim,f_tricubic_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_tricubic_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
	else if (interpol_method==2) {
	  IntpTricubic3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,
			 zMin,zSpacing,energies_temp[band],false);
	  struct f_gaussian_params p=
	    {&interpolator,sigma,energy_samples[energy]};
	  if (h) {
	    hcubature(fdim,f_tricubic_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_tricubic_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
	else {
	  IntpAkimaUniform3<double>
	    interpolator(xBound,yBound,zBound,xMin,
			 xSpacing,yMin,ySpacing,
			 zMin,zSpacing,energies_temp[band]);
	  struct f_gaussian_params p=
	    {&interpolator,sigma,energy_samples[energy]};
	  if (h) {
	    hcubature(fdim,f_akima_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	  else {
	    pcubature(fdim,f_akima_gaussian,&p,dim,vmin,vmax,
		      max_it,abs_err,rel_err,norm,&value,&error);
	  }
	}
      }
      // calculate regular dos
      int index=band*num_energy_samples+energy;
      dos[index]=spin_degen[band]*value/volume_bz/volume;
      // calculate integrated dos
      idos+=dos[index];
      int_dos[index]=idos;
    }
  }
  // clean all the mess
  for (int band=0;band<num_bands;band++) {
    for (int k=0;k<zBound;k++) {
      for (int j=0;j<yBound;j++) {
	delete[] energies_temp[band][k][j];
	if (band==0) {
	  delete[] delta_approx_energy_slit[k][j];
	}
      }
      delete[] energies_temp[band][k];
      if (band==0) {
	delete[] delta_approx_energy_slit[k];
      }
    }
    delete[] energies_temp[band];
    if (band==0) {
      delete[] delta_approx_energy_slit;
    }
  }
  delete[] energies_temp;
  }


void multiply_3x3_matrix(std::vector<std::vector<double> >
                         &input_matrix_a,
                         std::vector<std::vector<double> >
                         &input_matrix_b,
                         std::vector<std::vector<double> >
                         &output_matrix) {

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      output_matrix[i][j]=0.0;
      for(int k=0;k<3;k++){
        output_matrix[i][j]+=input_matrix_a[i][k]*input_matrix_b[k][j];
      }
    }
  }
}

void invert_3x3_matrix(std::vector<std::vector<double> >
                       &input_matrix,
                       std::vector<std::vector<double> >
                       &output_matrix) {

  double detinv=1.0/determinant_3x3_matrix(input_matrix);
  adjoint_and_scale_3x3_matrix(input_matrix,detinv,output_matrix);

}

void adjoint_and_scale_3x3_matrix(std::vector<std::vector<double> >
                                  &input_matrix,
                                  double scale,
                                  std::vector<std::vector<double> >
                                  &output_matrix) {

  output_matrix[0][0] = (scale) * (input_matrix[1][1] *
                                   input_matrix[2][2] -
                                   input_matrix[1][2] *
                                   input_matrix[2][1]);
  output_matrix[1][0] = (scale) * (input_matrix[1][2] *
                                   input_matrix[2][0] -
                                   input_matrix[1][0] *
                                   input_matrix[2][2]);
  output_matrix[2][0] = (scale) * (input_matrix[1][0] *
                                   input_matrix[2][1] -
                                   input_matrix[1][1] *
                                   input_matrix[2][0]);
  output_matrix[0][1] = (scale) * (input_matrix[0][2] *
                                   input_matrix[2][1] -
                                   input_matrix[0][1] *
                                   input_matrix[2][2]);
  output_matrix[1][1] = (scale) * (input_matrix[0][0] *
                                   input_matrix[2][2] -
                                   input_matrix[0][2] *
                                   input_matrix[2][0]);
  output_matrix[2][1] = (scale) * (input_matrix[0][1] *
                                   input_matrix[2][0] -
                                   input_matrix[0][0] *
                                   input_matrix[2][1]);
  output_matrix[0][2] = (scale) * (input_matrix[0][1] *
                                   input_matrix[1][2] -
                                   input_matrix[0][2] *
                                   input_matrix[1][1]);
  output_matrix[1][2] = (scale) * (input_matrix[0][2] *
                                   input_matrix[1][0] -
                                   input_matrix[0][0] *
                                   input_matrix[1][2]);
  output_matrix[2][2] = (scale) * (input_matrix[0][0] *
                                   input_matrix[1][1] -
                                   input_matrix[0][1] *
                                   input_matrix[1][0]);
}

double determinant_3x3_matrix(std::vector<std::vector<double> >
                              &matrix) {
  double determinant =
    matrix[0][0] * (matrix[1][1]*matrix[2][2] -
                    matrix[1][2] * matrix[2][1]);
  determinant -=
    matrix[0][1] * (matrix[1][0]*matrix[2][2] -
                    matrix[1][2] * matrix[2][0]);
  determinant +=
    matrix[0][2] * (matrix[1][0]*matrix[2][1] -
                    matrix[1][1] * matrix[2][0]);
  return determinant;
}

double maximumValue(double *array,int size) {
  double max = array[0];
  for(int i=1; i<size; i++) {
    if(array[i] > max) {
      max = array[i];
    }
  }
  return max;
}
double minimumValue(double *array,int size) {
  double min = array[0];
  for(int i=1; i<size; i++) {
    if(array[i] < min) {
      min = array[i];
    }
  }
  return min;
}
