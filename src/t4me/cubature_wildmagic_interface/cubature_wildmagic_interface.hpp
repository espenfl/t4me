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
#include <WildMagic5/Wm5IntpAkimaUniform1.h>
#include <WildMagic5/Wm5Math.h>
#include <cubature.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <chrono>
#include <limits>

class Timer {
public:
  Timer() : beg_(clock_::now()) {}
  void reset() { beg_ = clock_::now(); }
  double elapsed() const {
    return std::chrono::duration_cast<second_>
      (clock_::now() - beg_).count(); }

private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1> > second_;
  std::chrono::time_point<clock_> beg_;
};


struct f_gaussian_params {
  void *interpolator;
  double sigma;
  double energy_ref;
};

struct f_tdf_params {
  void *energy_interpolator;
  void *velocity_interpolator_l;
  void *velocity_interpolator_m;
  double beta;
  double chempot;
  double *tau0;
  int num_scatterings;
  double *tau_energy_trans;
};

struct f_tdf_params_num_scatter {
  void *energy_interpolator;
  void *velocity_interpolator_l;
  void *velocity_interpolator_m;
  void *tau_interpolator;
  double beta;
  double chempot;
  double energy_trans;
};

struct f_tdf_params_const_scatter {
  void *energy_interpolator;
  void *velocity_interpolator_l;
  void *velocity_interpolator_m;
  double beta;
  double chempot;
};

struct f_tdf_params_deriv {
  void *energy_interpolator;
  double beta;
  double chempot;
  double *tau0;
  int num_scatterings;
  double *tau_energy_trans;
  int *derivl;
  int *derivm;
};

struct f_tdf_params_deriv_num_scatter {
  void *energy_interpolator;
  void *tau_interpolator;
  double beta;
  double chempot;
  double energy_trans;
  int *derivl;
  int *derivm;
};

struct f_tdf_params_deriv_const_scatter {
  void *energy_interpolator;
  double beta;
  double chempot;
  int *derivl;
  int *derivm;
};

void cubature_wildmagic_execute_integration(int *num_points,
					    double *domainx,
					    double *domainy,
					    double *domainz,
					    double *data,
					    double &value,
					    double &error,
					    int itype);

void calc_transport_tensors_cubature(int *num_points,
				     double *domainx,
				     double *domainy,
				     double *domainz,
				     double *energy, int num_bands,
				     double *velocity, double *tau,
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
				     int scattering_type);

void calc_total_dos(double *dos_decomp, int num_decomps,
		    int num_energy_samples,double *dos);

void calc_dos_cubature(int *num_points,
		       double *domainx, double *domainy,
		       double *domainz,
		       double *energies, double *energy_samples,
		       int *spin_degen,
		       int num_energy_samples, int num_bands,
		       double sigma, int interpol_type,
		       int smear_then_interpolate,
                       double volume, double volume_bz,
		       unsigned max_it, double abs_err,
		       double rel_err, int integral_type,
                       int info,
		       double *dos, double *int_dos);

int hcubature(unsigned, integrand, void *,unsigned, const double *,
	      const double *, size_t, double, double, error_norm,
	      double *, double *);

int f_tdf_0_trilinear(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_1_trilinear(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_2_trilinear(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_0_trilinear_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_1_trilinear_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_2_trilinear_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_0_trilinear_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_1_trilinear_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_2_trilinear_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_0_trilinear_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_1_trilinear_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_2_trilinear_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_0_trilinear_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_1_trilinear_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_2_trilinear_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_0_tricubic(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_1_tricubic(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_2_tricubic(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_0_tricubic_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_1_tricubic_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_2_tricubic_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_0_tricubic_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_1_tricubic_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_2_tricubic_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_0_tricubic_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_1_tricubic_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_2_tricubic_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_0_tricubic_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_1_tricubic_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_2_tricubic_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_0_akima(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_1_akima(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_2_akima(unsigned,const double *, void *,
		  unsigned, double *);

int f_tdf_0_akima_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_1_akima_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_2_akima_num_scatter(unsigned,const double *, void *,
			      unsigned, double *);

int f_tdf_0_akima_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_1_akima_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_2_akima_deriv(unsigned,const double *, void *,
			unsigned, double *);

int f_tdf_0_akima_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_1_akima_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_2_akima_deriv_num_scatter(unsigned,const double *, void *,
				    unsigned, double *);

int f_tdf_0_akima_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_1_akima_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_tdf_2_akima_deriv_const_scatter(unsigned,const double *,
				      void *, unsigned, double *);

int f_akima(unsigned,const double *, void *, unsigned, double *);

int f_akima_gaussian(unsigned,const double *, void *,
		     unsigned, double *);

int f_tricubic(unsigned,const double *, void *,
	       unsigned, double *);

int f_tricubic_gaussian(unsigned,const double *, void *,
			unsigned, double *);

int f_trilinear(unsigned,const double *, void *,
		unsigned, double *);

int f_trilinear_gaussian(unsigned,const double *, void *,
			 unsigned, double *);

int gaussian(int, int, int, double *, double *, double, double);

double fetch_combined_scattering_parabolic(double, double *,
					   double *);

double generate_scattering_array(double *, int, double *, double **,
				 int, double *, double *);

void multiply_3x3_matrix(std::vector<std::vector<double> >
                         &input_matrix_a,
                         std::vector<std::vector<double> >
                         &input_matrix_b,
                         std::vector<std::vector<double> >
                         &output_matrix);


double determinant_3x3_matrix(std::vector<std::vector<double> >
                              &matrix);

void adjoint_and_scale_3x3_matrix(std::vector<std::vector<double> >
                                  &input_matrix,
                                  double scale,
                                  std::vector<std::vector<double> >
                                  &output_matrix);

void invert_3x3_matrix(std::vector<std::vector<double> >
                       &input_matrix,
                       std::vector<std::vector<double> >
                       &output_matrix);

double maximumValue(double *,int);

double minimumValue(double *,int);
