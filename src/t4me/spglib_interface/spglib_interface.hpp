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
#include <numeric>
#include <cmath>
#include <vector>
extern "C" {
#include <spglib.h>
#include <tetrahedron_method.h>
}

namespace constants {
  const double KB=8.6173324;
  const double PI=3.1415926535897932;
  const double G=7.7480917310;
  const double HBAR=6.582119514;
}

void get_reciprocal_mesh_interface(int *mesh, double *lattice, double
				   *positions, int *anumbers,
				   unsigned int num_atoms,
				   int *is_shifted, int *mesh_points,
				   int *mapping,
                                   char *intsym,
				   int is_time_reversal,
				   double symprec);

void calc_dos_weights(double *energies_temp, int *spg_grid_address_temp,
                      int *mapping_bz_to_ibz_temp,
                      int *mapping_ibz_to_bz_temp,
                      int *ibz_weights_temp,
                      int *mesh_temp,
                      double *rec_basis_temp, double *energy_samples,
                      int num_energy_samples,int num_bands,
                      int num_kpoints_ibz,
                      int *spin_degen_temp,
                      double volume, double volume_bz,
                      int weight_type,
                      double smearing,
                      double *dos, double *int_dos);

void calc_transport_tensors_weights(double *energies_temp, double *velocities_temp,
				  double *scattering_temp,
				  double *temperatures_temp,
				  double *chempots_temp,
				  int *spg_grid_address_temp,
				  int *mapping_bz_to_ibz_temp,
				  int *mapping_ibz_to_bz_temp,
                                  int *ibz_weights_temp,
				  int *mesh_temp, double *rec_basis_temp,
				  int num_bands, int num_kpoints_ibz,
				  int num_temp_steps,
				  int num_chempot_steps,
				  int int_method,
                                  int energy_samples, int weight_type,
                                  double smearing,
				  double chempot_cutoff, double volume,
				  int *spin_degen_temp,
				  double *cond,
				  double *seebeck,
				  double *lorenz);

void e_integral(std::vector<std::vector<std::vector<double> > >
		&sigma_k2,
		std::vector<std::vector<std::vector<std::vector<double>
		> > > &sigma_k3,
		std::vector<std::vector<std::vector<double> > > &scattering,
                double *energy_samples, double temperature,
		double chempot, int int_method,
		std::vector<std::vector<double> > &sigma_e2,
		std::vector<std::vector<std::vector<double> > >
		&sigma_e3,
		std::vector<std::vector<double> > &chi_e,
		std::vector<std::vector<double> > &kappa_e);

double e_integrand(double f, double energy, double chempot, double
		   beta,int i);

double dedf(double energy, double chempot, double beta);

void spectral_sum(std::vector<std::vector<std::vector<double> > >
		  &weights,
		  std::vector<std::vector<std::vector<double>
		  > > &int_weights,
		  std::vector<std::vector<std::vector<double> > >
		  &velocities,
		  std::vector<std::vector<std::vector<double> > > &scattering,
                  std::vector<int> &mapping_bz_to_ibz,
		  std::vector<std::vector<std::vector<std::vector<double> > > >
		  &sigma2,
		  std::vector<std::vector<std::vector<std::vector<std::vector<double>
		  > > > > &sigma3);

void smearing_weights(std::vector<std::vector<double> > &energies,
		      std::vector<int> &ibz_weights, double *energy_samples,
		      int num_energy_samples, int num_kpoints,
		      double sigma, bool ibz,
		      std::vector<int> &spin_degen,
		      std::vector<std::vector<std::vector<double> > >
		      &weights,
		      std::vector<std::vector<std::vector<double> > >
		      &int_weights);

void tetra_weights(std::vector<std::vector<double> > &energies,
		   std::vector<std::vector<int> > &spg_grid_address,
                   std::vector<int> &spg_grid_cantor_address,
                   std::vector<int> &mapping_bz_to_ibz,
                   std::vector<int> &mapping_ibz_to_bz,
                   std::vector<int> &ibz_weights,
                   std::vector<int> &mesh,
		   double rec_basis[3][3], double *energy_samples,
		   int num_bands, int num_kpoints_ibz,
		   bool ibz, std::vector<int> &spin_degen,
		   std::vector<std::vector<std::vector<double> > >
		   &kpoint_weights,
		   std::vector<std::vector<std::vector<double> > >
		   &int_weights);

static int grid_address_to_index(std::vector<int> &v, std::vector<int> &mesh);

void back_to_cell(int &gp, int &sampling);

int cantor_transform_3_i(std::vector<int> &v);

double gaussian(double energy, double energy_ref, double sigma);

void fetch_energy_samples(std::vector<double> &chempots,int
			  num_chempot_samples,double energy_cutoff,
			  double *energy_samples, int
			  num_energy_samples);

void invert_3x3_matrix(std::vector<std::vector<double> >
		       &input_matrix,
		       std::vector<std::vector<double> >
		       &output_matrix);

void multiply_3x3_matrix(std::vector<std::vector<double> >
                         &input_matrix_a,
                         std::vector<std::vector<double> >
                         &input_matrix_b,
                         std::vector<std::vector<double> >
                         &output_matrix);

void adjoint_and_scale_3x3_matrix(std::vector<std::vector<double> >
				  &input_matrix,
				  double scale,
				  std::vector<std::vector<double> >
				  &output_matrix);

double determinant_3x3_matrix(std::vector<std::vector<double> >
			      &matrix);
