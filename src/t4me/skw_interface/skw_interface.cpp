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

#include "skw_interface.hpp"

int interpolate_skw_interface(double *energies, int num_bands, double *kpoints, int num_kpoints, int *ksampling, double *lattice, double *positions, int *species, int num_atoms, double star_radius_factor, double *ienergies, double *ivelocities, double *ikpoints, int num_ikpoints, int *iksampling) {

  // convert everything to STL vectors
  // first the energies
  std::vector<std::vector<double> > energies_temp(num_bands, std::vector<double> (num_kpoints));
  for (int band=0;band<num_bands;band++) {
    for (int kpoint=0;kpoint<num_kpoints;kpoint++) {
      energies_temp[band][kpoint]=*((double *)energies + band*num_kpoints + kpoint);
    }
  }
  // then the kpoints
  std::vector<std::vector<double> > kpoints_temp(num_kpoints, std::vector<double> (3));
  for (int kpoint=0;kpoint<num_kpoints;kpoint++) {
    for (int dir=0;dir<3;dir++) {
      kpoints_temp[kpoint][dir]=*((double *)kpoints + kpoint*3+dir);
    }
  }
  std::vector<int> ksampling_temp(3);
  for (int dir=0;dir<3;dir++) {
    ksampling_temp[dir]=ksampling[dir];
  }
  // then the lattice
  std::vector<std::vector<double> > lattice_temp(3, std::vector<double> (3));
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      lattice_temp[i][j]=*((double *)lattice + i*3 + j);
    }
  }
  // then positions (really stupid, np->vector->array, but does not cause much
  // overhead)
  std::vector<std::vector<double> > positions_temp(num_atoms, std::vector<double> (3));
  std::vector<int> species_temp(num_atoms);
  for (int atom=0;atom<num_atoms;atom++) {
    species_temp[atom]=*((int *)species+atom);
    for (int dir=0;dir<3;dir++) {
      positions_temp[atom][dir]=*((double *)positions + atom*3 + dir);
    }
  }

  // now fetch the point group operations
  // sym_ops holds the point group operations returned by spglib
  std::vector<std::vector<std::vector<int> > > sym_ops;
  fetch_sym_ops(lattice_temp, positions_temp, species_temp, sym_ops);

  // generate lattice vectors r
  // holds the minimum set of r (maximum symmetry)
  std::vector<std::vector<int> > r;
  // holds the unique r folded out by the point group operations
  std::vector<std::vector<std::vector<int> > > r_sym_span;
  // holds the lenght of the r vectors (for the r vector)
  std::vector<double> rl;
  generate_lat_vec(lattice_temp,sym_ops,star_radius_factor,ksampling_temp,r,r_sym_span,rl);

  // generate the fourier components, or star functions, S
  std::vector<std::vector<double> > s;
  generate_s(r,kpoints_temp,sym_ops,s);

  // generate the roughness function rho
  std::vector<double> c;
  std::vector<double> rho;
  generate_rho(r,rl,c,rho);

  // generate the H matrix, also stor S times S divded by rho (ssr)
  std::vector<double> h;
  std::vector<std::vector<double> > ssr;
  generate_h_and_ssr(s,rho,h,ssr);

  // generate the energy difference, i.e. the input to the linear
  // equation solver below (stored in the lambda array)
  std::vector<double> lambda;
  generate_energies_diff(energies_temp,lambda);

  // generate the lambda solutions
  int info=generate_lambdas(energies_temp,h,lambda);
  if (info) {
    // serious error during the generation of lambdas (LAPACK/BLAS)
    return 1;
  }

  // generate the epsilon coefficients
  std::vector<std::vector<double> > epsilons;
  generate_epsilons(energies_temp,s,ssr,lambda,epsilons);

  // check that the interpolated energies return the orginial datapoint,
  // within TOL defined in skw.h
  check_interpolated_energies(energies_temp,s,epsilons);

  // everything is set up, good to go, interpolate energies and velocities
  // and fetch the kpoint set
  std::vector<std::vector<double> > ienergies_temp;
  std::vector<std::vector<std::vector<double> > > ivelocities_temp;
  std::vector<std::vector<double> > ikpoints_temp;
  std::vector<int> iksampling_temp(3);
  interpolate(epsilons,r_sym_span,lattice_temp,ienergies_temp,ivelocities_temp,ikpoints_temp,iksampling_temp);

  // copy content from the vector STL in order to pass data back to python
  // first the energies
  num_ikpoints=ikpoints_temp.size();
  for (int band=0;band<num_bands;band++) {
    for (int kpoint=0;kpoint<num_ikpoints;kpoint++) {
      ienergies[band*num_ikpoints+kpoint]=ienergies_temp[band][kpoint];
    }
  }
  // then the velocities
  for (int band=0;band<num_bands;band++) {
    for (int dir=0; dir < 3; dir++) {
      for (int kpoint=0;kpoint<num_ikpoints;kpoint++) {
	ivelocities[band*3*num_ikpoints + dir*num_ikpoints+kpoint]=ivelocities_temp[band][dir][kpoint];
      }
    }
  }
  // then the kpoints
  for (int kpoint=0;kpoint<num_ikpoints;kpoint++) {
    for (int dir=0;dir<3;dir++) {
      ikpoints[kpoint*3 + dir] = ikpoints_temp[kpoint][dir];
    }
  }
  // then the iksampling
  for (int dir=0;dir<3;dir++) {
    iksampling[dir]=iksampling_temp[dir];
  }

  return 0;

}
