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

#include <math.h>
#include "skw.h"

/* Dummy main, as we only use this in library mode */

int main(int argc, const char* argv[]) {

}

/* Fetches the symmetry operations using spglib

   Since spglib has a static C interface, we need to copy the 
   lattice and positions array from STL vectors to C arrays before
   calling the spglib routines
  
 */

void fetch_sym_ops(std::vector<std::vector<double> > &lattice, std::vector<std::vector<double> > &positions, std::vector<int> &species, std::vector<std::vector<std::vector<int> > > &sym_ops) {
  // declare some static arrays for spglib
  int num_atoms=positions.size();
  int max_size=192;
  int rotation[max_size][3][3];
  double translation[max_size][3];
  double lattice_temp[3][3];
  char symbol[11];
  double positions_temp[num_atoms][3];
  // copy lattice and positions 
  // (transpose the lattice) due to the C character of spglib
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      lattice_temp[i][j]=lattice[j][i];
    }
  }
  // really stupid copy process, but since we insist of using vectors, this is how it
  // is going to be
  for (int atom=0; atom<num_atoms;atom++) {
    for (int dir=0;dir<3;dir++) {
      positions_temp[atom][dir]=positions[atom][dir];
    }
  }
  // call spglib to obtain the sym_ops
  int num_ops=spg_get_symmetry(rotation, translation, max_size, lattice_temp, positions_temp,&species[0],num_atoms,constants::SYMPREC);
  // call spglib to get the space group etc.
  int space_group=spg_get_international(symbol, lattice_temp, positions_temp, &species[0], num_atoms, constants::SYMPREC);
  // copy the symmetry to the sym_ops vector
  sym_ops.resize(num_ops,std::vector<std::vector<int> > (3, std::vector<int> (3)));
  for (int ops=0;ops<num_ops;ops++) {
    DEBUG_DUMP_SYMOPS(std::cout << "sym_ops #: " << ops << std::endl);
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
	sym_ops[ops][i][j]=rotation[ops][i][j];
      }
      DEBUG_DUMP_SYMOPS(std::cout << sym_ops[ops][i][0] << ' ' << sym_ops[ops][i][1] << ' ' << sym_ops[ops][i][2] << std::endl);
    }
  }
  DEBUG_DUMP_SYMOPS(std::cout << "space group: " << space_group << " found" << std::endl);
}

/*
   Generate the lattice vectors R.

   Due to the symmetry of a crystal it is possible to include integer multiples of R
   and use this, in principle infinite overdetermination to determine e.g. the coefficients
   of a set of linear equations.

   In this context the lattice vectors are used to interpolate the band dispersion by using
   EFL:   a plane wave expansion of these lattice vectors and (NEW AND ORGINAL?) the kpoints. 
   This interpolation scheme is often termed Fourier interpolation, star interpolation etc.
   We use the translation symmetry to our advantage. In addition the point symmetry is used
   to exclude symmetry equivalent generated lattice vectors.

   REF: Pickett et al. PRB 38, 2721 (1988) and similar references
   
   Assume the plane wave expansion (or "star function" as it is frequently called)

   S_m(\vec{k})=\sum_J exp(i(J\vec{R}_m)\cdot\vec{k})/m,

   for n point group operations J (matrix) on the lattice vector

   \vec{R}_m=h\vec{a}+k\vec{b}+l\vec{c},

   where h, k and l is a combination of integers (for a given m) and \vec{a}, \vec{b} and 
   \vec{c} is the lattice vectors of the unit cell (in this case).

   By using the plane wave expansion we can expand the dispersion relations as

   \epsilon(\vec{k})=\sum_{m=1}^M \epsilon_m (\vec{k}) S_m(\vec{k}).

   For M data points. This relation needs to be minimized upon some constrictions.

   \epsilon(\vec{k})=\epsilon'(\vec{k}), where \epsilon' is the input dispersion

   This relation needs to be minimized. A natural choice is to minimize a roughness. This is
   introduced though the roughness function (after the reference above, other choices are possible)
   
   R_o=(1/N)\sum_k(C_0\epsilon(\vec{k})+C_1|\nabla \epsilon(\vec{k})|^2+C_2|\nabla^2 \epsilon(\vec{k})|^2+...)

   R_o=\sum_{m=1}^M |\epsilon_m|^2 \rho(R_m),

   where N is the number of k points in the BZ. R_m is the length of the lattice vector \vec{R}_m.
   C_i are some user adjustable constants to be determined.
   
   From the reference above, better fit is obtained by dropping the m=1 terms.

   Minimizing the roughness function R_o by a Lagrange multiplier \lambda yields the following
   sets of linear equations

   1. \epsilon_m*\rho_m=\sum_i^N\lambda_iS_m(\vec{k}_i) for m>1 and
   
   2. 0=\sum_{i=1}^N\lambda_i together with the constraint of

   3. \epsilon(\vec{k})=\epsilon'(\vec{k}) (more a demand), where \epsilon is the input dispersion

   Following the reference above one can choose e.g. the final k point as a reference (this can be
   any of the input k points) we get

   \delta \epsilon_j=\epsilon(\vec{k}_j)-\epsilon(\vec{k}_N)=\sum_{i=1}^{N-1} H_{ji}\lambda_i*,
   
   where

   H_{ji}=\sum_{m=2}^M (S_m(\vec{k}_j)-S_m(\vec{k}_N))(S_m*(\vec{k}_i)-S_m*(\vec{k}_N))/\rho_m

   Before solving these equations, the R_m grid need to be determined. The choice of R_m is fully up to the user. 
   Here we force the point group operations on R_m in order to avoid symmetry equivalent points. We span the hkl
   integers as a ellipsoid determined by the length of the lattice vectors.
   EFL: THIS DOES NOT WORK FOR E.G. VERY DISTORTED STRUCTURES. REVERT TO SPHERE?

   When the R_m grid is determined, S can be genrated, H can be laid out and we can solve the
   linear equations to obtain the \lambda_i above. When this is done, we have the following 
   equations to determined the energy at the m points.

   \epsilon_m=\sum_{i=1}^{N-1}\lambda_i*(S_m*(\vec{k}_i)-S_m*(\vec{k}_N))/\rho_m, for m>1

   and
   
   \epsilon_1=\epsilon(\vec{k}_N)-\sum_{m=2}^{M}\epsilon_m S_m(\vec{k}_N)

   The symmetry tag (set to 0) turn of exclusion of symmetry points (e.g. if we want to run
   on the full BZ)

*/

void generate_lat_vec(std::vector<std::vector<double> > &lattice, std::vector<std::vector<std::vector<int> > > &sym_ops, double &star_radius_factor, std::vector<int> &ksampling, std::vector<std::vector<int> > &r, std::vector<std::vector<std::vector<int> > > &r_sym_span,std::vector<double> &rl) {
  int rmax[3];
  std::vector<double> temp_vector(3,0.0), ttemp_vector(3,0.0);
  std::vector<std::vector<int> > r_tmp;
  std::vector<int> r_index;
  std::vector<double> r_length;
  std::vector<std::vector<int> > sym_vectors;
  std::vector<std::vector<double> > lattice_transposed(3, std::vector<double> (3));
  std::vector<double> lattice_vector_length(3);
  int num_kpoints=1;
  for (int dir=0; dir<3; dir++) {
    num_kpoints*=ksampling[dir];
  }

  // calculate length of the lattice vectors
  for (int dir=0;dir<3;dir++) {
    lattice_vector_length[dir]=length_of_vector_d(lattice[dir]);
  }  
  // locate the longest lattice vector and set the steps
  int longest_r_index=std::distance(std::begin(lattice_vector_length), std::max_element(lattice_vector_length.begin(), lattice_vector_length.end()));
  int k_points_long_r=floor(ksampling[longest_r_index]/2.0);
  // no reason to use less points than the original set, so set this
  // now set the rest, normalized to the longest lattice vector stepping
  // we require that we at least integrate as many k-points
  // as supplied, since we develop this in a sphere we need to account for
  // the difference
  // also the number of kpoints along the longest direction determines the
  // others by the scaling relations above
  //
  // since we generate data in a sphere, we also account for the difference between
  // the volume of a cube and a sphere a factor of (pi/6)^(1/3) and ceil this value
  //
  // then we multiply by star_radius_factor which expands the "virtual" radius and
  // basically multiplies the number of kpoints above by its factor along each
  // direction
  double vol_fact=pow(6.0/M_PI,1.0/3.0);
  double factor=pow(star_radius_factor,1.0/3.0);
  rmax[longest_r_index]=ceil(factor*vol_fact*k_points_long_r);
  if (longest_r_index==0) {
    rmax[1]=ceil(factor*vol_fact*k_points_long_r*lattice_vector_length[0]/lattice_vector_length[1]);
    rmax[2]=ceil(factor*vol_fact*k_points_long_r*lattice_vector_length[0]/lattice_vector_length[2]);
  }
  else if (longest_r_index==1) {
    rmax[0]=ceil(factor*vol_fact*k_points_long_r*lattice_vector_length[1]/lattice_vector_length[0]);
    rmax[2]=ceil(factor*vol_fact*k_points_long_r*lattice_vector_length[1]/lattice_vector_length[2]);
  }
  else {
    rmax[0]=ceil(factor*vol_fact*k_points_long_r*lattice_vector_length[2]/lattice_vector_length[0]);
    rmax[1]=ceil(factor*vol_fact*k_points_long_r*lattice_vector_length[2]/lattice_vector_length[1]);
  }
  INFO_DUMP(std::cout << "SKW: starting lattice vector expansion: " << rmax[0] << ' ' << rmax[1] << ' ' << rmax[2] << ", in each direction, in total: " << (2*rmax[0]+1)*(2*rmax[1]+1)*(2*rmax[2]+1) << std::endl);

  int index=0;
  // tranpose the lattice
  transpose_matrix_3x3(lattice,lattice_transposed);
  // loop and store points within the ellipsoid to set the index of the
  // lattice point, also calculate length and store for later sorting
  double outer_radius=rmax[longest_r_index]*lattice_vector_length[longest_r_index];
  for (int i=-rmax[0];i<=rmax[0];i++) {
    for (int j=-rmax[1];j<=rmax[1];j++) {
      for (int k=-rmax[2];k<=rmax[2];k++) {
	std::vector<double> ttemp_vector={(double)i,(double)j,(double)k};
	multiply_vector_matrix_3x3_d(lattice_transposed,ttemp_vector,temp_vector);
	// if inside radius, else, skip point
	double radius=length_of_vector_d(temp_vector);
	if (radius < outer_radius) {  
	  r_tmp.push_back(std::vector<int>(3));
	  r_tmp[index][0]=i;
	  r_tmp[index][1]=j;
	  r_tmp[index][2]=k;
	  r_length.push_back(radius);
	  DEBUG_DUMP_LATPOINTS(std::cout << "adding lattice vector (radius of " << radius << " < " << factor*lattice_vector_length[longest_r_index] << " : " << i << ' ' << j << ' ' << k << std::endl);
	  index++;
	}
      }
    } 
  }
  INFO_DUMP(std::cout << "SKW: total number of points remaining after radius cutout: " << index << std::endl);

  // now sort on radius before checking symmetrized vectors (simple, but this routine requires C++11 std) and check shells for symmetry
  for (auto i:sort_indexes(r_length)) {
    r_index.push_back(i);
  }
  // set first shell point and loop over shells
  r.push_back(r_tmp[r_index[0]]);
  rl.push_back(r_length[r_index[0]]);
  for (int i=1;i<r_index.size();i++) {
    // default, store vector
    bool add_vector=true;
    // if similar length, check symmetries and if the lattice vector already
    // exist in the array, otherwise also add
    if (abs(r_length[r_index[i]]-rl[rl.size()-1]) < constants::ZERO) {
      // lattice vector in same shell as previous vector, call symmetry
      // routines to check all lattice vectors spanned by point operations
      sym_vectors.clear();
      locate_symmetry_vectors(r_tmp[r_index[i]],sym_ops,sym_vectors);
      // compare all symmetrized lattice vectors to the storred lattice 
      // vectors
      for (int j=0;j<r.size() && add_vector;j++) {
	for (int k=0;k<sym_vectors.size();k++) {
	  int diffx=abs(r[j][0]-sym_vectors[k][0]);
	  int diffy=abs(r[j][1]-sym_vectors[k][1]);
	  int diffz=abs(r[j][2]-sym_vectors[k][2]);
	  // found lattice vector
	  if ((diffx==0) && (diffy==0) && (diffz==0)) {
	    add_vector=false;
	    break;
	  }
	}
      }
    }
    if (add_vector==true) {
      r.push_back(r_tmp[r_index[i]]);
      rl.push_back(r_length[r_index[i]]);
      INFO_DUMP(std::cout << "SWK: storing lattice vector: " << std::setw(2) << r_tmp[r_index[i]][0] << ' ' << std::setw(2) << r_tmp[r_index[i]][1] << ' ' << std::setw(2) << r_tmp[r_index[i]][2] << ", of length=" << r_length[r_index[i]] << std::endl);
    }
  }
  INFO_DUMP(std::cout << "SKW: total number of points remaining after symmetrization: " << r.size() << std::endl); 

  // now build jagged array for R vectors including symmetrization
  for (int latvec=0;latvec<r.size();latvec++) {
    sym_vectors.clear();
    locate_symmetry_vectors(r[latvec],sym_ops,sym_vectors);
    r_sym_span.push_back(sym_vectors);
  }
  
  INFO_DUMP(std::cout << "SKW: using a star radius factor of: " << star_radius_factor << " a total number of " << r.size() << " kpoints are going to be generated times symmetry equivalent vectors, in total " << index << " laid out in a radius. For reference, "  << num_kpoints << " kpoints were supplied." << std::endl);
  
}

// generate rho
void generate_rho(std::vector<std::vector<int> > &r, std::vector<double> &r_length, std::vector<double> &c, std::vector<double> &rho) {
  // set default c if not provided
  if (c.empty()) {
    set_c(c);
  }	
  // set rho[0] to a large numnber
  rho.push_back(1e12);
  // loop all lattice vectors
  for (int i=1;i<r.size();i++) {
    // r[0] is always zero by construction, so the next non-zero length r is r[1]
    double x2=pow(r_length[i]/r_length[1],2.0);
    rho.push_back(pow(1.0-c[1]*x2,2.0)+c[2]*pow(x2,3.0));
    DEBUG_DUMP_RHO(std::cout << x2 << ' ' << rho[i] << ' ' << r[i][0] << ' '<< r[i][1] << ' ' << r[i][2] << std::endl);
  }
}

void set_c(std::vector<double> &c) {
  // current default is 
  // c_0=0.0;
  // c_1=0.75;
  // c_2=0.75;
  c = {0.0,0.75,0.75};
}

// generate S
void generate_s(std::vector<std::vector<int> > &r, std::vector<std::vector<double> > &kpoints, std::vector<std::vector<std::vector<int> > > &sym_ops, std::vector<std::vector<double> > &s) {
  std::vector<std::vector<int> > sym_vectors;
  // allocate s
  s.resize(r.size(), std::vector<double> (kpoints.size()));
  // loop lattice vectors
  for (int latvec=0;latvec<r.size();latvec++) {
    // fetch all non-equivalent vectors by point group operations
    sym_vectors.clear();
    locate_symmetry_vectors(r[latvec],sym_ops,sym_vectors);
    DEBUG_DUMP_S(std::cout << "for r: " << r[latvec][0] << ' ' << r[latvec][1] << ' ' << r[latvec][2] << ", " << sym_vectors.size() << " non-equivalent vectors where found by " << sym_ops.size() << " point group operations" << std::endl);
    for (int kpoint=0;kpoint<kpoints.size();kpoint++) {
      std::complex<double> sum=(0.0,0.0);
      for (int sym=0;sym<sym_vectors.size();sym++) {
	std::complex<double> itwopi(0.0,constants::TWOPI);
	sum+=exp(itwopi*dot_vectors(kpoints[kpoint],sym_vectors[sym]));
      }
      s[latvec][kpoint]=sum.real()/sym_vectors.size();
      DEBUG_DUMP_S(std::cout << "latvec " << latvec << ": " << r[latvec][0] << ' ' << r[latvec][1] << ' ' << r[latvec][2] << ", kpoint " << kpoint << ":, " << kpoints[kpoint][0] << ' ' << kpoints[kpoint][1] << ' ' << kpoints[kpoint][2] << ", S: " << s[latvec][kpoint] << std::endl);
    }
  }
  INFO_DUMP(std::cout << "SKW: A total of " << kpoints.size() << " input kpoints were used to construct the coefficients." << std::endl);
}

void generate_h_and_ssr(std::vector<std::vector<double> > &s, std::vector<double> &rho, std::vector<double> &h, std::vector<std::vector<double> > &ssr) {
  int num_kpoints=s[0].size();
  std::vector<std::vector<double> > s1(s.size(),std::vector<double> (num_kpoints-1));
  std::vector<std::vector<double> > s2(num_kpoints,std::vector<double> (s.size()));
  std::vector<std::vector<double> > h_temp(num_kpoints-1,std::vector<double> (num_kpoints-1));

  /* THIS SETS UP H USING DGEMM IN LAPACK. TOO COMPLICATED?

   */
  /*
  std::vector<std::vector<double> > s1(num_kpoints,std::vector<double> (s.size()));
  std::vector<std::vector<double> > s2(s.size(),std::vector<double> (num_kpoints));
  h.resize(num_kpoints*num_kpoints);
  // first do S difference (no complex conjugate), store column order
  for (int i=0;i<num_kpoints;i++) {
    for (int latvec=0;latvec<s.size();latvec++) {
      s1[i][latvec]=s[latvec][i]-s[latvec][num_kpoints-1];
      //std::cout << "s1: " << i << ' ' << latvec << ' ' << s1[i][latvec] << std::endl;
    }
  }
  // now do second S difference with complex conjugate (transpose, cycle index)
  for (int latvec=0;latvec<s.size();latvec++) {
    for (int i=0;i<num_kpoints;i++) {
      s2[latvec][i]=s1[i][latvec]/rho[latvec];
      //std::cout << "s2: " << latvec << ' ' << i << ' ' << s2[latvec][i] << std::endl;
    } 
  }
  // use LAPACK dgemm
  char trans='N';
  double alpha=1.0;
  double beta=1.0;
  // LAPACK interface only accepts flat arrays
  std::vector<double> h_test(num_kpoints*num_kpoints);
  std::vector<double> s1_flat(num_kpoints*s.size());
  std::vector<double> s2_flat(num_kpoints*s.size());
  flatten_nxn_matrix(s1,s1_flat,0);
  flatten_nxn_matrix(s2,s2_flat,0);
  int num_latvecs=s.size();
  dgemm(&trans,&trans,&num_kpoints,&num_kpoints,&num_latvecs,&alpha,&s1_flat[0],&num_kpoints,&s2_flat[0],&num_latvecs,&beta,&h[0],&num_kpoints);
  */

  /* THIS SETS UP H WITH SIMPLE LOOPS. TOO EASY/SLOW?
     
     USING LAST POINT AS REFERENCE, GENERALIZE TO USE OTHER POINTS?
     
   */
  // set start of latvec, SKW recommends start_latvec=1, e.g. drop R_0.

  // first do S difference
  for (int latvec=0;latvec<s.size();latvec++) {
    for (int i=0;i<num_kpoints-1;i++) {
      s1[latvec][i]=s[latvec][i]-s[latvec][num_kpoints-1];
      DEBUG_DUMP_H(std::cout << "S1, latvec: " << latvec << ", kpoints: " << i << ", values: " << s1[latvec][i] << std::endl); 
    }
  }
  // now do second S difference (transpose, cycle index)
  for (int i=0;i<num_kpoints-1;i++) {
    for (int latvec=0;latvec<s.size();latvec++) {
      s2[i][latvec]=s1[latvec][i]/rho[latvec];
      DEBUG_DUMP_H(std::cout << "S2 (conj of S1), latvec: " << latvec << ", kpoints: " << i << ", values: " << s2[i][latvec] << std::endl);  
    } 
  }
  // generate H matrix
  for (int j=0;j<num_kpoints-1;j++) {
    for (int i=0;i<num_kpoints-1;i++) {
      double sum=0.0;
      for (int latvec=0;latvec<s.size();latvec++) {
	sum+=s1[latvec][j]*s2[i][latvec];
      }
      h_temp[j][i]=sum;
      DEBUG_DUMP_H(std::cout << "H(" << j << "," << i << "): " << sum << std::endl);
    }
  }
  // now flatten h and store it in column, fortran order
  h.resize((num_kpoints-1)*(num_kpoints-1));
  flatten_nxn_matrix(h_temp,h,0);

  // now set ssr (s star rho) for use in equation 12a
  ssr=s2;
}

void generate_energies_diff(std::vector<std::vector<double> > &energies, std::vector<double> &energies_diff) {
  int num_bands=energies.size();
  int num_kpoints=energies[0].size()-1;
  energies_diff.resize(num_bands*num_kpoints);
  for (int band=0,index=0;band<num_bands;band++) {
    for (int kpoint=0;kpoint<num_kpoints;kpoint++) {
      energies_diff[index]=energies[band][kpoint]-energies[band][num_kpoints];
      index++;
    }
  }
}

int generate_lambdas(std::vector<std::vector<double> > &energies, std::vector<double> &h, std::vector<double> &lambda) {
  char trans='N';
  int info;
  const int num_bands=energies.size();
  const int num_kpoints=energies[0].size()-1;
  std::vector<int> ipiv(num_kpoints);
  // assume h and lambda in fortran order
  // set up LU of H
  INFO_DUMP(std::cout << "SKW: Setting up LU factorization of H matrix." << std::endl);
  dgetrf(&num_kpoints,&num_kpoints,&h[0],&num_kpoints,&ipiv[0],&info);
  if (info!=0) {
    std::cout << "SKW: ERROR, INFO TAG from DGETRF: " << info << ". Something is not right with the H matrix. Please check your input files etc. Exiting." << std::endl;
    return 1;
  }
  // solve eq 10 in SKW 1988 paper
  INFO_DUMP(std::cout << "SKW: Solving equation 10 of SKW paper (PRB 38, 2722 (1988))" << std::endl);
  dgetrs(&trans,&num_kpoints,&num_bands, &h[0], &num_kpoints, &ipiv[0],&lambda[0],&num_kpoints,&info);
  if (info!=0) {
    std::cout << "SKW: ERROR, INFO TAG from DGETRS: " << info << ". Please check your input files etc. Exiting." << std::endl;
    return 1;
  }
  return 0;
}

void generate_epsilons(std::vector<std::vector<double> > &energies, std::vector<std::vector<double> > &s, std::vector<std::vector<double> > &ssr, std::vector<double> &lambda, std::vector<std::vector<double> > &epsilons) {
  // loop bands
  for (int band=0;band<energies.size();band++) {
    int num_kpoints=energies[band].size()-1;
    epsilons.push_back(std::vector<double> (s.size()));
    // for latvecs except the first, equation 12a we have
    double sum2=0.0;
    for (int latvec=1;latvec<s.size();latvec++) {
      double sum=0.0;
      for (int kpoint=0;kpoint<num_kpoints;kpoint++) {
	sum+=ssr[kpoint][latvec]*lambda[kpoint+num_kpoints*band];
      }
      epsilons[band][latvec]=sum;
      // this is for the first latvec, equation 12b
      sum2+=sum*s[latvec][num_kpoints];
    }
    epsilons[band][0]=energies[band][num_kpoints]-sum2;
  }
}

void check_interpolated_energies(std::vector<std::vector<double> > &energies, std::vector<std::vector<double> > &s, std::vector<std::vector<double> > &epsilon) {
  // generate interpolated energies
  for (int band=0;band<energies.size();band++) {
    for (int kpoint=0;kpoint<energies[band].size();kpoint++) {
      double sum=0.0;
      for (int latvec=0;latvec<s.size();latvec++) {
	sum+=epsilon[band][latvec]*s[latvec][kpoint];
      }
      if ((energies[band][kpoint]-sum)/energies[band][kpoint] > constants::TOL) {
	std::cout << "SKW WARNING: The interpolated energies differ by more than a relative error of " << constants::TOL << "(" << (energies[band][kpoint]-sum)/energies[band][kpoint] << ") compared to the original energies for kpoint: " << kpoint << ", and band: " << band << ". Try to use more lattice vectors in the expansion to avoid this problem." << std::endl;
      }

    }
  }
}

void interpolate(std::vector<std::vector<double> > &epsilons, std::vector<std::vector<std::vector<int> > > &r, std::vector<std::vector<double> > &lattice, std::vector<std::vector<double> > &energies, std::vector<std::vector<std::vector<double> > > &velocities, std::vector<std::vector<double> > &kpoints, std::vector<int> &ksampling) {
  // calculate size of fft dimensions (crude, take max-min dimension in each direction)
  // TODO: Consider checking different elements along each axis
  // TODO: thus there is probably some overhead here...
  int test=0;
  for (int latvec=0;latvec<r.size();latvec++) {
    test+=r[latvec].size();
  }
  int n0_min=r[0][0][0], n1_min=r[0][0][1], n2_min=r[0][0][2];
  int n0_max=n1_min, n1_max=n1_min, n2_max=n2_min;
  for (int latvec=0;latvec<r.size();latvec++) {
    for (int symop=0;symop<r[latvec].size();symop++) {
      if (r[latvec][symop][0] < n0_min) {
	n0_min=r[latvec][symop][0];
      }
      if (r[latvec][symop][1] < n1_min) {
	n1_min=r[latvec][symop][1];
      }
      if (r[latvec][symop][2] < n2_min) {
	n2_min=r[latvec][symop][2];
      }
      if (r[latvec][symop][0] > n0_max) {
	n0_max=r[latvec][symop][0];
      }
      if (r[latvec][symop][1] > n1_max) {
	n1_max=r[latvec][symop][1];
      }
      if (r[latvec][symop][2] > n2_max) {
	n2_max=r[latvec][symop][2];
      }      
    }
  }
  // assume symmetric layout around zero, add zero
  int n0=n0_max-n0_min+1, n1=n1_max-n1_min+1, n2=n2_max-n2_min+1;
  int n=n0*n1*n2;

  int num_bands=epsilons.size();
  
  // resize vectors
  kpoints.resize(n,std::vector<double> (3));
  energies.resize(num_bands,std::vector<double> (n));
  velocities.resize(num_bands, std::vector<std::vector<double> > (3, std::vector<double> (n)));
  
  // allocate some FFTW arrays
  fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
  fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

  // now fill input array with the epsilons for one test band
  // later add loop for band
  for (int band=0;band<num_bands;band++) {
    // do energies and velocities, energies slotted at dir=3
    for (int dir=3;dir>=0;dir--) {
      for (int i=0;i<n;i++) {
	// not really necessary for the real part, still do it
	in[i][0]=0.0;
	in[i][1]=0.0;
      }
      for (int latvec=0;latvec<r.size();latvec++) {
	int num_sym_ops=r[latvec].size();
	for (int symop=0;symop<num_sym_ops;symop++) {
	  // fetch index
	  int i=r[latvec][symop][0];
	  int j=r[latvec][symop][1];
	  int k=r[latvec][symop][2];
	  // fetch cartesian values (for the velocities etc.)
	  std::vector<double> r_cart(3);
	  std::vector<double> r_temp={(double)i, (double)j, (double)k};
	  multiply_vector_matrix_3x3_d(lattice,r_temp,r_cart);
	  // make sure negative half plane is stored at the end
	  if (i<0) {
	    i+=n0;
	  }
	  if (j<0) {
	    j+=n1;
	  }
	  if (k<0) {
	    k+=n2;
	  }
	  int index=k+n2*(j+n1*i);
	  if (dir==3) {
	    in[index][0]=epsilons[band][latvec]/num_sym_ops;
	    in[index][1]=0.0;
	  }
	  else {
	    // velocities are complex due to the derivative of
	    // the plane wave expansion
	    // r is in direct coordinates, transform to cartesian
	    in[index][0]=0.0;
	    in[index][1]=r_cart[dir]*epsilons[band][latvec]/num_sym_ops;
	  }
	}
      }
      fftw_plan p = fftw_plan_dft_3d(n0,n1,n2,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
      fftw_execute(p);
      
      int starti=ceil((double)n0/2);
      int startj=ceil((double)n1/2);
      int startk=ceil((double)n2/2);

      for (int i=starti,index=0,ii=0;i<n0+starti;i++) {
	if (i>n0-1) {
	  ii=i-n0;
	}
	else {
	  ii=i;
	}
	for (int j=startj,jj=0;j<n1+startj;j++) {
	  if (j>n1-1) {
	    jj=j-n1;
	  }
	  else {
	    jj=j;
	  }
	  for (int k=startk,kk=0;k<n2+startk;k++) {
	    if (k>n2-1) {
	      kk=k-n2;
	    }
	    else {
	      kk=k;
	    }
	    int kpoint=kk+n2*(jj+n1*ii);
	    if (dir==3) {
	      // assume delta spacing for the lattice vector is always one
	      // same kpoint for each band
	      if (band==0) {
		kpoints[index][0]=(double)(i-n0)/n0;
		kpoints[index][1]=(double)(j-n1)/n1;
		kpoints[index][2]=(double)(k-n2)/n2;
	      }
	      energies[band][index]=out[kpoint][0];
	      DEBUG_DUMP_FFT(std::cout << "SKW: interpolated energy grid, band: " << band << ", kpoint: " << std::setw(10) << kpoints[index][0] << ' ' << std::setw(10) << kpoints[index][1] << ' ' << std::setw(10) << kpoints[index][2] << ", energy: " << energies[band][index] << std::endl);
	    }
	    else {
	      velocities[band][dir][index]=out[kpoint][0];
	      DEBUG_DUMP_FFT(std::cout << "SKW: interpolated velocity grid, band: " << band << ", kpoint: " << std::setw(10) << kpoints[index][0] << ' ' << std::setw(10) << kpoints[index][1] << ' ' << std::setw(10) << kpoints[index][2] << ", velocities(x,y,z): " << velocities[band][0][index] << ' ' << velocities[band][1][index] << ' ' << velocities[band][2][index] << std::endl);
	    }
	    index++;
	  }
	}
      }
      fftw_destroy_plan(p);
    }
  }
  fftw_free(in);
  fftw_free(out);
  
  // store ksamplings and return
  ksampling[0]=n0;
  ksampling[1]=n1;
  ksampling[2]=n2;
  INFO_DUMP(std::cout << "SKW: A total of " << kpoints.size() << " kpoints for each band were extracted from the interpolation routine. " << std::endl);
}

void locate_symmetry_vectors(std::vector<int> &vec, std::vector<std::vector<std::vector<int> > > &sym_ops, std::vector<std::vector<int> > &sym_vectors) {
  // first add input vector
  sym_vectors.push_back(vec);
  DEBUG_DUMP_SYMOPS(std::cout << __FUNCTION__ << ": checking lattice vector " << vec[0] << " " << vec[1] << " " << vec[2] << " for symmetry, applying " << sym_ops.size() << " operations" << std::endl);
  for (int sym=0;sym<sym_ops.size();sym++) {
    std::vector<int> sym_vector(3,0);    
    multiply_vector_matrix_3x3_i(sym_ops[sym],vec,sym_vector);
    // loop previous vectors and check if it exists
    int num_sym_vectors=sym_vectors.size();
    bool found_vector=false;
    for (int i=0;i<num_sym_vectors;i++) {
      if ((abs(sym_vector[0]-sym_vectors[i][0]) == 0) && (abs(sym_vector[1]-sym_vectors[i][1]) == 0) && (abs(sym_vector[2]-sym_vectors[i][2]) == 0)) {
	found_vector=true;
	break;
      }
    }
    if (!found_vector) {
      sym_vectors.push_back(sym_vector);
      DEBUG_DUMP_SYMOPS(std::cout << __FUNCTION__ << ": did not find rotated vector in database, adding vector " << sym_vector[0] << ' ' << sym_vector[1] << ' ' << sym_vector[2]  << ", now in total " << sym_vectors.size() << " exist in the database" << std::endl);
    }
  }
}


void transpose_matrix_3x3(std::vector<std::vector<double> > &matrix, std::vector<std::vector<double> > &transposed_matrix) {
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      transposed_matrix[j][i]=matrix[i][j];
    }
  }
}

void multiply_vector_matrix_3x3_d(std::vector<std::vector<double> > &matrix, std::vector<double> &vec, std::vector<double> &return_vec) {
  for (int i=0;i<3;i++) {
    double sum=0.0;
    for (int j=0;j<3;j++) {
      sum+=matrix[i][j]*vec[j];
    }
    return_vec[i]=sum;
  }
}

void multiply_vector_matrix_3x3_i(std::vector<std::vector<int> > &matrix, std::vector<int> &vec, std::vector<int> &return_vec) {
  for (int i=0;i<3;i++) {
    int sum=0;
    for (int j=0;j<3;j++) {
      sum+=matrix[i][j]*vec[j];
    }
    return_vec[i]=sum;
  }
}

void multiply_matrix_matrix_3x3_i(std::vector<std::vector<int> > &matrix1, std::vector<std::vector<int> > &matrix2, std::vector<std::vector<int> > &return_matrix) {
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      double sum=0.0;
      for (int k=0;k<3;k++) {
	sum+=matrix1[i][k]*matrix2[k][j];
      }
      return_matrix[i][j]=sum;
    }
  }
}

double dot_vectors(std::vector<double> &vec1, std::vector<int> &vec2) {
  double sum=0.0;
  for (int i=0;i<3;i++) {
    sum+=vec1[i]*(double)vec2[i];
  }
  return sum;
}

double length_of_vector_d(std::vector<double> &vec) {
  return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

double length_of_vector_i(std::vector<int> &vec) {
  return sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
}

void flatten_nxn_matrix(std::vector<std::vector<double> > &matrix, std::vector<double> &matrix_flat, int order) {
  // row, c order
  if (order) {
    for (int i=0, index=0; i<matrix.size();i++) {
      for (int j=0;j<matrix[i].size();j++) {
	matrix_flat[index]=matrix[i][j];
	index++;
      }
    }
  }
  // column, fortran order
  else {
    for (int i=0;i<matrix.size();i++) {
      for (int j=0;j<matrix[i].size();j++) {
	matrix_flat[j*matrix.size()+i]=matrix[i][j];
      }
    }
  }
}

void de_flatten_nxn_matrix(std::vector<double> &matrix_flat, std::vector<std::vector<double> > &matrix, int order) {
  // row, c order
  if (order) {
    for (int i=0, index=0; i<matrix.size();i++) {
      for (int j=0;j<matrix[i].size();j++) {
	matrix[i][j]=matrix_flat[0];
	index++;
      }
    }
  }
  // column, fortran order
  else {
    for (int i=0;i<matrix.size();i++) {
      for (int j=0;j<matrix[i].size();j++) {
	matrix[i][j]=matrix_flat[j*matrix.size()+i];
      }
    }
  } 
}
