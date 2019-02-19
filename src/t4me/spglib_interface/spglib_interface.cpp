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

#include "spglib_interface.hpp"

/* Interface to spglib

   We can not call spglib directly using numpy arrays through
   Cython. This interface thus copies the numpy content into new static
   C arrays and pass them to the spglib input arrays (where needed)

 */


void calc_dos_weights(double *energies_temp, int *spg_grid_address_temp,
                      int *mapping_bz_to_ibz_temp,
                      int *mapping_ibz_to_bz_temp,
                      int *ibz_weights_temp,
                      int *mesh_temp, double *rec_basis_temp,
                      double *energy_samples, int num_energy_samples,
                      int num_bands, int num_kpoints_ibz,
                      int *spin_degen_temp,
                      double volume, double volume_bz,
                      int weight_type,
                      double smearing,
                      double *dos, double *int_dos) {

  // set to false if the full BZ sum is wanted, in fact, there are
  // subtle differences between these two approaches, default is to run
  // with IBZ on, since the computational performance is better
  bool ibz=true;

  std::vector<int> mesh(3);
  for (auto sampling=0;sampling<3;sampling++) {
    mesh[sampling] = mesh_temp[sampling];
  }
  int num_kpoints=std::accumulate(std::begin(mesh),std::end(mesh), 1, std::multiplies<int>());

  std::vector<int> mapping_bz_to_ibz(num_kpoints);
  std::vector<int> mapping_ibz_to_bz(num_kpoints_ibz);
  std::vector<int> ibz_weights(num_kpoints_ibz);
  for (auto kpoint=0; kpoint < num_kpoints; kpoint++) {
    mapping_bz_to_ibz[kpoint] = mapping_bz_to_ibz_temp[kpoint];
    if (kpoint < num_kpoints_ibz) {
      mapping_ibz_to_bz[kpoint] = mapping_ibz_to_bz_temp[kpoint];
      ibz_weights[kpoint] = ibz_weights_temp[kpoint];
    }
  }

  std::vector<std::vector<double> >
    energies(num_bands,
             std::vector<double> (num_kpoints_ibz));

  std::vector<int> spin_degen(num_bands);
  for (auto band=0;band<num_bands;band++) {
    spin_degen[band] = spin_degen_temp[band];
    for (auto kpoint=0;kpoint<num_kpoints_ibz;kpoint++) {
      energies[band][kpoint] = *((double *)energies_temp +
                                 band*num_kpoints_ibz + kpoint);
    }
  }

  // needs to be array due to call to SPGLIB
  double rec_basis[3][3];
  for (auto i=0;i<3;++i) {
    for (auto j=0;j<3;j++) {
      rec_basis[i][j] = *((double *)rec_basis_temp + i*3 + j);
    }
  }

  std::vector<std::vector<int> >
    spg_grid_address(num_kpoints, std::vector<int> (3));
  for (auto kpoint=0;kpoint<num_kpoints;kpoint++) {
    for (auto dir=0;dir<3;dir++) {
      spg_grid_address[kpoint][dir] =
	*((int *)spg_grid_address_temp + kpoint*3 + dir);
    }
  }

  // determine max cantor array size
  int max_index = cantor_transform_3_i(mesh);
  std::vector<int> spg_grid_cantor_address(max_index);

  // cantor does not support negative values (at least not
  // the function implemented here)
  // from SPGLIB the points are always
  // -floor(mesh[i]/2), floor(mesh[i]/2) (odd)
  // -floor(mesh[i]/2), floor(mesh[i]/2)+1 (even with right adjust, default old version)
  // -floor(mesh[i]/2)-1, floor(mesh[i]/2) (even with left adjust)
  // so right now we only shift all points by +floor(mesh[i]/2)
  // and calculate the cantor value
  for (auto kpoint=0;kpoint<num_kpoints;kpoint++) {
    std::vector<int> gp(3);
    for (auto dir=0;dir<3;dir++) {
      gp[dir] = spg_grid_address[kpoint][dir]+floor(mesh[dir]/2);
    }
    int cantor_index = cantor_transform_3_i(gp);
    // we want to be able to fetch the ibz index straight from the cantor index
    spg_grid_cantor_address[cantor_index] = mapping_bz_to_ibz[kpoint];
  }

  // NOW EVERYTHING (WELL, ALMOST) IS STL AND WE CONTINUE WITH THE
  // INTERESTING STUFF

  // fetch weights
  std::vector<std::vector<std::vector<double> > >
    weights(num_energy_samples, std::vector<std::vector<double> >
	    (num_bands, std::vector<double> (num_kpoints_ibz)));

  std::vector<std::vector<std::vector<double> > >
    int_weights(num_energy_samples, std::vector<std::vector<double> >
		(num_bands, std::vector<double> (num_kpoints_ibz)));

  if (weight_type == 0) {
    tetra_weights(energies,spg_grid_address,spg_grid_cantor_address, mapping_bz_to_ibz,mapping_ibz_to_bz,
                  ibz_weights,mesh,rec_basis,energy_samples,num_bands,
                  num_kpoints_ibz, ibz,spin_degen,
                  weights,int_weights);
  }
  else if (weight_type == 2) {
    smearing_weights(energies,ibz_weights,energy_samples,
                     num_energy_samples,num_kpoints,smearing,ibz,
                     spin_degen,weights,int_weights);
  }
  else {
    std::cout << "The supplied weight_type is not supported." << std::endl;
  }

  // sum the weights over the IBZ
  double scaling = 1.0/volume;
  if (ibz) {
    for (int energy=0; energy < num_energy_samples; energy++) {
      for (int band=0; band<num_bands; band++) {
	dos[band*num_energy_samples+energy]=0.0;
	int_dos[band*num_energy_samples+energy]=0.0;
	for (int kpoint=0; kpoint < num_kpoints_ibz; kpoint++) {
	  dos[band*num_energy_samples+energy] +=
	    weights[energy][band][kpoint]*scaling;
	  int_dos[band*num_energy_samples+energy] +=
	    int_weights[energy][band][kpoint]*scaling;
	}
      }
    }
  }
  // sum the weights over the BZ
  else {
    for (int energy = 0; energy < num_energy_samples; energy++) {
      for (int band=0; band<num_bands; band++) {
	dos[band*num_energy_samples+energy]=0.0;
	int_dos[band*num_energy_samples+energy]=0.0;
	for (int kpoint = 0; kpoint < num_kpoints; kpoint++) {
	  dos[band*num_energy_samples+energy] +=
	    weights[energy][band][mapping_bz_to_ibz[kpoint]]*scaling;
	  int_dos[band*num_energy_samples+energy] +=
	    int_weights[energy][band][mapping_bz_to_ibz[kpoint]]*scaling;
	}
      }
    }
  }
}

void calc_transport_tensors_weights(double *energies_temp,
				  double *velocities_temp,
				  double *scattering_temp,
				  double *temperatures_temp,
				  double *chempots_temp,
				  int *spg_grid_address_temp,
				  int *mapping_bz_to_ibz_temp,
				  int *mapping_ibz_to_bz_temp,
                                  int *ibz_weights_temp,
				  int *mesh_temp, double *rec_basis_temp,
				  int num_bands, int num_kpoints_ibz,
				  int num_temp_samples,
				  int num_chempot_samples,
				  int int_method,
				  int num_energy_samples,
                                  int weight_type,
                                  double smearing,
				  double energy_cutoff,
				  double volume,
				  int *spin_degen_temp,
				  double *cond,
				  double *seebeck,
				  double *lorenz) {


  // Now convert everything to STL vectors. This might seem strange
  // but is just a consequence of Cythons limited Numpy -> STL mapping
  // In the future, consider writing custom interface using the
  // Numpy C-API. Also, in the future, SPGLIB stuff will most likely
  // move to C++ and thus also use STL vectors/arrays
  // ALSO WE DO NOT REALLY NEED TO CONVERT EVERYTHING, BUT WE JUST DO
  // IT FOR CONSISTENCY HERE...

  // For the future also consider std::unique_ptr
  std::vector<int> mesh(3);
  for (auto sampling=0;sampling<3;sampling++) {
    mesh[sampling] = mesh_temp[sampling];
  }
  int num_kpoints=std::accumulate(std::begin(mesh),std::end(mesh), 1, std::multiplies<int>());

  std::vector<double> temperatures(num_temp_samples);
  for (auto temperature=0;temperature<num_temp_samples;temperature++) {
    temperatures[temperature] = temperatures_temp[temperature];
  }

  std::vector<double> chempots(num_chempot_samples);
  for (auto chempot=0;chempot<num_chempot_samples;chempot++) {
    chempots[chempot] = chempots_temp[chempot];
  }

  std::vector<int> mapping_bz_to_ibz(num_kpoints);
  std::vector<int> mapping_ibz_to_bz(num_kpoints_ibz);
  std::vector<int> ibz_weights(num_kpoints_ibz);
  for (auto kpoint=0; kpoint < num_kpoints; kpoint++) {
    mapping_bz_to_ibz[kpoint] = mapping_bz_to_ibz_temp[kpoint];
    if (kpoint < num_kpoints_ibz) {
      mapping_ibz_to_bz[kpoint] = mapping_ibz_to_bz_temp[kpoint];
      // we have to run the spectral sum over the whole BZ, thus
      // set the weights to one
      //ibz_weights[kpoint] = ibz_weights_temp[kpoint];
      ibz_weights[kpoint] = 1.0;
    }
  }

  std::vector<std::vector<double> >
    energies(num_bands,
             std::vector<double> (num_kpoints_ibz));

  std::vector<std::vector<std::vector<double> > > scattering;

  std::vector<std::vector<std::vector<double> > >
    velocities(num_bands,std::vector<std::vector<double> >
               (num_kpoints, std::vector<double> (3)));

  // we use an energy interval that is energy_cutoff eV on both sides of
  // the chempot range with num_energy_samples sampling
  double energy_samples[num_energy_samples];
  fetch_energy_samples(chempots,num_chempot_samples,energy_cutoff,
		       energy_samples,num_energy_samples);

  // resize scattering array
  scattering.resize(num_temp_samples, std::vector<std::vector<double> > (num_bands, std::vector<double> (num_kpoints)));

  std::vector<int> spin_degen(num_bands);
  for (auto band=0;band<num_bands;band++) {
    spin_degen[band] = spin_degen_temp[band];
    for (auto kpoint=0;kpoint<num_kpoints_ibz;kpoint++) {
      energies[band][kpoint] = *((double *)energies_temp +
				      band*num_kpoints_ibz + kpoint);
    }
    for (auto kpoint=0;kpoint<num_kpoints;kpoint++) {
      // set scattering for kpoints
      for (auto temp=0;temp<num_temp_samples;temp++) {
        scattering[temp][band][kpoint]=scattering_temp[temp*num_bands*num_kpoints+band*num_kpoints+kpoint];
      }
      for (auto dir=0;dir<3;dir++) {
	velocities[band][kpoint][dir] =
	  *((double *)velocities_temp +
	    band*num_kpoints*3+num_kpoints*dir+kpoint);
      }
    }
  }

  // needs to be array due to call to SPGLIB
  double rec_basis[3][3];
  for (auto i=0;i<3;++i) {
    for (auto j=0;j<3;j++) {
      rec_basis[i][j] = *((double *)rec_basis_temp + i*3 + j);
    }
  }

  std::vector<std::vector<int> >
    spg_grid_address(num_kpoints, std::vector<int> (3));
  for (auto kpoint=0;kpoint<num_kpoints;kpoint++) {
    for (auto dir=0;dir<3;dir++) {
      spg_grid_address[kpoint][dir] =
	*((int *)spg_grid_address_temp + kpoint*3 + dir);
    }
  }

  // determine max cantor array size
  int max_index = cantor_transform_3_i(mesh);
  std::vector<int> spg_grid_cantor_address(max_index);

  // cantor does not support negative values (at least not
  // the function implemented here)
  // from SPGLIB the points are always
  // -floor(mesh[i]/2), floor(mesh[i]/2) (odd)
  // -floor(mesh[i]/2), floor(mesh[i]/2)+1 (even with right adjust, default old version)
  // -floor(mesh[i]/2)-1, floor(mesh[i]/2) (even with left adjust)
  // so right now we only shift all points by +floor(mesh[i]/2)
  // and calculate the cantor value
  for (auto kpoint=0;kpoint<num_kpoints;kpoint++) {
    std::vector<int> gp(3);
    for (auto dir=0;dir<3;dir++) {
      gp[dir] = spg_grid_address[kpoint][dir]+floor(mesh[dir]/2);
    }
    int cantor_index = cantor_transform_3_i(gp);
    // we want to be able to fetch the ibz index straight from the cantor index
    spg_grid_cantor_address[cantor_index] = mapping_bz_to_ibz[kpoint];
  }

  // NOW EVERYTHING (WELL, ALMOST) IS STL AND WE CONTINUE WITH THE
  // INTERESTING STUFF

  // fetch weights
  std::vector<std::vector<std::vector<double> > >
    weights(num_energy_samples, std::vector<std::vector<double> >
	    (num_bands, std::vector<double> (num_kpoints_ibz)));

  std::vector<std::vector<std::vector<double> > >
    int_weights(num_energy_samples, std::vector<std::vector<double> >
		(num_bands, std::vector<double> (num_kpoints_ibz)));

    // first, fetch the weights for the IBZ
  if (weight_type == 0) {
    // now we perform the tetrahedron integration basically: 1. fetch
    // kpoint weights for the IBZ cell 2. multiply weights with the
    // function X, see PRB 49, 16223, 1994
    tetra_weights(energies, spg_grid_address,
                  spg_grid_cantor_address,
                  mapping_bz_to_ibz,
                  mapping_ibz_to_bz, ibz_weights,
                  mesh,rec_basis,energy_samples,
                  num_bands,num_kpoints_ibz,false,spin_degen,
                  weights,int_weights);
  }
  else if (weight_type == 2) {
    // just a Gaussian smearing function
    smearing_weights(energies,ibz_weights,
                     energy_samples,num_energy_samples,num_kpoints,
                     smearing,false,spin_degen,weights,int_weights);
  }
  else {
    std::cout << "The supplied weight_type " << weight_type << " is not supported. Exiting." << std::endl;
    exit(1);
  }

  // do the integral (basically just the multiplication of weights and
  // input function) we then obtain the spectral function at a given
  // energy for all the transport coefficients
  std::vector<std::vector<std::vector<std::vector<double> > > >
    sigma_k2(num_temp_samples, std::vector<std::vector<std::vector<double> > >
	     (num_energy_samples, std::vector<std::vector<double> >
	      (3,std::vector<double> (3))));

  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >
    sigma_k3(num_temp_samples, std::vector<std::vector<std::vector<std::vector<double> > > >
	     (num_energy_samples, std::vector<std::vector<std::vector<double> > >
	      (3, std::vector<std::vector<double> >
	       (3, std::vector<double> (3)))));
  spectral_sum(weights,int_weights,velocities,scattering,
	       mapping_bz_to_ibz,sigma_k2,sigma_k3);

  // now we need to perform the energy integral to obtain the transport
  // coefficient tensors for each temperature and chemical potential we
  // are now in principle done with the tetrahedron integration and only
  // need to perform an integral in energy space
  std::vector<std::vector<double> > sigma(3, std::vector<double> (3));
  std::vector<std::vector<std::vector<double> > >
    sigma_e3(3, std::vector<std::vector<double> >
	     (3, std::vector<double> (3)));
  std::vector<std::vector<double> > chi(3,std::vector<double> (3));
  std::vector<std::vector<double> > kappa(3,std::vector<double> (3));
  std::vector<std::vector<double> >
    sigma_inv(3,std::vector<double>
	      (3));
  std::vector<std::vector<double> > seebeck_temp(3,std::vector<double> (3));
  std::vector<std::vector<double> > lorenz_highdeg(3,std::vector<double> (3));
  std::vector<std::vector<double> > lorenz_cor(3,std::vector<double> (3));
  std::vector<std::vector<double> > tmp(3,std::vector<double> (3));
  // set units
  double sigmaunit=1e11*constants::PI*constants::G /
    (2*constants::HBAR*constants::KB*volume);
  double seebeckunit=1e6; double lorenzunit=1e8;
  for (int temp=0;temp<num_temp_samples;temp++) {
    for (int chempot=0;chempot<num_chempot_samples;chempot++) {
      e_integral(sigma_k2[temp],sigma_k3[temp],scattering,energy_samples,
		 temperatures[temp],chempots[chempot],int_method,
		 sigma,sigma_e3,chi,kappa);
      // invert the sigma matrix
      invert_3x3_matrix(sigma,sigma_inv);
      // now calculate the seebeck tensors and add units
      // also store the conductivity
      multiply_3x3_matrix(sigma_inv,chi,seebeck_temp);
      for (int dir1=0;dir1<3;dir1++) {
	for (int dir2=0;dir2<3;dir2++) {
          cond[temp*num_chempot_samples*9+9*chempot+3*dir1+dir2]=
	    sigmaunit*sigma[dir1][dir2]/temperatures[temp];
          seebeck[temp*num_chempot_samples*9+9*chempot+3*dir1+dir2]=-
	    seebeckunit*seebeck_temp[dir1][dir2]/temperatures[temp];
        }
      }
      // now calculate the lorenz tensor and add units
      // first store the highly degenerate part (first term)
      multiply_3x3_matrix(kappa,sigma_inv,lorenz_highdeg);
      for (int dir1=0;dir1<3;dir1++) {
	for (int dir2=0;dir2<3;dir2++) {
	  lorenz[temp*num_chempot_samples*9+9*chempot+3*dir1+dir2]=
	    lorenzunit*lorenz_highdeg[dir1][dir2]/pow(temperatures[temp],2.0);
	}
      }
      // now calculate the correction term and subtract
      multiply_3x3_matrix(chi,seebeck_temp,tmp);
      multiply_3x3_matrix(tmp,sigma_inv,lorenz_cor);
      for (int dir1=0;dir1<3;dir1++) {
	for (int dir2=0;dir2<3;dir2++) {
	  lorenz[temp*num_chempot_samples*9+9*chempot+3*dir1+dir2]-=
	    lorenzunit*lorenz_cor[dir1][dir2]/pow(temperatures[temp],2.0);
	}
      }

    }
  }
}
void e_integral(std::vector<std::vector<std::vector<double> > > &sigma_k2,
		std::vector<std::vector<std::vector<std::vector<double> > > > &sigma_k3,
		std::vector<std::vector<std::vector<double> > > &scattering,
		double *energy_samples,
		double temperature, double chempot, int int_method,
		std::vector<std::vector<double> > &sigma_e2,
		std::vector<std::vector<std::vector<double> > > &sigma_e3,
		std::vector<std::vector<double> > &chi_e,
		std::vector<std::vector<double> > &kappa_e) {

  double beta=1e5/(temperature*constants::KB);
  // simple sum rule integration (trapeziodal)
  if (int_method==0) {
    int num_energy_samples=sigma_k2.size();
    double spacing=(energy_samples[num_energy_samples-1]-
		    energy_samples[0])/(num_energy_samples-1);
    // interior
    for (int energy=1;energy<num_energy_samples-1;energy++) {
      for (int dir1=0;dir1<3;dir1++) {
	for (int dir2=0;dir2<3;dir2++) {
	  sigma_e2[dir1][dir2]+=
	    e_integrand(sigma_k2[energy][dir1][dir2],
            		energy_samples[energy],chempot,beta,0);
	  chi_e[dir1][dir2]+=
	    e_integrand(sigma_k2[energy][dir1][dir2],
			energy_samples[energy],chempot,beta,1);
	  kappa_e[dir1][dir2]+=
	    e_integrand(sigma_k2[energy][dir1][dir2],
			energy_samples[energy],chempot,beta,2);
	  for (int dir3=0;dir3<3;dir3++) {
	    sigma_e3[dir1][dir2][dir3]+=
	      e_integrand(sigma_k3[energy][dir1][dir2][dir3],
			  energy_samples[energy],chempot,beta,0);
	  }
	}
      }
    }
    // borders
    for (int dir1=0;dir1<3;dir1++) {
      for (int dir2=0;dir2<3;dir2++) {
    	sigma_e2[dir1][dir2]+=
	  (e_integrand(sigma_k2[0][dir1][dir2],
		       energy_samples[0],
		       chempot,beta,0)+
	   e_integrand(sigma_k2[num_energy_samples-1][dir1][dir2],
		       energy_samples[num_energy_samples-1],
		       chempot,beta,0))/2.0;
    	chi_e[dir1][dir2]+=
	  (e_integrand(sigma_k2[0][dir1][dir2],
		       energy_samples[0],
		       chempot,beta,1)+
	   e_integrand(sigma_k2[num_energy_samples-1][dir1][dir2],
		       energy_samples[num_energy_samples-1],
		       chempot,beta,1))/2.0;
    	kappa_e[dir1][dir2]+=
	  (e_integrand(sigma_k2[0][dir1][dir2],
		       energy_samples[0],chempot,beta,2)+
	   e_integrand(sigma_k2[num_energy_samples-1][dir1][dir2],
		       energy_samples[num_energy_samples-1],
		       chempot,beta,2))/2.0;
    	for (int dir3=0;dir3<3;dir3++) {
	  sigma_e3[dir1][dir2][dir3]+=
	    e_integrand(sigma_k3[0][dir1][dir2][dir3],
                        energy_samples[0],chempot,beta,0)+
            e_integrand(sigma_k3[num_energy_samples-1][dir1][dir2][dir3],
                        energy_samples[num_energy_samples],
                        chempot,beta,0);
    	}
      }
    }
    // scale
    for (int dir1=0;dir1<3;dir1++) {
      for (int dir2=0;dir2<3;dir2++) {
    	sigma_e2[dir1][dir2]=spacing*sigma_e2[dir1][dir2];
    	chi_e[dir1][dir2]=spacing*chi_e[dir1][dir2];
    	kappa_e[dir1][dir2]=spacing*kappa_e[dir1][dir2]; }
    } }
  // composite Simpson rule
  else if (int_method==1) {
    exit(1);
  }
  // and a interpolative adaptive Simpson
  else if (int_method==2) {
    exit(1);
  }
  // and a interpolative adaptive Cubature
  else if (int_method==3) {
    exit(1);
  }
  else {
    std::cout << "The supplied int_method is unavilable. Please supply a \
    valid int_method. Exiting." << std::endl; exit(1);
  }
}

double e_integrand(double f, double energy, double chempot, double
		   beta,int i) {
  return f*pow(energy-chempot,i)/dedf(energy,chempot,beta); }

double dedf(double energy, double chempot, double beta) {
  return 1.0+cosh((energy-chempot)*beta); }

void spectral_sum(std::vector<std::vector<std::vector<double> > >
		  &weights,
		  std::vector<std::vector<std::vector<double> > >
		  &int_weights,
		  std::vector<std::vector<std::vector<double> > >
		  &velocities,
		  std::vector<std::vector<std::vector<double> > > &scattering,
		  std::vector<int> &mapping_bz_to_ibz,
		  std::vector<std::vector<std::vector<std::vector<double> > > >
		  &sigma2,
		  std::vector<std::vector<std::vector<std::vector<std::vector<double> > > > >
		  &sigma3) {

  for (int temp=0;temp<scattering.size();temp++) {
    for (unsigned int energy=0;energy<weights.size();energy++) {
      for (unsigned int band=0;band<velocities.size();band++) {
	for (unsigned int kpoint=0;kpoint<velocities[0].size();kpoint++) {
	  for (int dir1=0;dir1<3;dir1++) {
	    for (int dir2=0;dir2<3;dir2++) {
              sigma2[temp][energy][dir1][dir2]+=
                scattering[temp][band][kpoint]*
                weights[energy][band][mapping_bz_to_ibz[kpoint]]*
                velocities[band][kpoint][dir1]*
                velocities[band][kpoint][dir2];
              for (int dir3;dir3<3;dir3++) {
                // not currently implemented
                sigma3[temp][energy][dir1][dir2][dir3]+=0.0;
              }
	    }
	  }
        }
      }
    }
  }
}

void smearing_weights(std::vector<std::vector<double> > &energies,
		      std::vector<int> &ibz_weights,
		      double *energy_samples, int num_energy_samples,
		      int num_kpoints, double sigma,
		      bool ibz, std::vector<int> &spin_degen,
		      std::vector<std::vector<std::vector<double> > >
		      &weights,
		      std::vector<std::vector<std::vector<double> > >
		      &int_weights) {

  // scaling
  double scaling=1.0/num_kpoints;
  // units of these weights are 1/energy_samples units
  for (unsigned int energy=0;energy<num_energy_samples;energy++) {
    for (unsigned int band=0;band<energies.size();band++) {
      for (unsigned int kpoint=0;kpoint<energies[band].size();kpoint++) {
	weights[energy][band][kpoint]=
	  spin_degen[band]*scaling*ibz_weights[kpoint]*
	  gaussian(energies[band][kpoint],
		   energy_samples[energy],sigma);
	int_weights[energy][band][kpoint]=0.0;
      }
    }
  }
}

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
		   std::vector<std::vector<std::vector<double> >
		   > &weights,
		   std::vector<std::vector<std::vector<double> > >
		   &int_weights) {

  unsigned int num_kpoints=mesh[0]*mesh[1]*mesh[2];
  unsigned int num_energy_samples=weights.size();

  // fetch relative grid addresses of the tetrahedrons around one point
  int relative_grid_address[24][4][3];
  thm_get_relative_grid_address(relative_grid_address,
					   rec_basis);

  // scaling
  double scaling=1.0/num_kpoints;
  // some temporaries
  double local_energies[24][4];
  std::vector<int> gp(3);
  // now find the integration weights
  for (unsigned int energy = 0;
       energy < num_energy_samples; energy++) {
    for (unsigned int band=0; band<num_bands; band++) {
      for (unsigned int kpoint = 0;
	   kpoint < num_kpoints_ibz; kpoint++) {
        for (unsigned int tetra = 0; tetra < 24; tetra++) {
          for (unsigned int corner = 0; corner < 4; corner++) {
            for (unsigned int dir = 0; dir < 3; dir++) {
              gp[dir] = spg_grid_address[mapping_ibz_to_bz[kpoint]][dir] +
                relative_grid_address[tetra][corner][dir];
              // cantor does not support negative numbers
              // -floor(mesh[i]/2), floor(mesh[i]/2) (odd)
              // -floor(mesh[i]/2), floor(mesh[i]/2)+1 (even with right adjust, default old version)
              // -floor(mesh[i]/2)-1, floor(mesh[i]/2) (even with left adjust)
              // so right now we only shift all points by +floor(mesh[i]/2)
              // and calculate the cantor value
	    }
            // shift back to cell and shift to avoid zeros
            for (unsigned int dir = 0; dir < 3; dir++) {
              back_to_cell(gp[dir], mesh[dir]);
              gp[dir] = gp[dir]+floor(mesh[dir]/2);
            }
            int cantor_index = cantor_transform_3_i(gp);
            // spg_grid_cantor_address[cantor_index] gives us the ibz index equivalent
            // to the gp point
            local_energies[tetra][corner] =
	      energies[band][spg_grid_cantor_address[cantor_index]];
          }
	}
	// the weights returned from the
	// spg_get_tetrahedra_integration_weights are in 1/energy_samples
	// units, e.g. 1/eV or 1/hbaromega, while the irreducible
	// weights are unitless. Thus the units of the total weights
	// would also be 1/energy_samples units
        weights[energy][band][kpoint] = spin_degen[band] * scaling *
	  ibz_weights[kpoint] *
	  thm_get_integration_weight(energy_samples[energy],
						local_energies, 'I');
	int_weights[energy][band][kpoint] =
	  spin_degen[band] * scaling * ibz_weights[kpoint] *
	  thm_get_integration_weight(energy_samples[energy],
						local_energies, 'J');
      }
    }
  }
}
void back_to_cell(int &gp, int &sampling) {
  double ratio = (double)gp / sampling;
  if (ratio > 0.5) {
    gp -= sampling;
  }
  else if (ratio < -0.5) {
    gp += sampling;
  }
}


static int grid_address_to_index(std::vector<int> &g, std::vector<int> &mesh) {
  return 0;
}

int cantor_transform_3_i(std::vector<int> &v) {
  int i = v[0];
  i+=((v[0]+v[1])*(v[0]+v[1]+1))/2;
  i+=((v[0]+v[1]+v[2])*(v[0]+v[1]+v[2]+1)*(v[0]+v[1]+v[2]+2))/6;
  return i;
}

void get_reciprocal_mesh_interface(int *mesh, double *lattice,
				   double *positions,
				   int *anumbers,
				   unsigned int num_atoms,
				   int *is_shifted,
				   int *mesh_points, int *mapping,
                                   char *intsym,
				   int is_time_reversal,
				   double symprec) {

  const unsigned int num_k_points=mesh[0]*mesh[1]*mesh[2];
  unsigned int i,j;

  // spglib does not use pointers or vectors, so convert 2D vectors to
  // std. c arrays (1d vector is compatible with a 1d array, so leave
  // those intact)
  int temp_mesh_points[num_k_points][3];
  for(i=0;i<num_k_points;++i) {
    for (j=0;j<3;j++) {
      temp_mesh_points[i][j] =
	*((int *)mesh_points + i*3 + j);
    }
  }
  double temp_lattice[3][3];
  for(i=0;i<3;++i) {
    for (j=0;j<3;j++) {
      temp_lattice[i][j]=
	*((double *)lattice + i*3 + j);
    }
  }
  double temp_positions[num_atoms][3];
  for(i=0;i<num_atoms;++i) {
    for (j=0;j<3;j++) {
      temp_positions[i][j]=
	*((double *)positions + i*3 + j);
    }
  }

  // fetch the international symmetry group
  spg_get_international(intsym, temp_lattice,
                             temp_positions, anumbers,
                             num_atoms, symprec);

  // then get the kpoint mesh and associated ibz through the mapping
  // table
  spg_get_ir_reciprocal_mesh(temp_mesh_points,mapping,mesh,
			     is_shifted,is_time_reversal,
			     temp_lattice,temp_positions,
			     anumbers,num_atoms,symprec);
  // copy pointers back
  for(i=0;i<num_k_points;++i) {
    for (j=0;j<3;j++) {
      *((int *)mesh_points + i*3 + j)=temp_mesh_points[i][j];
    }
  }
}

double gaussian(double energy, double energy_ref, double sigma) {
  return
    exp(-0.5*pow((energy_ref-energy)/sigma,2.0))/
    (sigma*sqrt(2*constants::PI));
}

void fetch_energy_samples(std::vector<double> &chempots,int
			  num_chempot_samples,double energy_cutoff,
			  double *energy_samples,
			  int num_energy_samples) {

  double start=chempots[0]-energy_cutoff;
  double end=chempots[num_chempot_samples-1]+energy_cutoff;

  double energy_step=(end-start)/
    (num_energy_samples-1);
  for (int i=0;i<num_energy_samples;i++) {
    energy_samples[i]=start+i*energy_step;
  }
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
