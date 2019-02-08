# Copyright 2016 Espen Flage-Larsen
#
#    This file is part of T4ME.
#
#    T4ME is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    T4ME is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with T4ME.  If not, see <http://www.gnu.org/licenses/>.

# distutils: language=c++

import cython
import logging
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "spglib_interface.hpp":
    void get_reciprocal_mesh_interface(int * mesh,
                                       double * lattice,
                                       double * positions,
                                       int * anumbers,
                                       unsigned int num_atoms,
                                       int * is_shifted,
                                       int * mesh_points,
                                       int * mapping,
                                       char * intsym,
                                       int is_time_reversal,
                                       double symprec)

    void calc_dos_weights(double * energies,
                          int * spg_grid_address,
                          int * mapping_bz_to_ibz,
                          int * mapping_ibz_to_bz,
                          int * ibz_weights,
                          int * mesh,
                          double * rec_basis,
                          double * energy_samples,
                          int num_energy_samples,
                          int num_bands, int num_kpoints_ibz,
                          int * spin_degen,
                          double volume, double volume_bz,
                          int weight_type,
                          double smearing,
                          double * dos, double * int_dos)

    void calc_transport_tensors_weights(double * energies,
                                        double * velocities,
                                        double * scattering,
                                        double * temperatures,
                                        double * chempots,
                                        int * spg_grid_address,
                                        int * mapping_bz_to_ibz,
                                        int * mapping_ibz_to_bz,
                                        int * ibz_weights,
                                        int * mesh, double * rec_basis,
                                        int num_bands,
                                        int num_kpoints_ibz,
                                        int num_temp_steps,
                                        int num_chempot_steps,
                                        int int_method,
                                        int energy_samples,
                                        int weight_type,
                                        double smearing,
                                        double chempot_cutoff,
                                        double volume,
                                        int * spin_degen,
                                        double * cond,
                                        double * seebeck,
                                        double * lorenz)


@cython.boundscheck(False)
@cython.wraparound(False)
def get_reciprocal_mesh(
        np.ndarray[int, ndim=1, mode="c"] mesh not None,
        np.ndarray[double, ndim=2, mode="c"] lattice not None,
        np.ndarray[double, ndim=2, mode="c"] positions not None,
        np.ndarray[int, ndim=1, mode="c"] anumbers not None,
        np.ndarray[int, ndim=1, mode="c"] is_shifted not None,
        np.ndarray[int, ndim=2, mode="c"] mesh_points not None,
        np.ndarray[int, ndim=1, mode="c"] mapping not None,
        is_time_reversal=True, symprec=1e-5):

    cdef char c_intsym[11]
    num_atoms = anumbers.shape[0]
    get_reciprocal_mesh_interface(& mesh[0],
                                   & lattice[0, 0],
                                   & positions[0, 0],
                                   & anumbers[0],
                                   num_atoms,
                                   & is_shifted[0],
                                   & mesh_points[0, 0],
                                   & mapping[0],
                                   & c_intsym[0],
                                   is_time_reversal,
                                   symprec)
    return c_intsym.decode()


def calc_density_of_states_interface(
        np.ndarray[double, ndim=2, mode="c"] energy not None,
        np.ndarray[int, ndim=2, mode="c"] grid_address not None,
        np.ndarray[int, ndim=1, mode="c"] mapping_bz_to_ibz not None,
        np.ndarray[int, ndim=1, mode="c"] mapping_ibz_to_bz not None,
        np.ndarray[int, ndim=1, mode="c"] ibz_weights not None,
        np.ndarray[int, ndim=1, mode="c"] mesh not None,
        np.ndarray[double, ndim=2, mode="c"] rec_basis not None,
        np.ndarray[double, ndim=1, mode="c"] energy_samples not None,
        num_energy_samples, num_bands, num_kpoints_ibz,
        spin_degen,
        volume, volume_bz,
        weight_type,
        smearing,
        np.ndarray[double, ndim=2, mode="c"] dos not None,
        np.ndarray[double, ndim=2, mode="c"] int_dos not None):

    # convert bool to int (cython<->with bool does not work very well)
    cdef np.ndarray[int, ndim = 1, mode = "c"] spin_degen_int = \
        spin_degen.astype(dtype="intc")
    calc_dos_weights( & energy[0, 0],
                     & grid_address[0, 0],
                     & mapping_bz_to_ibz[0],
                     & mapping_ibz_to_bz[0],
                     & ibz_weights[0],
                     & mesh[0],
                     & rec_basis[0, 0],
                     & energy_samples[0],
                     num_energy_samples, num_bands, num_kpoints_ibz,
                     & spin_degen_int[0],
                     volume, volume_bz,
                     weight_type,
                     smearing,
                     & dos[0, 0], &
                     int_dos[0, 0])


def calc_transport_tensors_weights_interface(
        np.ndarray[double, ndim=2, mode="c"] energies not None,
        np.ndarray[double, ndim=3, mode="c"] velocities not None,
        np.ndarray[double, ndim=3, mode="c"] scattering not None,
        np.ndarray[double, ndim=1, mode="c"] temperatures not None,
        np.ndarray[double, ndim=1, mode="c"] chempots not None,
        np.ndarray[int, ndim=2, mode="c"] grid_address not None,
        np.ndarray[int, ndim=1, mode="c"] mapping_bz_to_ibz not None,
        np.ndarray[int, ndim=1, mode="c"] mapping_ibz_to_bz not None,
        np.ndarray[int, ndim=1, mode="c"] ibz_weights not None,
        np.ndarray[int, ndim=1, mode="c"] mesh not None,
        np.ndarray[double, ndim=2, mode="c"] rec_basis not None,
        num_bands, num_kpoints_ibz, num_temp_steps,
        num_chempot_steps, int_method, energy_samples,
        weight_type, smearing,
        energy_cutoff, volume,
        spin_degen,
        np.ndarray[double, ndim=4, mode="c"] cond not None,
        np.ndarray[double, ndim=4, mode="c"] seebeck not None,
        np.ndarray[double, ndim=4, mode="c"] lorenz not None):

    # convert bool to int (cython<->with bool does not work very well)
    cdef np.ndarray[int, ndim = 1, mode = "c"] spin_degen_int = \
        spin_degen.astype(dtype="intc")
    calc_transport_tensors_weights(& energies[0, 0],
                                    & velocities[0, 0, 0],
                                    & scattering[0, 0, 0],
                                    & temperatures[0],
                                    & chempots[0],
                                    & grid_address[0, 0],
                                    & mapping_bz_to_ibz[0],
                                    & mapping_ibz_to_bz[0],
                                    & ibz_weights[0],
                                    & mesh[0],
                                    & rec_basis[0, 0],
                                    num_bands, num_kpoints_ibz,
                                    num_temp_steps,
                                    num_chempot_steps, int_method,
                                    energy_samples, weight_type,
                                    smearing,
                                    energy_cutoff, volume,
                                    & spin_degen_int[0],
                                    & cond[0, 0, 0, 0],
                                    & seebeck[0, 0, 0, 0],
                                    & lorenz[0, 0, 0, 0])
