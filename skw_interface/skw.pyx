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

import cython
import sys
import logging
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector

cdef extern from "skw_interface.cpp":
    int interpolate_skw_interface(double * energies, int num_bands, double * kpoints, int num_kpoints, int * ksampling, double * lattice, double * positions, int * species, int num_atoms, double star_radius_factor, double * ienergies, double * ivelocities, double * ikpoints, int num_ikpoints, int * iksampling)


@cython.boundscheck(False)
@cython.wraparound(False)
def interpolate(np.ndarray[double, ndim=2, mode="c"] energies not None, np.ndarray[double, ndim=2, mode="c"] kpoints not None, np.ndarray[int, ndim=1, mode="c"] ksampling not None, np.ndarray[double, ndim=2, mode="c"] lattice not None, np.ndarray[double, ndim=2, mode="c"] positions not None, np.ndarray[int, ndim=1, mode="c"] species not None, double star_radius_factor):
    # configure logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)
    # this is slightly bigger than the expected array size, disect later
    # we want the determination to be handled in the C routine,
    # a good estimate is
    cdef int est_num_ikpoints = int(np.power(np.ceil(star_radius_factor * 6.0 *
                                                     np.amax(ksampling) / np.pi), 3.0))
    cdef int num_atoms = species.shape[0]
    cdef int num_kpoints = kpoints.shape[0]
    cdef int num_bands = energies.shape[0]
    # use one dimensional linear arrays here since they are overpadded in size
    # it is easier to disect away the end on linear arrays
    cdef np.ndarray[double, ndim = 1, mode = "c"] ikpoints = np.zeros(est_num_ikpoints * 3, dtype="double")
    cdef np.ndarray[double, ndim = 1, mode = "c"] ienergies = np.zeros(num_bands * est_num_ikpoints, dtype="double")
    cdef np.ndarray[double, ndim = 1, mode = "c"] ivelocities = np.zeros(num_bands * 3 * est_num_ikpoints, dtype="double")
    cdef np.ndarray[int, ndim = 1, mode = "c"] iksampling = np.zeros(3, dtype="intc")
    cdef int info
    info = interpolate_skw_interface(& energies[0, 0], num_bands, & kpoints[0, 0], num_kpoints, & ksampling[0], & lattice[0, 0], & positions[0, 0], & species[0], num_atoms, star_radius_factor, & ienergies[0], & ivelocities[0], & ikpoints[0], est_num_ikpoints, & iksampling[0])
    if info:
        logger.error(
            "SKW: Serious error during the generation of lambdas (LAPACK/BLAS). Exiting.")
        sys.exit(1)
    # disect padding of the arrays and return, true size in iksampling
    cdef int num_ikpoints = np.prod(iksampling)
    return np.reshape(ikpoints[0:num_ikpoints * 3], (num_ikpoints, 3)), np.reshape(ienergies[0:num_ikpoints * num_bands], (num_bands, num_ikpoints)), np.reshape(ivelocities[0:num_ikpoints * 3 * num_bands], (num_bands, 3, num_ikpoints)), iksampling
