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

import sys
import cython
import logging
import numpy as np
cimport numpy as np
from libcpp cimport bool
import t4me.constants as constants

cdef extern from "cubature_wildmagic_interface.hpp":
    void cubature_wildmagic_execute_integration(int * num_points,
                                                double * domainx,
                                                double * domainy,
                                                double * domainz,
                                                double * data,
                                                double & value,
                                                double & error,
                                                int itype)

    void calc_transport_tensors_cubature(int * num_points,
                                         double * domainx,
                                         double * domainy,
                                         double * domainz,
                                         double * energy,
                                         int num_bands,
                                         double * velocity,
                                         double * tau,
                                         double * energy_samples,
                                         int * spin_degen,
                                         int num_scatterings,
                                         int num_energy_samples,
                                         double * tau_energy_trans,
                                         double * chempot,
                                         int num_chempot,
                                         double * temperatures,
                                         int num_temps,
                                         int num_kpoints_ibz,
                                         int * k_sampling,
                                         double * rec_basis,
                                         int tensor_type,
                                         int interpol_type,
                                         int integral_type,
                                         int max_it,
                                         double abs_err,
                                         double rel_err,
                                         int iso,
                                         int info,
                                         double * cond,
                                         double * seebeck,
                                         double * lorenz,
                                         int gen_velocities,
                                         int scattering_type)

    void calc_dos_cubature(int * num_points,
                           double * domainx,
                           double * domainy,
                           double * domainz,
                           double * energies,
                           double * energy_samples,
                           int * spin_degen,
                           int num_energy_samples,
                           int num_bands,
                           double sigma,
                           int interpol_method,
                           int smear_then_interpolate,
                           double volume, double volume_bz,
                           unsigned max_it,
                           double abs_err,
                           double rel_err, int integral_type,
                           int info,
                           double * dos, double * int_dos)


@cython.boundscheck(False)
@cython.wraparound(False)
def cubature_trilinear_execute_interface(
        np.ndarray[int, ndim=1, mode="c"] num_points not None,
        np.ndarray[double, ndim=1, mode="c"] domainx not None,
        np.ndarray[double, ndim=1, mode="c"] domainy not None,
        np.ndarray[double, ndim=1, mode="c"] domainz not None,
        np.ndarray[double, ndim=1, mode="c"] data):

    logging.debug("Wildmagic: Running trilinear interpolation \
     and cubature integration")
    cdef double value = 0.0
    cdef double error = 0.0
    cubature_wildmagic_execute_integration( & num_points[0],
                                           & domainx[0],
                                           & domainy[0],
                                           & domainz[0],
                                           & data[0],
                                           value, error, 0)

    return value, error


def cubature_tricubic_exact_execute_interface(
        np.ndarray[int, ndim=1, mode="c"] num_points not None,
        np.ndarray[double, ndim=1, mode="c"] domainx not None,
        np.ndarray[double, ndim=1, mode="c"] domainy not None,
        np.ndarray[double, ndim=1, mode="c"] domainz not None,
        np.ndarray[double, ndim=1, mode="c"] data):

    logging.debug("Wildmagic: Running tricubic (exact) \
     interpolation and cubature integration")
    cdef double value = 0.0
    cdef double error = 0.0
    cubature_wildmagic_execute_integration( & num_points[0],
                                           & domainx[0],
                                           & domainy[0],
                                           & domainz[0],
                                           & data[0],
                                           value, error, 1)

    return value, error


def cubature_tricubic_bspline_execute_interface(
        np.ndarray[int, ndim=1, mode="c"] num_points not None,
        np.ndarray[double, ndim=1, mode="c"] domainx not None,
        np.ndarray[double, ndim=1, mode="c"] domainy not None,
        np.ndarray[double, ndim=1, mode="c"] domainz not None,
        np.ndarray[double, ndim=1, mode="c"] data):

    logging.debug("Wildmagic: Running tricubic (bspline) \
     interpolation and cubature integration")
    cdef double value = 0.0
    cdef double error = 0.0
    cubature_wildmagic_execute_integration( & num_points[0],
                                           & domainx[0],
                                           & domainy[0],
                                           & domainz[0],
                                           & data[0],
                                           value, error, 2)

    return value, error


def cubature_akima_execute_interface(
        np.ndarray[int, ndim=1, mode="c"] num_points not None,
        np.ndarray[double, ndim=1, mode="c"] domainx not None,
        np.ndarray[double, ndim=1, mode="c"] domainy not None,
        np.ndarray[double, ndim=1, mode="c"] domainz not None,
        np.ndarray[double, ndim=1, mode="c"] data):

    logging.debug("Wildmagic: Running uniform Akima \
     interpolation and cubature integration")
    cdef double value = 0.0
    cdef double error = 0.0
    cubature_wildmagic_execute_integration( & num_points[0],
                                           & domainx[0],
                                           & domainy[0],
                                           & domainz[0],
                                           & data[0],
                                           value, error, 3)

    return value, error


def calc_density_of_states_interface(
        np.ndarray[int, ndim=1, mode="c"] num_points not None,
        np.ndarray[double, ndim=1, mode="c"] domainx not None,
        np.ndarray[double, ndim=1, mode="c"] domainy not None,
        np.ndarray[double, ndim=1, mode="c"] domainz not None,
        np.ndarray[double, ndim=2, mode="c"] energies not None,
        np.ndarray[double, ndim=1, mode="c"] energy_samples not None,
        np.ndarray[int, ndim=1, mode="c"] spin_degen not None,
        num_energy_samples, num_bands,
        sigma,
        interpol_method, smear_then_interpolate,
        volume, volume_bz,
        max_it, abs_err, rel_err, h,
        nfo,
        np.ndarray[double, ndim=2, mode="c"] dos not None,
        np.ndarray[double, ndim=2, mode="c"] int_dos not None):

    logging.debug("Cubature: Running density of states interface")
    cdef int integral_type = h
    cdef int info = nfo
    calc_dos_cubature( & num_points[0],
                      & domainx[0], & domainy[0], & domainz[0],
                      & energies[0, 0], & energy_samples[0],
                      & spin_degen[0],
                      num_energy_samples, num_bands, sigma,
                      interpol_method, smear_then_interpolate,
                      volume, volume_bz,
                      max_it, abs_err, rel_err, integral_type,
                      info,
                      & dos[0, 0], & int_dos[0, 0])


def calc_transport_tensors_cubature_interface(
        np.ndarray[int, ndim=1, mode="c"] num_points not None,
        np.ndarray[double, ndim=1, mode="c"] domainx not None,
        np.ndarray[double, ndim=1, mode="c"] domainy not None,
        np.ndarray[double, ndim=1, mode="c"] domainz not None,
        np.ndarray[double, ndim=2, mode="c"] energy not None,
        np.ndarray[double, ndim=3, mode="c"] velocities not None,
        np.ndarray[double, ndim=3, mode="c"] tau not None,
        np.ndarray[double, ndim=1, mode="c"] energy_samples not None,
        np.ndarray[int, ndim=1, mode="c"] spin_degen not None,
        np.ndarray[double, ndim=2, mode="c"] energy_trans not None,
        np.ndarray[double, ndim=1, mode="c"] chempot not None,
        np.ndarray[double, ndim=1, mode="c"] temperatures not None,
        num_kpoints_ibz,
        np.ndarray[int, ndim=1, mode="c"] k_sampling not None,
        np.ndarray[double, ndim=2, mode="c"] rec_basis not None,
        tensor_type, interpol_method, interpol_type,
        max_it, abs_err, rel_err, h, iso,
        nfo,
        np.ndarray[double, ndim=4, mode="c"] sigma not None,
        np.ndarray[double, ndim=4, mode="c"] seebeck not None,
        np.ndarray[double, ndim=4, mode="c"] lorenz not None,
        gen_velocities, scattering_type):

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)
    logger.debug("Wildmagic: Running transport tensor interface")
    cdef int interp_type = 0
    if interpol_method != "wildmagic":
        logger.error(
            "The supplied transport_interpol_method is not supported. Currently only 'wildmagic' is supported for the cubature integration. Exiting.")
        sys.exit(1)
    if interpol_type not in constants.wildmagic_methods:
        logger.error("The specified transport_interpol_type is not recognized by the Wildmagic method, please consult the Wildmagic documentation for valid flags. Possible flags for itype_sub are: " +
                     ', '.join(map(str, constants.wildmagic_methods)) + ". Exiting.")
        sys.exit(1)

    if interpol_type == "trilinear":
        interp_type = 0
    elif interpol_type == "tricubic_exact":
        interp_type = 1
    elif interpol_type == "tricubic_bspline":
        interp_type = 2
    elif interpol_type == "akima":
        interp_type = 3

    cdef int num_scatterings = energy_trans.shape[1]
    cdef int num_energy_samples = energy_samples.shape[0]
    cdef int num_bands = energy.shape[0]
    cdef int num_chempot = chempot.shape[0]
    cdef int num_temps = temperatures.shape[0]
    cdef int integral_type = h
    cdef int info = nfo
    calc_transport_tensors_cubature( & num_points[0],
                                    & domainx[0], & domainy[0],
                                    & domainz[0],
                                    & energy[0, 0], num_bands,
                                    & velocities[0, 0, 0],
                                    & tau[0, 0, 0],
                                    & energy_samples[0],
                                    & spin_degen[0],
                                    num_scatterings,
                                    num_energy_samples,
                                    & energy_trans[0, 0],
                                    & chempot[0], num_chempot,
                                    & temperatures[0], num_temps,
                                    num_kpoints_ibz,
                                    & k_sampling[0], & rec_basis[0, 0],
                                    tensor_type,
                                    interp_type, integral_type,
                                    max_it, abs_err, rel_err, iso,
                                    info,
                                    & sigma[0, 0, 0, 0],
                                    & seebeck[0, 0, 0, 0],
                                    & lorenz[0, 0, 0, 0],
                                    gen_velocities, scattering_type)
