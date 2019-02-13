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

cdef extern from "wildmagic_interface.hpp":
    void wildmagic_execute_interpolation(int * num_points, double * domainx, double * domainy, double * domainz, double * data, double * ix, double * iy, double * iz, int ip, int num_bands, double * idata, int itype)

    void wildmagic_gradient_execute_interpolation(int * num_points, double * domainx, double * domainy, double * domainz, double * data, double * ix, double * iy, double * iz, int ip, int num_bands, double * idata, double * igradx, double * igrady, double * igradz, int itype)


@cython.boundscheck(False)
@cython.wraparound(False)
def trilinear_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None):
    logging.debug("Wildmagic: Running trilinear interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_execute_interpolation(& num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], 0)


def trilinear_gradient_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None, np.ndarray[double, ndim=2, mode="c"] igradx not None, np.ndarray[double, ndim=2, mode="c"] igrady not None, np.ndarray[double, ndim=2, mode="c"] igradz not None):
    logging.debug(
        "Wildmagic: Running trilinear (also fetching gradients) interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_gradient_execute_interpolation( & num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], & igradx[0, 0], & igrady[0, 0], & igradz[0, 0], 0)


def tricubic_exact_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None):
    logging.debug("Wildmagic: Running tricubic (exact) interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_execute_interpolation(& num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], 1)


def tricubic_exact_gradient_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None, np.ndarray[double, ndim=2, mode="c"] igradx not None, np.ndarray[double, ndim=2, mode="c"] igrady not None, np.ndarray[double, ndim=2, mode="c"] igradz not None):
    logging.debug(
        "Wildmagic: Running tricubic (exact, also fetching gradients) interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_gradient_execute_interpolation( & num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], & igradx[0, 0], & igrady[0, 0], & igradz[0, 0], 1)


def tricubic_bspline_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None):
    logging.debug("Wildmagic: Running tricubic (bspline) interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_execute_interpolation(& num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], 2)


def tricubic_bspline_gradient_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None, np.ndarray[double, ndim=2, mode="c"] igradx not None, np.ndarray[double, ndim=2, mode="c"] igrady not None, np.ndarray[double, ndim=2, mode="c"] igradz not None):
    logging.debug(
        "Wildmagic: Running tricubic (bspline, also fetching gradients) interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_gradient_execute_interpolation( & num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], & igradx[0, 0], & igrady[0, 0], & igradz[0, 0], 2)


def akima_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None):
    logging.debug("Wildmagic: Running uniform Akima interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_execute_interpolation(& num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], 3)


def akima_gradient_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data, np.ndarray[double, ndim=1, mode="c"] ix not None, np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None, np.ndarray[double, ndim=2, mode="c"] igradx not None, np.ndarray[double, ndim=2, mode="c"] igrady not None, np.ndarray[double, ndim=2, mode="c"] igradz not None):
    logging.debug(
        "Wildmagic: Running uniform Akima (also fetching gradients) interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    wildmagic_gradient_execute_interpolation( & num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], & igradx[0, 0], & igrady[0, 0], & igradz[0, 0], 3)
