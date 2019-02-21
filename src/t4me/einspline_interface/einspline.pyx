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
import numpy as np
cimport numpy as np
import logging

cdef extern from "einspline_interface.hpp":
    void einspline_execute_uniform(int * num_points, double * boundaryx,
                                   double * boundaryy, double * boundaryz, double * data, double * ix,
                                   double * iy, double * iz, int ip, int num_bands, double * idata,
                                   double * igradx, double * igrady, double * igradz, int grad, char * )
    void einspline_execute_nonuniform(int * num_points, double * gridx, double * gridy, double * gridz, double * data, double * ix, double * iy, double * iz, int ip, int num_bands, double * idata, double * igradx, double * igrady, double * igradz, int grad, char * )


def einspline_execute_interface(np.ndarray[int, ndim=1, mode="c"] num_points not None, np.ndarray[double, ndim=1, mode="c"] domainx not None, np.ndarray[double, ndim=1, mode="c"] domainy not None, np.ndarray[double, ndim=1, mode="c"] domainz not None, np.ndarray[double, ndim=2, mode="c"] data not None, np.ndarray[double, ndim=1, mode="c"] ix not None,  np.ndarray[double, ndim=1, mode="c"] iy not None, np.ndarray[double, ndim=1, mode="c"] iz not None, np.ndarray[double, ndim=2, mode="c"] idata not None, np.ndarray[double, ndim=2, mode="c"] igradx not None, np.ndarray[double, ndim=2, mode="c"] igrady not None, np.ndarray[double, ndim=2, mode="c"] igradz not None, int grad, boundary):
    logging.debug("Running Einspline interpolation")
    cdef int ip = np.intc(ix.shape[0])
    cdef int num_bands = np.intc(data.shape[0])
    # passing strings is a bit iffy, particularly for Py 2 and 3 support
    # make sure we have uppercase and then byte strings
    boundary = boundary.upper()
    bc_bytes = str.encode(boundary)
    cdef char * bc = bc_bytes
    uniform = True
    if uniform:
        logging.debug("Running uniform Einspline interpolation")
        einspline_execute_uniform( & num_points[0], & domainx[0], & domainy[0],
                                  & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0],
                                  ip, num_bands, & idata[0, 0], &
                                  igradx[0, 0], & igrady[0, 0], & igradz[0, 0], grad, &
                                  bc[0])
    else:
        logging.debug("Running non-uniform Einspline interpolation")
        einspline_execute_nonuniform( & num_points[0], & domainx[0], & domainy[0], & domainz[0], & data[0, 0], & ix[0], & iy[0], & iz[0], ip, num_bands, & idata[0, 0], & igradx[0, 0], & igrady[0, 0], & igradz[0, 0], grad, & bc[0])
