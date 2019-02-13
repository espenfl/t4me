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

import sys
import math

# The complete Fermi-Dirac integral F_j(x) is given by,
# F_j(x)   := (1/\Gamma(j+1)) \int_0^\infty dt (t^j / (\exp(t-x) + 1))
# Note that the Fermi-Dirac integral is sometimes defined without the
# normalisation factor in other texts.


cdef extern from "gsl/gsl_sf_fermi_dirac.h":
    # These routines compute the complete Fermi-Dirac integral with an index
    # of -1. This integral i   # s given by F_{-1}(x) = e^x / (1 + e^x).
    double gsl_sf_fermi_dirac_m1(double x)

    # These routines compute the complete Fermi-Dirac integral with an index
    # of 0. This integral is   # given by F_j(x) = \ln(1 + e^x).
    double gsl_sf_fermi_dirac_0(double x)

    # These routines compute the complete Fermi-Dirac integral with an index
    # of 1, F_1(x) = \int_0^\  # infty dt (t /(\exp(t-x)+1)).
    double gsl_sf_fermi_dirac_1(double x)

    # These routines compute the complete Fermi-Dirac integral with an index
    # of 2, F_2(x) = (1/2) \i  # nt_0^infty dt (t^2 /(\exp(t-x)+1)).
    double gsl_sf_fermi_dirac_2(double x)

    # These routines compute the complete Fermi-Dirac integral with an integer
    # index of j, F_j(x) =   # (1/\Gamma(j+1)) \int_0^\infty dt (t^j
    # /(\exp(t-x)+1)).
    double gsl_sf_fermi_dirac_int(int j, double x)

    # These routines compute the complete Fermi-Dirac integral F_{-1/2}(x).
    double gsl_sf_fermi_dirac_mhalf(double x)

    # These routines compute the complete Fermi-Dirac integral F_{1/2}(x).
    double gsl_sf_fermi_dirac_half(double x)

    # These routines compute the complete Fermi-Dirac integral F_{3/2}(x).
    double gsl_sf_fermi_dirac_3half(double x)

    # These routines compute the Gamma function \Gamma(x), subject to x not
    # being a negative integer or zero. The function is computed using the real
    # Lanczos method. The maximum value of x such that \Gamma(x) is not considered
    # an overflow is given by the macro GSL_SF_GAMMA_XMAX and is 171.0.
    double gsl_sf_gamma(double x)


def fermidiracint_m1(double x):
    print("F_{-1}(x) is not well defined due to the Gamma renormalization is infinite...exiting!")
    sys.exit(1)
    # return gsl_sf_gamma(x+1)*gsl_sf_fermi_dirac_m1(x)


def fermidiracint_0(double x):
    # gamma renorm = 1
    return gsl_sf_fermi_dirac_0(x)


def fermidiracint_1(double x):
    # gamma renorm = 1
    return gsl_sf_fermi_dirac_1(x)


def fermidiracint_2(double x):
    # gamma renorm = 2
    return 2.0 * gsl_sf_fermi_dirac_2(x)


def fermidiracint_int(int j, double x):
    # gamma renorm = j!
    return math.factorial(j) * gsl_sf_fermi_dirac_int(j, x)


def fermidiracint_mhalf(double x):
    # gamma renorm = sqrt(pi)
    return math.sqrt(math.pi) * gsl_sf_fermi_dirac_mhalf(x)


def fermidiracint_half(double x):
    # gamma renorm = sqrt(pi)/2
    return (math.sqrt(math.pi) / 2) * gsl_sf_fermi_dirac_half(x)


def fermidiracint_3half(double x):
    # gamma renorm = 3*sqrt(pi)/4
    return (3 * math.sqrt(math.pi) / 4) * gsl_sf_fermi_dirac_3half(x)


def gamma(double x):
    return gsl_sf_gamma(x)
