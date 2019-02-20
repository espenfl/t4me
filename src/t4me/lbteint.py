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

#!/usr/bin/python
"""Contains the routines to perform the Boltzmann transport integrals."""

# pylint: disable=useless-import-alias, too-many-arguments, invalid-name, too-many-statements

import sys
import math
import logging
import numpy as np
import scipy.integrate

import t4me.scattering as scattering  # pylint: disable=useless-import-alias
import t4me.constants as constants  # pylint: disable=useless-import-alias
import t4me.utils as utils  # pylint: disable=useless-import-alias


def scipy_k_integrals(eta, beta, effmass, e0, i, l, m, method="tplquad"):
    r"""
    Calculates the three dimensional wave vector integrals.

    Uses the SciPy function :func:`scipy.integrate.tplquad`.

    Parameters
    ----------
    eta : float
        The reduced chemical potential
    beta : float
        The :math:`\\beta` factor, :math:`(\\mathrm{k_b}T)^{-1}` in
        eV.
    effmass : ndarray
        | Dimension: (3)

        The effective mass along the three reciprocal unit vectors in
        units of the free electron mass.
    e0 : float
        The energy shift, e.g. :math:`E=\\hbar^2 k^2/2m + E_0`, where
        :math:`E_0` is the energy shift in eV.
    i : int
        The order of the transport tensor.
    l : {0,1,2}
        The first index of the transport tensor
    m : {0,1,2}
        The second index of the transport tensor
    method : {"tplquad"}, optional
        The SciPy three dimensional integration method using
        :func:`scipy.integrate.tplquad`.

    Returns
    -------
    float
        The resulting integral over the wave vectors.

    See Also
    --------
    scipy.integrate.tplquad

    """

    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running scipy_k_integrals.")
    func = {"tplquad": scipy.integrate.tplquad}
    if method not in func:
        logger.error("The supplied method is not recognized. Exiting.")
        sys.exit(1)
    return func[method](
        analytic_k_space_integrand,
        -1.0,
        1.0,
        lambda kx: -1.0,
        lambda kx: 1.0,
        lambda kx, ky: -1.0,
        lambda kx, ky: 1.0,
        args=(eta, beta, effmass, e0, i, l, m),
        epsabs=1e-2)[0]


def scipy_k_integrals_discrete(tr,
                               integrand_type,
                               energies,
                               velocities,
                               scatter,
                               chempot,
                               beta,
                               order,
                               spin_fact,
                               method="trapz"):
    r"""
    Calculates the three dimensional integrals over the k-points for discrete data.

    Uses SciPy integration functions for discrete data.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    energies: ndarray
        Contains the band energies in eV for each k-point.
    velocities: ndarray
        Contains the derivative if `energies` without the `\\hbar` factors for each k-point.
    scatter:
        Contains the relaxation time at each k-point.
    chempot : float
        The chemical potential in eV
    beta : float
        The :math:`\\beta` factor, :math:`(\\mathrm{k_b}T)^{-1}` in eV.
    spin_fact : int
        The spin factor, 1 for non-spin degeneracy and 2 for spin degeneracy.
    method : {"trapz", "simps", "romb"}, optional
        The SciPy three dimensional integration method for the
        :func:`scipy.integrate.trapz`, :func:`scipy.integrate.simps` and
        the :func:`scipy.integrate.romb` functions, respectively. Defaults
        to "trapz".

    Returns
    -------
    integral : float
        The resulting integral over the wave vectors.

    """

    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running scipy_k_integrals for discrete data.")
    func = {
        "trapz": scipy.integrate.trapz,
        "simps": scipy.integrate.simps,
        "romb": scipy.integrate.romb
    }
    if method not in func:
        logger.error("The supplied method is not recognized. Exiting.")
        sys.exit(1)
    # fetch step size in direct coordinates and set Jacobian, which is
    # needed since we perform the integration in direct coordinates
    kx, ky, kz = tr.lattice.fetch_kmesh_step_size(direct=True)
    jacobian = np.linalg.det(tr.lattice.runitcell)
    ksampling = tr.lattice.ksampling
    # now if we want romberg, we need to check for grid samples
    if method == "romb":
        if not utils.is_power_of_two(ksampling[0] - 1):
            logger.error("User requests romberg integration, but "
                         "the samplings in the first direction is not "
                         "2^k - 1. Exiting.")
            sys.exit(1)
        if not utils.is_power_of_two(ksampling[1] - 1):
            logger.error("User requests romberg integration, but "
                         "the samplings in the second direction is not "
                         "2^k - 1. Exiting.")
            sys.exit(1)
        if not utils.is_power_of_two(ksampling[2] - 1):
            logger.error("User requests romberg integration, but "
                         "the samplings in the third direction is not "
                         "2^k - 1. Exiting.")
            sys.exit(1)

    # now prepare the data that is to be integrated
    if integrand_type == "normal":
        integrand = concatenate_integrand(energies, velocities, scatter,
                                          spin_fact, chempot, beta, order)
        integrand_shaped = integrand.reshape(3, 3, ksampling[0], ksampling[1],
                                             ksampling[2])
    else:
        logger.error(
            "The supplied integrand_type: %s is not supported. Exiting.",
            integrand_type)
        sys.exit(1)

    integral = func[method](
        func[method](
            func[method](integrand_shaped, dx=kz, axis=4), dx=ky, axis=3),
        dx=kx,
        axis=2)

    # add Jacobian (we integrate in direct coordinates)
    integral = jacobian * integral

    return integral


def scipy_k_integrals_discrete2(tr,
                                energies,
                                velocities,
                                scatter,
                                chempot,
                                beta,
                                spin_fact,
                                order,
                                method="trapz"):
    r"""
    Calculates the three dimensional integrals over the k-points for discrete data.

    Uses integration functions for discrete data.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    chempot : float
        The chemical potential in eV
    beta : float
        The :math:`\\beta` factor, :math:`(\\mathrm{k_b}T)^{-1}` in eV.
    spin_fact : int
        The spin factor, 1 for non-spin degeneracy and 2 for spin degeneracy.
    kx, ky, kz : float, float, float
        The spacing in inverse AA between the points along each direction.
    order : float
        The order of the energy minus chemical potential term in the
        denominator.
    method : {"trapz", "simps", "romb"}, optional
        The SciPy three dimensional integration method for the
        :func:`scipy.integrate.trapz`, :func:`scipy.integrate.simps` and
        the :func:`scipy.integrate.romb` functions, respectively. Defaults
        to "trapz".

    """

    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running scipy_k_integrals for discrete data.")
    func = {
        "trapz": scipy.integrate.trapz,
        "simps": scipy.integrate.simps,
        "romb": scipy.integrate.romb
    }
    if method is None:
        method = tr.param.transport_integration_method
    if method not in func:
        logger.error("The supplied method is not recognized. Exiting.")
        sys.exit(1)

    # fetch step size in direct coordinates and set Jacobian, which is
    # needed since we perform the integration in direct coordinates
    kx, ky, kz = tr.lattice.fetch_kmesh_step_size(direct=True)
    jacobian = np.linalg.det(tr.lattice.runitcell)
    ksampling = tr.lattice.ksampling

    # set integrand
    integrand = concatenate_integrand_band(energies, velocities, scatter,
                                           spin_fact, chempot, beta, order)

    # reshape integrand
    integrand_shaped = integrand.reshape(3, 3, ksampling[0], ksampling[1],
                                         ksampling[2])

    # now if we want romberg, we need to check for grid samples
    if method == "romb":
        logger.debug("Running SciPy Romberg integration for discrete " "data.")
        if not utils.is_power_of_two(ksampling[0] - 1):
            logger.error("User requests romberg integration, but "
                         "the samplings in the first direction is not "
                         "2^k - 1. Exiting.")
            sys.exit(1)
        if not utils.is_power_of_two(ksampling[1] - 1):
            logger.error("User requests romberg integration, but "
                         "the samplings in the second direction is not "
                         "2^k - 1. Exiting.")
            sys.exit(1)
        if not utils.is_power_of_two(ksampling[2] - 1):
            logger.error("User requests romberg integration, but "
                         "the samplings in the third direction is not "
                         "2^k - 1. Exiting.")
            sys.exit(1)

    elif method == "trapz":
        logger.debug("Running SciPy trapeziodal integration for "
                     "discrete data.")
    elif method == "simps":
        logger.debug("Running SciPy Simpson integration for discrete " "data.")

    integral = func[method](
        func[method](
            func[method](integrand_shaped, dx=kz, axis=4), dx=ky, axis=3),
        dx=kx,
        axis=2)

    # add Jacobian (we integrate in direct coordinates)
    integral = jacobian * integral

    return integral


def scipy_e_integrals(transport,
                      integrand,
                      e_min,
                      e_max,
                      w0,
                      eta,
                      beta,
                      energy_trans,
                      effmass,
                      order,
                      spin_fact,
                      method="quad"):
    r"""
    Calculates the one dimensional energy integrals.

    Uses the SciPy function :func:`scipy.integrate.quad`.

    Parameters
    ----------
    transport : object
        A `Transport()` object
    integrand : {"normal","hall","dos"}
        Selects the type of integrand to be used. "normal" selects
        :func:`integrandpar`. "hall" selects :func:`integrandpart2`,
        "dos" selects :func:`integrandpardos`.
    e_min : float
        The lower integration limit in eV.
    e_max : float
        The higher integration limit in eV.
    w0 : ndarray
        | Dimension: (12)

        Contains the scattering rate prefactor (inverse of
        relaxation time) for the different scattering
        mechanisms in units of inverse fs.
    eta : float
        The reduced chemical potential.
    beta : float
        The :math:`\\beta` factor, :math:`(\\mathrm{k_b}T)^{-1}`
        in eV.
    effmass : ndarray
        | Dimension: (3)

        The effective mass along the three reciprocal unit
        vectors in units of the free electron mass.
    e0 : float
        The energy shift, e.g. :math:`E=\\hbar^2 k^2/2m + E_0`,
        where :math:`E_0` is the energy shift in eV.
    i : int
        The order of the transport tensor.
    l : {0,1,2}
        The first index of the transport tensor.
    m : {0,1,2}
        The second index of the transport tensor.
    spin_fact : int
        The spin degeneracy factor. 1 for non-spin degeneracy,
        2 for spin degeneracy.
    method : {"quad"}, optional
        The SciPy three dimensional integration method using
        :func:`scipy.integrate.quad`.

    Returns
    -------
    float
        The resulting integral over energy.

    """
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running scipy_e_integrals.")
    func = {"quad": scipy.integrate.quad}
    if method not in func:
        logger.error("The supplied method is not recognized. Exiting.")
        sys.exit(1)
    integrands = {
        "normal": integrandpar,
        "hall": integrandpart2,
        "dos": integrandpardos
    }
    if integrand not in integrands:
        logger.error("The supplied integrand is not valid. Exiting.")
        sys.exit(1)
    integral = func[method](
        integrands[integrand],
        e_min,
        e_max,
        args=(transport, w0, eta, beta, energy_trans, effmass, order))[0]

    return spin_fact * integral


def fermiintclosed(order, eta, spin_fact):  # noqa: MC0001
    r"""
    Returns the value of the closed expressions for the Fermi integrals.

    Parameters
    ----------
    order : integer
        The Fermi integral order (two times :math:`r`
        to avoid half integers).
    eta : float
        The chemical potential given in reduced
        form (:math:`\\mu/kT`, dimensionless).
    spin_fact : int
        The spin degeneracy. 1 for non-spin degeneracy
        and 2 for spin degeneracy.

    Returns
    -------
    float
         The value of the Fermi integral.

    Notes
    -----
    Utilizes the GSL :cite:`gsl` library and a few inlined
    function from the literature. Consult Ref.
    :cite:`halen_1985_joap_assatfioo` and Ref.
    :cite:`halen_1986_joap_essatfioo`
    as a suplement. The Gamma factor renormalization which
    is included in the GSL returns should be removed and
    this is done in the interface Cython file.

    .. math:: F_i=\\int_0^{\\infty}\\epsilon^i d\\epsilon /
              (1+\\exp(\\epsilon-\\eta))

    .. warning:: in order to avoid half numbers, order
                 should be given as two times the actual
                 order of the integral.

    .. rubric:: References

    .. bibliography:: references.bib
        :style: unsrt
        :filter: docname in docnames

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running fermiintclosed.")

    # lazy import og gsl
    import t4me.gsl as gsl  # pylint: disable=import-error, no-name-in-module

    # check for bogus r values
    if order < -2:
        logger.error(
            "Bogus r passed to fermiintclosed, r < -2, r=%s. Exiting.", order)
        sys.exit(1)
    if order % 2 != 0 and order > 8:
        logger.error(
            "Bogus r passed to fermiintclosed, r > 4 and not integer, r=%s. Exiting.",
            order)
        sys.exit(1)
    if order == -2:
        integral = gsl.fermidiracint_m1(eta)
    elif order == -1:
        integral = gsl.fermidiracint_mhalf(eta)
    elif order == 1:
        integral = gsl.fermidiracint_half(eta)
    elif order == 3:
        integral = gsl.fermidiracint_3half(eta)
    elif order == 5:
        # after JAP 57, 5271, 1985 (also see erratum)
        gamma = 15.0 * np.sqrt(constants.pi) / 8.0
        a1 = np.array(
            [1.0, 0.088392, 0.021407, 0.007917, 0.003723, 0.001716, 0.000451])
        a2 = np.array(
            [0.085972, 1.23738, 1.07293, 0.362030, 38.7579, -750.718, 4378.70])
        a3 = np.array([
            0.927560, 0.866971, 0.383690, 0.098863, 0.017398, 0.000418,
            -0.000067
        ])
        if eta <= 0.0:
            f = 0.0
            for index, value in np.ndenumerate(a1):
                f = f + np.power(-1.0, index[0] + 2) * \
                    value * np.exp(eta * (index[0] + 1))

        elif 0.0 < eta <= 4.0:
            f = 0.0
            for index, value in np.ndenumerate(a3):
                f = f + value * np.power(eta, index[0])
        else:
            f = 0.0
            for index, value in np.ndenumerate(a2):
                f = f + value / np.power(eta, 2 * (index[0]))
            f = f * np.power(eta, 3.5)
        integral = f * gamma

    elif order == 7:
        # after JAP 57, 5271, 1985 (also see erratum)
        gamma = 40320 * math.sqrt(math.pi) / 6144
        a1 = np.array(
            [1.0, 0.044203, 0.007157, 0.001976, 0.000719, 0.000317, 0.000106])
        a2 = np.array([
            0.019105, 0.494958, 2.13722, -0.503902, -6.99243, 96.6031, -426.046
        ])
        a3 = np.array([
            0.961478, 0.927751, 0.432494, 0.129617, 0.023308, 0.004067,
            -0.000051
        ])
        if eta <= 0.0:
            f = 0.0
            for index, value in np.ndenumerate(a1):
                f = f + np.power(-1.0, index[0] + 2) * \
                    value * np.exp(eta * (index[0] + 1))

        elif 0.0 < eta <= 4.0:
            f = 0.0
            for index, value in np.ndenumerate(a3):
                f = f + value * np.power(eta, index[0])
        else:
            f = 0.0
            for index, value in np.ndenumerate(a2):
                f = f + value / np.power(eta, 2 * (index[0]))
            f = f * np.power(eta, 3.5)
        integral = f * gamma

    elif order == 0:
        integral = gsl.fermidiracint_0(eta)
    elif order == 2:
        integral = gsl.fermidiracint_1(eta)
    elif order == 4:
        integral = gsl.fermidiracint_2(eta)
    elif order > 4 or order < -2:
        if order % 2 == 0:
            integral = gsl.fermidiracint_int(order / 2, eta)
    else:
        logger.error("No Fermi-Dirac integration routines "
                     "selected. Exiting.")
        sys.exit(1)

    return spin_fact * integral


def integrandpardos(eps, transport, w0, eta, beta, energy_trans, effmass, i):  # pylint: disable=unused-argument
    r"""
    The integrand for the density of states integral over energy.

    Parameters
    ----------
    eps : float
        The reduced carrier energy
    transport : object
        A `Transport()` object
    w0 : ndarray
        | Dimension: (12)

        Contains the scattering rate prefactor for the different
        scattering mechanisms in units of inverse fs. Not
        used in this routine, but it needs the dummy from the
        call argument.
    eta : float
        The reduced chemical potential
    beta : float
        The :math:`\\beta` factor in eV. Not
        used in this routine, but it needs the dummy from the
        call argument.
    energy_trans : ndarray
        | Dimension: (12)

        Contains the energy transitions (that is added to the
        energy in :math:`\\tau=\\tau_0E^{r-1/2}`, typically,
        :math:`E=E+\\hbar \\omega`, where :math:`\\hbar \\omega`
        is the size of the energy transition. Set it to zero
        for the non-relevant scattering mechanisms. Not
        used in this routine, but it needs the dummy from the
        call argument.
    effmass : float
        The effective mass in units of the electron mass. Not
        used in this routine, but it needs the dummy from the
        call argument.
    i : int
        The order of the transport integral to be evaluated.
        Not used in this routine, but it needs the dummy from the
        call argument.

    Returns
    -------
    float
       The integrand value for the density
       of states.

    Notes
    -----
    Calculates the density of states integrand

    .. math:: \\frac{\\epsilon^{1/2}}{1+\\exp{(\\epsilon-\\eta)}}

    """

    return math.sqrt(eps) / (1 + math.exp(eps - eta))


def integrandpar(eps, transport, w0, eta, beta, energy_trans, effmass, i):  # pylint: disable=unused-argument
    r"""
    Returns the integrand used in the analytic energy integration of the transport coefficients in :func:`integrate_e`

    Parameters
    ----------
    eps : float
        The reduced carrier energy.
    transport : object
        A `Transport()` object.
    w0 : ndarray
        | Dimension: (12)

        Contains the scattering rate prefactor for the different
        scattering mechanisms in units of inverse fs.
    eta : float
        The reduced chemical potential.
    beta : float
        The :math:`\\beta` factor in eV.
    energy_trans : ndarray
        | Dimension: (12)

        Contains the energy transitions (that is added to the
        energy in :math:`\\tau=\\tau_0E^{r-1/2}`, typically,
        :math:`E=E+\\hbar \\omega`, where :math:`\\hbar \\omega`
        is the size of the energy transition. Set it to zero
        for the non-relevant scattering mechanisms.
    effmass : float
       The effective mass in units of the electron mass.
    i : int
        The order of the transport integral to be evaluated.

    Returns
    -------
    float
        The integrand value.

    Notes
    -----

    The total scattering is calculated based on the well
    known scattering models for parabolic energy dispersions
    :math:`\\tau=\\tau_0\\epsilon^{r-1/2}`,
    where :math:`r` is the scattering factor.

    """

    return scattering.combined_scattering(
        transport, np.absolute(eps), w0, energy_trans) * \
        pow(eps, 1.5) * pow(eps - eta, i) / \
        (1 + math.cosh(eps - eta))


def integrandpart2(eps, transport, w0, eta, beta, energy_trans, effmass, i):  # pylint: disable=unused-argument
    r"""
    Returns the integrand used in the analytic energy integration of the transport distribution function with a quadratic :math:`\\tau` term

    Parameters
    ----------
    eps : float
        The reduced carrier energy.
    transport : object
        A `Transport()` object.
    w0 : ndarray
        | Dimension: (12)

        Contains the scattering rate prefactor for the different
        scattering mechanisms in units of inverse fs.
    eta : float
        The reduced chemical potential.
    beta : float
        The :math:`\\beta` factor in eV.
    energy_trans : ndarray
        | Dimension: (12)

        Contains the energy transitions (that is added to the
        energy in :math:`\\tau=\\tau_0E^{r-1/2}`,
        typically, :math:`E=E+\\hbar \\omega`,
        where :math:`\\hbar \\omega` is the size of the
        energy transition. Set it to zero for the
        non-relevant scattering mechanisms.
    effmass : float
        The effective mass in units of the electron mass.
    i : int
        The order of the transport integral to be evaluated.

    Returns
    -------
    float
        The integrand value.

    Notes
    -----
    The total scattering is calculated based on the well
    known scattering models for parabolic energy dispersions
    :math:`\\tau=\\tau_0\\epsilon^{r-1/2}`,
    where :math:`r` is the scattering factor.
    Difference from :func:`integrandpar`: here tau^2
    is used in the integrand (for the
    calculation of the Hall factor).

    """

    return np.power(scattering.combined_scattering(
        transport, np.absolute(eps), w0, energy_trans), 2.0) * \
        pow(eps, 1.5) * pow(eps - eta, i) / \
        (1 + math.cosh(eps - eta))


def analytic_k_space_integrand(kz, ky, kx, eta, beta, effmass, e0, i, l, m):
    r"""
    Returns the integrand for the anlytic reciprocal space integration of the transport tensor.

    Parameters
    ----------
    kz : float
        The :math:`k_z` in cartesian coordinates.
    ky : float
        The :math:`k_y` in cartesian coordinates.
    kx : float
        The :math:`k_x` in cartesian coordinates.
    eta : float
        The reduced chemical potential.
    beta : float
        The :math:`\\beta` factor, :math:`\\beta=(k_bT)^{-1}`
        in eV.
    effmass : float
        The effective mass in units of the free electron mass.
    e0 : float
        The energy shift, e.g. :math:`E=\\hbar^2 k^2/2m + E_0`,
        where :math:`E_0` is the energy shift in eV.
    i : int
        The order of the transport tensor.
    l : {0,1,2}
        The first index of the transport tensor.
    m : {0,1,2}
        The second index of the transport tensor.

    Returns
    -------
    float
        The integrand value.

    """
    energy = analytic_k_space_energy(kx, ky, kz, effmass, e0)
    return analytic_k_space_velocity(kx, ky, kz, effmass, l) * \
        analytic_k_space_velocity(kx, ky, kz, effmass, m) * \
        pow(energy - eta / beta, i) / \
        (1 + math.cosh(energy * beta - eta))


def analytic_k_space_energy(kx, ky, kz, effmass, e_shift):
    """
    Returns the parabolic energy dispersion.

    Parameters
    ----------
    transport : object
        A `Transport()` object
    kx : float
        The :math:`k_x` in cartesian coordinates.
    ky : float
        The :math:`k_y` in cartesian coordinates.
    kz : float
        The :math:`k_z` in cartesian coordinates.
    effmass : ndarray
        | Dimension: (3)

        The effective mass along :math:`k_x`,
        :math:`k_y` and :math:`k_z`, respectively.

    Returns
    -------
    float
        The energy value in eV.

    .. warning:: This routine only accepts the diagonal
                 elements of the effective mass tensor

    """

    # Prefactor of hbar^2/(2m) not included
    # (beware for non-parabolic expressions)
    # Need to shift the energy in order to
    # avoid E=0 for the summation of the 1/E for
    # the scattering. Should not be a problem
    return constants.bandunit * (kx * kx / effmass[0] + ky * ky / effmass[1] +
                                 kz * kz / effmass[2]) + e_shift


def analytic_k_space_velocity(kx, ky, kz, effmass, i):
    r"""
    Returns the parabolic velocity dispersion.

    Parameters
    ----------
    kx : float
        The :math:`k_x` in cartesian coordinates.
    ky : float
        The :math:`k_y` in cartesian coordinates.
    kz : float
        The :math:`k_z` in cartesian coordinates.
    effmass : ndarray
        | Dimension: (3)

        The effective mass along :math:`k_x`,
        :math:`k_y` and :math:`k_z`, respectively.
    i : {0,1,2}
        The direction to evaluate the velocity
        (0 is along :math:`k_x` etc.).

    Returns
    -------
    float
        The velocity in eVAA.

    .. warning:: This routine only accepts the diagonal
                 elements of the effective mass tensor.
                 The :math:`\\hbar/m_e` factor is not
                 returned and need to be introduced
                 externally.

    """

    # Prefactor of hbar/m not included
    # (beware for non-parabolic expressions)
    if i == 0:
        return 2.0 * constants.bandunit * kx / effmass[0]
    if i == 1:
        return 2.0 * constants.bandunit * ky / effmass[1]
    if i == 2:
        return 2.0 * constants.bandunit * kz / effmass[2]

    return None


def concatenate_integrand(energies, velocities, scatter, spin_fact, chempot,
                          beta, order):
    """Concatenates the integrand in the Boltzmann transport integral."""
    denom = 1.0 + np.cosh((energies - chempot) * beta)
    integrand = spin_fact * np.power(energies[
        np.newaxis, np.newaxis, :] - chempot, order) * \
        velocities[:, np.newaxis, :] * velocities[np.newaxis, :, :] * \
        scatter[np.newaxis, np.newaxis, :] / \
        denom[np.newaxis, np.newaxis, :]
    return integrand


def concatenate_integrand_band(energies, velocities, tau, spin_fact, chempot,
                               beta, order):
    """Concatenates the integrand in the Boltzmann transport integral and sums the bands."""
    denom = 1.0 + np.cosh((energies - chempot) * beta)
    integrand = spin_fact[:, np.newaxis, np.newaxis, np.newaxis] * \
        np.power(energies[:, np.newaxis, np.newaxis, :] - chempot, order) * \
        velocities[:, :, np.newaxis, :] * velocities[:, np.newaxis, :, :] * \
        tau[:, np.newaxis, np.newaxis, :] / \
        denom[:, np.newaxis, np.newaxis, :]
    # sum bands
    integrand = np.sum(integrand, axis=0)
    return integrand
