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
"""Contains routines that sets up and selects the integral functions to be called."""

# pylint: disable=useless-import-alias, too-many-arguments, invalid-name, too-many-statements, too-many-lines, no-name-in-module

import sys
import math
import logging
import numpy as np

import t4me.lbteint as lbteint
import t4me.constants as constants
import t4me.scattering as scattering
import t4me.utils as utils
from t4me.bandstructure import parabolic_effective_mass


def parabolice(tr, eta, temperature, bs, tau0, method):
    r"""
    A wrapper for all the parabolic Fermi integrals.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    eta : float
        Contains the reduced chemical potential.
    temperature : float
        Contains the temperature in K.
    bs : object
        `Bandstructure()` object containing the band structure.
    method : {"numeric", "closed"}
        How to evaluate the Fermi integrals.

        | "numeric": solve the Fermi integrals using
        | numerical integration
        | "closed": solve the closed Fermi integrals
        | (exact analytic expressions)

    Returns
    -------
    tupple: ndarray, ndarray, ndarray, ndarray, ndarray, nadarray
        | Dimension: (3,3), (3,3), (3,3), (3,3), (3,3), (3,3)

        The electrical conductivity, Seebeck coefficient, Lorenz number,
        Hall coefficient (big R, where the small Hall factor is
        divided by the charge carrier concentration) and charge
        carrier concentration in units of

        .. math:: S/m, \\mu \mathrm{V}/\mathrm{K},
                  10^{-8}\mathrm{V}^2/\mathrm{K}^2,
                  \mathrm{cm}^3/\mathrm{C},
                  10^{21} \mathrm{cm}^{-3}`,

        respectively.

    Notes
    -----
    All integrals in this function is evaluated over energy
    and not k-points (one of the most usual
    procedures for solving the Boltzmann transport integrals)

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running parabolice.")
    func = {"numeric": parabolic_numeric}
    try:
        func[method]
    except KeyError:
        logger.error(
            "The method: '%s' does not exist for calculating the parabolic Fermi integrals. Please use 'numeric'. Exiting. ",
            method)
        sys.exit(1)
    sigma, seebeck, lorenz, hall, ccn, ccp = func[method](tr, eta, bs, tau0,
                                                          temperature)
    return sigma, seebeck, lorenz, hall, ccn, ccp


def parabolic_closed(tr, eta, bs, tau0_t, temperature):  # pylint: disable=too-many-locals
    r"""
    Calculates the parabolic Fermi integrals for the transport coefficients.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    eta : float
        The reduced chemical potential
    bs : object
        A `Bandstructure()` object containing the dispersion
        relations of the N included bands
    tau0_t : (N,M) ndarray
        The relaxation time approximation (RTA) prefactors (tau0) in
        units of fs for the N bands and M defined scattering
        mechanisms. This is converted to the r factor in order
        to use the closed Fermi integrals from tau0
    temperature : float
        The temperature in K.

    Returns
    -------
    tupple : ndarray, ndarray, ndarray, ndarray, ndarray, ndarray
        | Dimension: (3,3), (3,3), (3,3), (3,3), (3,3), (3,3)

        The electrical conductivity, Seebeck coefficient, Lorenz number,
        Hall coefficient (big R, where the small Hall factor is
        divided by the charge carrier concentration) and charge
        carrier concentration for n and ptype carriers in units
        of

        .. math:: \\mathrm{S}/\\mathrm{m}, \\mu \\mathrm{V}/\\mathrm{K},
                  10^{-8}\\mathrm{V}^2/\\mathrm{K}^2,
                  \\mathrm{cm}^3/\\mathrm{C},
                  10^{21} \\mathrm{cm}^{-3}, respectively.

    Notes
    -----
    The closed equations for the Fermi integrals can easily be
    developed from the following equaion from :func:`parabolic_numeric`

    .. math:: \\Sigma^i_{lm}=-\\frac{4e^2\\sqrt{m}}{3\\sqrt{2}\\pi^2
              \\hbar^3}\\int \\tau(E) E^{3/2}(E-\\mu)^i
              \\frac{\\partial f_0}{\\partial E} dE.

    We now chose the :math:`i` for the given tensor and enforce a
    scattering model of the form :math:`\\tau=\\tau_0E^{r-1/2}`,
    where :math:`\\tau_0` is constant in energy. Finally we use the
    product rule and expand the integral. This yields

    .. math:: \\Sigma^0=\\frac{4e^2\\sqrt{m}}{3\\sqrt{2}\\pi^2
              \\hbar^3}(r+1)F_r(E),

    .. math:: \\Sigma^1=\\frac{k_n}{e_n}\\left(\\frac{(r+2)
              F_{r+1}(E)}{(r+1)F_r(E)-\\mu} \\right),

    .. math:: \\Sigma^2=\\left(\\frac{k_n}{e_n}\\right)^2\\left(
              \\frac{(r+3)F_{r+2}(E)}{(r+1)F_r(E)} -
              \\frac{(r+2)F_{r+1}(E)}{(r+1)F_r(E)} \\right),

    where :math:`F_i(E)` is the famous Fermi integrals

    .. math:: F_i(E)=\\int f_0 E^i dE.

    It is customary to introduce the dimensionless energy
    :math:`\\epsilon=E/k_bT` and chemical potential
    :math:`\\eta=\\mu/k_bT`. We have the general relation

    .. math:: F_i(E)=\\beta^{-(i+1)}F_i(\\epsilon),

    and the :math:`\\Sigma` coefficients can be easily transformed.
    Please consult :func:`parabolic_numeric` for the specific
    transport tensors and units used.

    The execution of this method needs currently needs the effective
    mass tensor to be isotropic and only one scattering mechanism can
    be used per band. This method is selected by setting
    `transport_method` to "closed" in the general configuration file.

    See Also
    --------
    setup_scattering
    find_r_from_tau0
    parabolic_numeric

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.info("Running the closed analytic Fermi-Dirac integrals "
                "to calculate the transport coefficients.")

    raise NotImplementedError


def parabolic_numeric(tr, eta, bs, tau0_t, temperature):  # pylint: disable=too-many-locals
    r"""
    The solution of the energy integrals for the BTE RTA for a parabolic dispersion.

    Parameters
    ----------
    tr : object
        A `Transport()` object.
    eta : float
        The reduced chemical potential.
    bs : object
        A `Bandstructure()` object containing the energy dispersions
        for the N bands.
    tau0_t : (N,M) ndarray
        The relaxation time approximation (RTA) prefactors (tau0)
        in units of fs for the N bands and M defined scattering mechanisms.
    temperature : float
        The temperature in K.

    Returns
    -------
    tupple, ndarray, ndarray, ndarray, ndarray, ndarray
        | Dimension: (3,3), (3,3), (3,3), (3,3), (3,3), (3,3)

        The electrical conductivity, Seebeck coefficient, Lorenz number,
        Hall coefficient (big R, where the small Hall factor is
        divided by the charge carrier concentration) and charge carrier
        concentration for n and ptype carriers in units of

        .. math:: S/m, \\mu \mathrm{V}/\mathrm{K},
                  10^{-8}\mathrm{V}^2/\mathrm{K}^2,
                  \mathrm{cm}^3/\mathrm{C},
                  10^{21} \mathrm{cm}^{-3}`,

        respectively.

    Notes
    -----
    This routine is the same as :func:`parabolic_closed` except
    here we solve the Fermi intergrals numerically. This allows to
    use concatenated scattering mechanisms for each band. Otherwise
    the approximations and requirements are similar. The method
    is selected by setting `transport_method` to "numeric".

    The transport coefficients are defined as

    .. math:: \\Sigma^i_{lm}=-s\\frac{e^2}{8\\pi^3}\\int
              \\tau(E(\\vec{k}))v_l(E(\\vec{k}))
              v_m(E(\\vec{k}))(E(\\vec{k})-\\mu)^i
              \\frac{\\partial f_0}{\\partial
              E(\\vec{k})} d\\vec{k}.

    Using the fact that :math:`d\\vec{k}=k^2\\sin(\\theta)
    d\\theta d\\phi` we get

    .. math:: \\Sigma^i_{lm}=-s\\frac{e^2}{2\\pi^2}\\int
              k^2 \\tau(E(\\vec{k}))v_l(E(\\vec{k}))
              v_m(E(\\vec{k}))(E(\\vec{k})-\\mu)^i
              \\frac{\\partial f_0}{\\partial E(\\vec{k})} dk.

    For parabolic bands :math:`E=\\hbar^2 k^2/2m`. We also use
    :math:`E=mv^2/2`. Furthermore we assume that our crystal is
    isotropic and cubic, such that :math:`\\Sigma=\\Sigma_{00}`.
    Then :math:`v_i^2=(v_{xx}+v_{yy}+v_{zz})/3=v^2/3` and the
    expression above simplifies to

    .. math:: \\Sigma^i=-s\\frac{\\sqrt{2}e^2\\sqrt{m}}{3\\pi^2
              \\hbar^3}\\int \\tau(E) E^{3/2}(E-\\mu)^i
              \\frac{\\partial f_0}{\\partial E} dE

    And we have a very manageble integral over energy. Here, we
    have assumed that :math:`\\tau` can be expressed in terms of
    energy instead of the wave vector. A similar procedure can be
    used to obtain energy integrals for other dispersion relations
    than parabolic. As opposed to :func:`parabolic_closed` we here
    want to solve the integrals numerically (in order to be able
    to use composite :math:`\\tau`) and thus want to simplify
    the mathematical operations in the integrands as much
    as possible. Since

    .. math:: \\partial f_0/\\partial
             E=-\\beta/(2(1+\\cosh(\\beta(E-\\mu)))),

    where :math:`\\beta=(k_bT)^{-1}` we get

    .. math:: \\Sigma^i=s\\frac{\\sqrt{2m}e^2\\beta}{6\\pi^2
              \\hbar^3}\\int \\frac{\\tau(E) E^{3/2}(E-\\mu)^i}
              {1+\\cosh(\\beta(E-\\mu))} dE

    Notice now that there is a factor of 1/2 difference in front of
    these integrals compared to the ones in :func:`parabolic_closed`
    due to the expansion of the derivative of the Fermi function
    with respect to energy. The factor :math:`s=2` accounts for
    spin degeneracy, otherwise :math:`s=1` and is set to True or False
    for each band with the parameter `spin_deg` in the bandstructure
    configuration file.

    It is customary to introduce the dimensionless energy
    :math:`\\epsilon=E/k_bT` and chemical potential
    :math:`\\eta=\\mu/k_bT`. Doing this we obtain

    .. math:: \\Sigma^i=s\\frac{\\sqrt{2\\mathrm{m}}\\mathrm{e}^2}
              {6\\mathrm{\\hbar}^3 \\pi^2 \\beta^{i+3/2}}\\int
              \\frac{(\\epsilon-\\eta)^i\\tau(\\epsilon/\\beta)
              \\epsilon^{3/2}}{1+\\cosh(\\epsilon-\\eta)}d\\epsilon

    The carrier density :math:`n` can be calculated using

    .. math:: n=s\\frac{1}{8\\pi^3}\\int f_0 d\\vec{k}

    Depending on what kind of scattering models that is chosen,
    the rescaling of :math:`\\tau` needs to be performed for a
    particular scattering model in order to pull the correct
    :math:`\\beta` factor outside the integral.

    In order to obtain the correct units the equation above for
    each transport coefficient have been implemented as

    Electrical conductivity:

    .. math:: \\frac{s\\sqrt{20}}{6\\pi}G\\frac{\\sqrt{m_c}}{\hbar_c
              \\tilde{\hbar}}k^{3/2}a^{1/2}T^{3/2}
              \\int_0^{\\infty}\\frac{\\tau(\\epsilon/\\beta)
              \\epsilon^{3/2}}{1+\\cosh(\\epsilon-\\eta)}
              d\\epsilon \\left[ \\frac{\\mathrm{S}}{\\mathrm{m}}\\right].

    One usually set :math:`s=2`.
    Seebeck coefficient:

    .. math:: \\alpha=10^2\\frac{k_n}{e_n}\\frac{\\int_0^{\\infty}
              \\frac{\\tau(\\epsilon/\\beta)(\\epsilon-\\eta)
              \\epsilon^{3/2}}{1+\\cosh(\\epsilon-\\eta)}
              d\\epsilon}{\\int_0^{\\infty}
              \\frac{\\tau(\\epsilon/\\beta)\\epsilon^{3/2}}
              {1+\\cosh(\\epsilon-\\eta)}d\\epsilon}
              \\left[ \\frac{\\mu \mathrm{V}}{\mathrm{K}} \\right].

    Lorenz number:

    .. math:: L=\\left(\\frac{k_n}{e_n}\\right)^2\\left(
              \\frac{\\int_0^{\\infty}\\frac{
              \\tau(\\epsilon/\\beta)(\\epsilon-\\eta)^2
              \\epsilon^{3/2}}{1+\\cosh(\\epsilon-\\eta)}
              d\\epsilon}{\\int_0^{\\infty} \\frac{
              \\tau(\\epsilon/\\beta)\\epsilon^{3/2}}
              {1+\\cosh(\\epsilon-\\eta)}d\\epsilon}-
              \\left(\\frac{\\int_0^{\\infty}\\frac{
              \\tau(\\epsilon/\\beta)(\\epsilon-\\eta)
              \\epsilon^{3/2}}{1+\\cosh(\\epsilon-\\eta)}
              d\\epsilon}{\\int_0^{\\infty}\\frac{
              \\tau(\\epsilon/\\beta)\\epsilon^{3/2}}
              {1+\\cosh(\\epsilon-\\eta)}d\\epsilon}\\right)^2
              \\right) \\left[10^{-8}\\frac{\\mathrm{V^2}}
              {\\mathrm{K^2}}\\right].

    Hall coeffcient:

    .. math:: R_H=10^{-3}\\frac{3\\pi^2}{4\\sqrt{5}
              \\mathrm{e_n}}\\left(\\frac{\\mathrm{m_c k}}
              {\\mathrm{\hbar_c^2}}\\right)^{-3/2}
              \\left(aT\\right)^{-3/2}\\frac{\\int_0^{\\infty}
              \\frac{\\tau(\\epsilon/\\beta)^2\\epsilon^{3/2}}
              {1+\\cosh(\\epsilon-\\eta)}d\\epsilon}{
              \\left(\\int_0^{\\infty}\\frac{
              \\tau(\\epsilon/\\beta)\\epsilon^{3/2}}
              {1+\\cosh(\\epsilon-\\eta)}d\\epsilon\\right)^2}
              \\left[\\frac{\\mathrm{cm}^3}{\\mathrm{C}}\\right].

    Carrier concentration:

    .. math:: n=\\frac{20\\sqrt{5}}{\\pi^2}\\left(
              \\frac{\\mathrm{m_c k}}{\\mathrm{\hbar_c^2}}\\right)^{3/2}
              \\left(aT\\right)^{3/2}\\int_0^{\\infty}\\frac{
              \\sqrt{\\epsilon}}{1+\\cosh(\\epsilon-\\eta)}
              d\\epsilon \\left[10^{21} \\mathrm{cm}^{-3}\\right].

    From these, the Hall carrier concentration :math:`n_h` and
    Hall factor, :math:`r_h` can be calculated

    .. math:: n_h=\\frac{1}{eR_H},

    and

    .. math:: r_h=\\frac{n}{n_h}.

    The multiband expressions are accordingly (using :math:`i`
    as the band index)

    .. math:: n=\\sum_i n_i

    .. math:: \\sigma=\\sum_i \\sigma_i

    .. math:: \\alpha=\\frac{\\sum_i \\alpha_i \\sigma_i}
              {\\sum_i \\sigma_i}

    .. math:: L=\\frac{\\sum_i L_i \\sigma_i}{\\sum_i \\sigma_i}+
              \\frac{\\sum_i \\sigma_i \\alpha_i^2}
              {\\sum_i \\sigma_i}-\\alpha^2

    The units on the rest of the variables are

        | :math:`\\tau [\mathrm{fs}]`
        | The relaxation time given in fs.

        | :math:`T [\\mathrm{K}]`
        | The temperature given in K.

        | :math:`a`
        | The effective mass factor, e.g. :math:`m^*=am_e`.

    Furthermore, the coefficients that sets the correct
    scaling is defined as

        | :math:`k = 8.6173324`
        | The coefficient of the Boltzmann constant.

        | :math:`e_n = 2.417989348`
        | The coefficient for the normalized (e/h) electron charge.

        | :math:`m_c = 0.510998928`
        | The coefficient of :math:`mc^2`

        | :math:`\hbar_c = 197.3269718`
        | The coefficient of :math:`\\hbar c`

        | :math:`\\tilde{\hbar} = 6.582119514`
        | The coefficient of :math:`\\hbar`

        | :math::`G = 7.7480917346`
        | The coefficient of the conductance quantum

    .. warning:: ONLY VALID FOR PARABOLIC ENERGY DISPERSIONS
                 AND SCATTERING MODELS

    .. todo:: ADD POSIBILITY TO USE DIFFERENT EFFECTIVE MASSES ALONG
              DIFFERENT DIRECTIONS

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.info("Running numerical calculation of the energy intergrals "
                "in the parabolic approximation to calculate the transport "
                "coefficients.")
    numbands = tr.included_bands.size
    effmass = np.zeros(numbands)
    beta = 1e5 / (constants.kb * temperature)
    # it is easier to use 1/tau0 in the routine that
    w0 = np.nan_to_num(1.0 / tau0_t)
    # sums the combined scattering during integrations
    # (less operations), so invert
    # loop bands to check initials
    for band in range(numbands):
        band_actual = tr.included_bands[band]
        # check for parabolic effective mass
        effmass_vec = bs.effmass[band_actual]
        parabolic_effective_mass(effmass_vec)
        # and make sure it is positive
        effmass[band] = np.abs(effmass_vec[0])
        w0[band_actual] = w0[band_actual] * \
            np.power(beta, tr.scattering_r_factor - 0.5)
    sign = {"v": 1.0, "c": -1.0}
    # calculate contribution for each band
    sigma = np.zeros(numbands)
    seebeck = np.zeros(numbands)
    lorenz = np.zeros(numbands)
    hall = np.zeros(numbands)
    cc = np.zeros(numbands)
    # calculate single band coefficients, set integration limits
    e_min = 0.0
    e_max = tr.param.maxeint
    for band in range(numbands):
        band_actual = tr.included_bands[band]
        spin_deg = tr.bs.spin_degen[band_actual]
        tr.scattering_tau0_select = tr.bs.select_scattering[band_actual]
        sigmasigma = lbteint.scipy_e_integrals(
            tr, "normal", e_min, e_max, w0[band_actual], eta[band_actual],
            beta, tr.tau_energy_trans[band_actual], effmass[band], 0.0,
            spin_deg)
        sigmachi = lbteint.scipy_e_integrals(tr, "normal", e_min, e_max,
                                             w0[band_actual], eta[band_actual],
                                             beta,
                                             tr.tau_energy_trans[band_actual],
                                             effmass[band], 1.0, spin_deg)
        sigmakappa = lbteint.scipy_e_integrals(
            tr, "normal", e_min, e_max, w0[band_actual], eta[band_actual],
            beta, tr.tau_energy_trans[band_actual], effmass[band], 2.0,
            spin_deg)
        sigmahall = lbteint.scipy_e_integrals(tr, "hall", e_min, e_max,
                                              w0[band_actual],
                                              eta[band_actual], beta,
                                              tr.tau_energy_trans[band_actual],
                                              effmass[band], 0.0, spin_deg)
        sigmados = lbteint.scipy_e_integrals(tr, "dos", e_min, e_max,
                                             w0[band_actual], eta[band_actual],
                                             beta,
                                             tr.tau_energy_trans[band_actual],
                                             effmass[band], 0.0, spin_deg)
        # conductivity
        bandvaluesigma = math.sqrt(effmass[band]) * sigmasigma
        # for the seebeck, the units in front of the integral
        # in Sigma cancels, including the effective mass (thus
        # we only use sigmachi and sigmasigma for
        # the calculation of the Seebeck coefficient),
        # but we need a sign (n and p type)
        bandvalueseebeck = sign[tr.bs.status[band_actual]] * \
            (sigmachi / sigmasigma)
        # and now the lorenz
        bandvaluelorenz = sigmakappa / sigmasigma - \
            pow(sigmachi / sigmasigma, 2.0)
        # and hall
        bandvaluehall = sign[tr.bs.status[band_actual]] * \
            np.power(effmass[band], -1.5) * sigmahall / \
            np.power(sigmasigma, 2.0)
        # and cc
        bandvaluecc = np.power(effmass[band], 1.5) * sigmados
        # set sigma pr band (adjust units later)
        sigma[band] = bandvaluesigma
        # set seebeck pr band (adjust units later)
        seebeck[band] = bandvalueseebeck
        # set lorenz pr band (adjust units later)
        lorenz[band] = bandvaluelorenz
        # set hall pr band (adjust units later)
        hall[band] = bandvaluehall
        # set cc pr band (adjust units later)
        cc[band] = bandvaluecc

    # now we need to calculate the multiband totals, first conductivity
    sigma_total = np.sum(sigma)
    # some temporaries
    seebecktimessigma = np.multiply(seebeck, sigma)
    seebeck2timessigma = np.multiply(np.power(seebeck, 2.0), sigma)
    lorenztimessigma = np.multiply(lorenz, sigma)
    halltimessigma2 = np.multiply(hall, np.power(sigma, 2.0))
    # then seebeck
    seebeck_total = np.sum(seebecktimessigma) / sigma_total
    # then lorenz
    lorenz_total = (np.sum(lorenztimessigma) +
                    np.sum(seebeck2timessigma)) / \
        sigma_total - np.power(np.sum(seebecktimessigma) /
                               sigma_total, 2.0)
    # then hall
    hall_total = np.sum(halltimessigma2) / np.power(sigma_total, 2.0)

    # finally the cc, need to separate between n and p-type
    ntype_index = np.where(tr.bs.status == "c")
    ptype_index = np.where(tr.bs.status == "v")
    cc_total_n = np.sum(cc[ntype_index])
    cc_total_p = np.sum(cc[ptype_index])

    # now set units
    # units have been set from the parabolic_closed expressions
    # please notice that since we here use the cosh expression for the
    # derivative of df/dE the units in parabolic_numeric need to be
    # a factor of 0.5 lower than the units in front of the Fermi
    # integrals themself. This does not apply for the carrier
    # concentration, since there, only f enters.
    # Also notice that for the Seebeck and Lorenz this cancels, so we do
    # not need to add the 0.5 factor, while for the Hall coefficient,
    # the factor is 2.0, due to the fact that there is an integral
    # squared in the denominator while the numerator is not squared
    num_prefact_sigma = 0.5 * np.sqrt(20.0) / (3 * constants.pi)
    sigma_units = num_prefact_sigma * constants.g0 * constants.sqmcsq * \
        np.power(constants.kb, 1.5) * \
        np.power(temperature, 1.5) / constants.hbar
    num_prefact_seebeck = 100.0
    seebeck_units = (num_prefact_seebeck * constants.knorm / constants.enorm)
    lorenz_units = np.power(constants.knorm / constants.enorm, 2.0)  # pylint: disable=assignment-from-no-return
    # an additional factor of two is included below in order
    # to comply with e.g. the data given by May for the
    # Hall factor. Not confirmed, but suspect this is from
    # the spin degeneracy which is squared in the
    # denominator and thus needs to be multiplied away
    # TODO: Check this in more details pylint: disable=fixme
    # also add factor of two from above
    num_prefact_hall = 2.0 * 2.0 * 3.0 * 1e4 * \
        np.power(constants.pi, 2.0) / (2.0 * math.sqrt(2.0))
    hall_units = num_prefact_hall * np.power(constants.hbar, 3.0) / \
        (constants.elcharge * np.power(constants.elmass *
                                       constants.jtoev *
                                       constants.kb *
                                       temperature, 1.5))
    num_prefact_cc = 10.0 * math.sqrt(5.0) / np.power(math.pi, 2.0)
    cc_units = num_prefact_cc * \
        np.power(temperature, 1.5) * constants.scmcsqthird

    # and add units to the integrals
    sigma_total = sigma_units * sigma_total
    seebeck_total = seebeck_units * seebeck_total
    lorenz_total = lorenz_units * lorenz_total
    hall_total = hall_units * hall_total
    cc_total_n = cc_units * cc_total_n
    cc_total_p = cc_units * cc_total_p

    # now create the 3x3 tensors and fill the diagonal
    sigma_tensor = np.zeros((3, 3))
    seebeck_tensor = np.zeros((3, 3))
    lorenz_tensor = np.zeros((3, 3))
    hall_tensor = np.zeros((3, 3))
    cc_tensor_n = np.zeros((3, 3))
    cc_tensor_p = np.zeros((3, 3))
    np.fill_diagonal(sigma_tensor, sigma_total)
    np.fill_diagonal(seebeck_tensor, seebeck_total)
    np.fill_diagonal(lorenz_tensor, lorenz_total)
    np.fill_diagonal(hall_tensor, hall_total)
    np.fill_diagonal(cc_tensor_n, cc_total_n)
    np.fill_diagonal(cc_tensor_p, cc_total_p)

    return sigma_tensor, seebeck_tensor, lorenz_tensor, \
        hall_tensor, cc_tensor_n, cc_tensor_p


def numerick(tr, chempots, temperatures, bs=None):  # pylint: disable=too-many-locals, too-many-local-branches # noqa: MC0001
    r"""
    Calculates the transport coefficients according to the tensors defined in :func:`full_k_space_analytic`

    Parameters
    ----------
    tr : object
        A `Transport()` object.
    chempots : ndarray
        | Dimension: (N)

        The N chemical potentials in eV on which to calculate
        the transport coefficients.
    temperature : float
        | Dimension: (M)

        The M temperatutes in K.
    bs : ndarray, optional
        | Dimension: (I,J)

        The energy dispersion in eV for the charge carriers for I bands
        and J k-points. If not given, it defaults to the `Bandstructure()`
        object stored in `tr`.

    Returns
    -------
    sigma, seebeck, lorenz : ndarray, ndarray, ndarray
        | Dimension: (M,N,3,3)

        Returns the electrical condcutivity, Seebeck coefficient
        and Lorenz tensor for M temperature and N chemical potential
        steps in units of :math:`\\mathrm{S}/\\mathrm{m}`,
        :math:`\\mu \\mathrm{V}/\\mathrm{K}`,
        :math:`\\mathrm{V^2}/\\mathrm{K^2}`.

    Notes
    -----
    This routine accepts a predetermined array containg the
    relaxation time sampled at different carrier energy steps.
    When the evaluation of the intergrals are executed the carrier
    energy dispersion, velocity and scattering array are interpolated on
    all energies stored in `bs`
    and then the integrals are evaluated statically
    (if `transport_integration_method` is "trapz", "simps", "romb",
    "tetra" or "smeared"). The former method can also be used in the case
    where the scattering array is pretetermined by setting
    `transport_use_scattering_ontfly` in the general
    configuration file to True.

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running numerick.")
    if bs is None:
        bs = tr.bs
    energies = bs.energies
    velocities = bs.velocities
    # check if we need to generate velocities
    gen_velocities = int(tr.bs.gen_velocities)
    # get going on the transport tensor calculations
    sigma = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
    seebeck = np.zeros((sigma.shape))
    lorenz = np.zeros((sigma.shape))
    # check the velocities
    if not bs.check_velocities(constants.zerocut):
        if gen_velocities and not \
           bs.param.dispersion_velocities_numdiff:
            logger.info("No band velocities were supplied and we have no "
                        "way to calculate them.")
            sys.exit(1)
    # use simple static integration of the input grids
    # accuracy cannot be controlled, but it is rather fast and easy
    if (tr.param.transport_integration_method == "simps"  # pylint: disable=too-many-nested-blocks
            or tr.param.transport_integration_method == "trapz"
            or tr.param.transport_integration_method == "romb"):
        logger.info("Running simps, trapz or romb from Scipy to "
                    "calculate the transport coefficients.")
        # first check that the scattering array have been set up on
        # the energy grid
        scattering.check_scattering(tr)

        # currently method 2 is the fastest, but accuracy in
        # sigma flatlines at somewhere between 1e-4 and 1e-5
        # investigate further before switching entirely to method 2
        method = 1
        if method == 1:
            numbands = tr.included_bands.size
            sigma_band = np.zeros((numbands, 3, 3))
            seebeck_band = np.zeros((numbands, 3, 3))
            lorenz_band = np.zeros((numbands, 3, 3))
            # loop temperature and chemical potential manually
            # (broadcast later?)
            for tindex, temp in np.ndenumerate(temperatures):
                beta = 1e5 / temp / constants.kb
                # now set units and add temperature scaling
                sigma_units = 1e11 * constants.g0 / \
                    (16 * np.power(constants.pi, 2.0) *
                     constants.hbar * constants.kb) / temp
                seebeck_units = -1e6 / temp
                lorenz_units = 1e8 / np.power(temp, 2.0)
                for cindex, chempot in np.ndenumerate(chempots):
                    for band in range(numbands):
                        band_actual = tr.included_bands[band]
                        spin_deg = tr.bs.spin_degen[band_actual]
                        energies = tr.bs.energies[band_actual]
                        velocities = tr.bs.velocities[band_actual]
                        scatter = tr.scattering_total_inv[tindex[0],
                                                          band_actual]
                        sigmasigma = lbteint.scipy_k_integrals_discrete(
                            tr,
                            "normal",
                            energies,
                            velocities,
                            scatter,
                            chempot,
                            beta,
                            0.0,
                            spin_deg,
                            method=tr.param.transport_integration_method)
                        sigmachi = lbteint.scipy_k_integrals_discrete(
                            tr,
                            "normal",
                            energies,
                            velocities,
                            scatter,
                            chempot,
                            beta,
                            1.0,
                            spin_deg,
                            method=tr.param.transport_integration_method)
                        sigmakappa = lbteint.scipy_k_integrals_discrete(
                            tr,
                            "normal",
                            energies,
                            velocities,
                            scatter,
                            chempot,
                            beta,
                            2.0,
                            spin_deg,
                            method=tr.param.transport_integration_method)
                        # conductivity
                        bandvaluesigma = sigmasigma
                        # for the seebeck, the units in front of the integral
                        # in the Sigmas cancels
                        # remember to do matrix and not elementwise
                        # for all tensors
                        # check if sigma is singular
                        sigmainv = utils.invert_matrix(sigmasigma)
                        bandvalueseebeck = np.dot(sigmainv, sigmachi)
                        # and now the lorenz (same with the units)
                        bandvaluelorenz = (
                            np.dot(sigmakappa, sigmainv) - np.dot(
                                np.dot(sigmachi, bandvalueseebeck), sigmainv))
                        # set sigma pr band (adjust units later)
                        sigma_band[band] = bandvaluesigma
                        # set seebeck pr band (adjust units later)
                        seebeck_band[band] = bandvalueseebeck
                        # set lorenz pr band (adjust units later)
                        lorenz_band[band] = bandvaluelorenz

                    # now we need to calculate the multiband totals and
                    # add explicit temperature dependency, first conductivity
                    sigma_total = np.sum(sigma_band, axis=0)
                    # some temporaries
                    # do not print invalid multiplies
                    with np.errstate(invalid='ignore'):
                        seebecktimessigma = np.multiply(
                            seebeck_band, sigma_band)
                        seebeck2timessigma = np.multiply(
                            np.power(seebeck_band, 2.0), sigma_band)
                        lorenztimessigma = np.multiply(lorenz_band, sigma_band)
                    # then seebeck (remove nans etc.)
                    seebeck_total = np.nan_to_num(
                        np.sum(seebecktimessigma, axis=0) / sigma_total)
                    # then lorenz (remove nans etc.)
                    lorenz_total = np.nan_to_num(
                        (np.sum(lorenztimessigma, axis=0) +
                         np.sum(seebeck2timessigma, axis=0)) / sigma_total -
                        np.power(
                            np.sum(seebecktimessigma, axis=0) /
                            sigma_total, 2.0))

                    # and add units to the integrals
                    sigma[tindex, cindex] = sigma_units * sigma_total
                    seebeck[tindex, cindex] = seebeck_units * seebeck_total
                    lorenz[tindex, cindex] = lorenz_units * lorenz_total

        elif method == 2:
            numbands = tr.included_bands.size
            sigma_band = np.zeros((numbands, 3, 3))
            seebeck_band = np.zeros((numbands, 3, 3))
            lorenz_band = np.zeros((numbands, 3, 3))
            # loop temperature and chemical potential manually
            # (broadcast later?)
            energies = energies[tr.included_bands]
            velocities = velocities[tr.included_bands]
            spin_fact = tr.bs.spin_degen[tr.included_bands]
            for tindex, temp in np.ndenumerate(temperatures):
                beta = 1e5 / temp / constants.kb
                # now set units and add temperature scaling
                sigma_units = 1e11 * constants.g0 / \
                    (16 * np.power(constants.pi, 2.0) *
                     constants.hbar * constants.kb) / temp
                seebeck_units = -1e6 / temp
                lorenz_units = 1e8 / np.power(temp, 2.0)
                if tr.param.parallel:
                    # use mpi4py to distribute calculations over
                    # chemical potentials
                    # DOES NOT WORK AT THE MOMENT.
                    raise NotImplementedError(
                        'Parallelization is not yet implemented.')
                # serial version
                for cindex, chempot in np.ndenumerate(chempots):
                    scatter = tr.scattering_total_inv[tindex, tr.
                                                      included_bands]
                    sigmasigma = lbteint.scipy_k_integrals_discrete2(
                        tr,
                        energies,
                        velocities,
                        scatter,
                        chempot,
                        beta,
                        spin_fact,
                        0.0,
                        method=tr.param.transport_integration_method)
                    sigmachi = lbteint.scipy_k_integrals_discrete2(
                        tr,
                        energies,
                        velocities,
                        scatter,
                        chempot,
                        beta,
                        spin_fact,
                        1.0,
                        method=tr.param.transport_integration_method)
                    sigmakappa = lbteint.scipy_k_integrals_discrete2(
                        tr,
                        energies,
                        velocities,
                        scatter,
                        chempot,
                        beta,
                        spin_fact,
                        2.0,
                        method=tr.param.transport_integration_method)
                    # conductivity
                    sigma_nounit = sigmasigma
                    # for the seebeck, the units in front of the integral
                    # in the Sigmas cancels
                    # remember to do matrix and not elementwise
                    # for all tensors
                    # check if sigma is singular
                    sigmainv = utils.invert_matrix(sigmasigma)
                    seebeck_nounit = np.dot(sigmainv, sigmachi)
                    # and now the lorenz (same with the units)
                    lorenz_nounit = (
                        np.dot(sigmakappa, sigmainv) -
                        np.dot(np.dot(sigmachi, seebeck_nounit), sigmainv))
                    # and add units to the integrals
                    sigma[tindex, cindex] = sigma_units * sigma_nounit
                    seebeck[tindex, cindex] = seebeck_units * \
                        seebeck_nounit
                    lorenz[tindex, cindex] = lorenz_units * lorenz_nounit

    # here we use weighted integration, accuracy cannot be controlled,
    # but it is rather fast and easy, currently linear tetrahedron and
    # smeared weights are implemented.
    elif ((tr.param.transport_integration_method == "tetra")
          or (tr.param.transport_integration_method == "smeared")):

        # lazy import
        import t4me.spglib_interface as spglib_interface  # pylint: disable=import-error, no-name-in-module
        # this method needs the IBZ. This is not available if we work
        # on the full grid.
        if tr.param.work_on_full_grid:
            logger.error("The tetrahedron and smeared integration method "
                         "can not work if the IBZ grid is not present. Please "
                         "choose one of the other integration techniques. "
                         "Exiting.")
            sys.exit(1)
        # first check that the scattering array have been set up on
        # the energy grid
        scattering.check_scattering(tr)
        # this parameter controls the energy integral, but
        # this is redudant as only trapezoidal is currently implemented
        int_method = 0
        energy_cutoff = tr.param.transport_integration_spectral_energy_cutoff
        energy_samples = tr.param.transport_integration_spectral_density
        volume = np.dot(
            tr.lattice.unitcell[0],
            np.cross(tr.lattice.unitcell[1], tr.lattice.unitcell[2]))
        scatter = tr.scattering_total_inv

        # now check if we actually want tetrahedron weight, histogram or
        # smeared
        if tr.param.transport_integration_method == "tetra":
            logger.info("Running linear tetrahedron method to calculate "
                        "the transport coefficients.")
            weight_type = 0
            smearing = 0.0
        # THIS SHOULD BE REMOVED AS THIS IN THE END BASICALLY
        # BASICALLY YIELDS THE TRAPEZOIDAL RULE IN THE K SPACE
        if tr.param.transport_integration_method == "histogram":
            weight_type = 1
            smearing = 0.0
        if tr.param.transport_integration_method == "smeared":
            logger.info("Running smeared integration method to calculate "
                        "the transport coefficients.")
            weight_type = 2
            smearing = tr.param.transport_integration_spectral_smearing
        # due to the fact that the velocities are different for supposedly
        # similar ibz-bz points we set the ibz_weights table to 1.0 and pass
        # that
        # ibz_weights_dummy = np.ones(
        #    tr.lattice.mapping_ibz_to_bz.size, dtype='intc')
        energies_ibz = np.take(tr.bs.energies,
                               tr.lattice.mapping_ibz_to_bz,
                               axis=1)
        spglib_interface.calc_transport_tensors_weights_interface(  # pylint: disable=c-extension-no-member
            energies_ibz, tr.bs.velocities, scatter, temperatures, chempots,
            tr.lattice.spg_kmesh,
            np.ascontiguousarray(tr.lattice.mapping_bz_to_ibz, dtype="intc"),
            np.ascontiguousarray(tr.lattice.mapping_ibz_to_bz, dtype="intc"),
            np.ascontiguousarray(tr.lattice.ibz_weights, dtype="intc"),
            tr.lattice.ksampling, tr.lattice.runitcell,
            tr.bs.energies.shape[0], tr.lattice.mapping_ibz_to_bz.size,
            tr.temperatures.shape[0], tr.chempots.shape[0], int_method,
            energy_samples, weight_type, smearing, energy_cutoff, volume,
            tr.bs.spin_degen, sigma, seebeck, lorenz)
    else:
        logger.error("The supplied parameter for transport_integration_method "
                     "is not supported. Consult the documentation. Exiting.")
        sys.exit(1)

    return sigma, seebeck, lorenz


def calculate_hall_carrier_concentration(hall):
    r"""
    Calculates the Hall carrier concentration :math:`n_h=(\\mathrm{e}R_h)^{-1}`.

    Parameters
    ----------
    hall : ndarray
        | Dimension: (N,M,3,3)

        The Hall tensor :math:`R_h` in units of
        :math:`\\mathrm{cm^3}/\\mathrm{C}` for
        N temperature and M chemical potential samplings.

    Returns
    -------
    ndarray
        | Dimension: (N,M,3,3)

        The Hall carrier concentration in units of
        :math:`10^{-21} \mathrm{cm}^{-3}`

    """

    return 0.01 * np.linalg.inv(hall[:, :]) / constants.elcharge


def calculate_hall_factor(n, nh):
    r"""
    Calculates the Hall factor :math:`r_h=n/n_h`

    Parameters
    ----------
    n : ndarray
        | Dimension: (N,M,3,3)

        The calcalated carrier concentration in for N
        temperature and M chemical potential
        samplings in units of :math:`10^{21} \\mathrm{cm}^{-3}`.

    n : ndarray
        | Dimension: (N,M,3,3)

        The Hall carrier concentration in for N temperature
        and M chemical potential
        samplings in units of :math:`10^{21} \\mathrm{cm}^{-3}`.

    Returns
    -------
    ndarray
        | Dimension: (N,M,3,3)

        The Hall factor.
    """

    inv_nh = np.linalg.inv(nh[:, :])
    return inv_nh * n
