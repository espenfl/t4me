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
"""Contains routines to set up the scattering of the charge carriers."""

# pylint: disable=useless-import-alias, too-many-arguments, invalid-name, too-many-statements, too-many-lines, global-statement, too-many-nested-blocks

import sys
import logging
import numpy as np
import scipy

import t4me.constants as constants
from t4me.bandstructure import parabolic_effective_mass


def scattering_dos(tr, dos, energies, select_scattering):  # pylint: disable=too-many-locals, too-many-branches # noqa: MC0001
    """
    Setup scattering mechnisms.

    Store values in the scattering arrays using the density of states data
    as the energy dependency.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    dos : ndarray
        | Dimension: (N,M)

        Array containing the partial density of states
        (1/eV/AA^3), where N is the band index
        and M is the energy index.
    energies : ndarray
        | Dimension: (M)

        Array containing the energy in eV at M samplings
        where the density of states is calculated.
    select_scattering : ndarray
        | Dimension: (12)

        Array containing integers. Set to 1 to select
        the scattering, 0 to exclude.
        The variables in `select_scattering` are set in
        the bandstructure configuration file, one value
        for each scattering and band. See notes below
        for the currrently available scattering mechnisms.

    Returns
    -------
    scattering_inv : ndarray
        | Dimension: (T,N,M,12)

        The scattering array in fs units
        in the current `Transport()` object for T temperature steps,
        N number of bands, M number of energy steps and 12 number of
        scattering mechanisms
    scattering_total_inv : ndarray
        | Dimension: (T, N, M)

        The total (all mechanisms summed) scattering array
        in fs units) in the current `Transport()` object for T
        temperature steps, N number of bands and M number of
        energy steps
    scattering_tau0 : ndarray
        | Dimension: (T, N, 12)

        The scattering prefactor array, tau0 in units of fs,
        in the current `Transport()` object for T temperature
        steps, N number of bands and 12 number of
        scattering mechanisms.

    Notes
    -----
    Currently only the following scattering mechanisms are supported:

    ========================= ====================
    `select_scattering` index scattering mechanism
    ========================= ====================
    1                         Acoustic phonon scattering from def. pot.
    2                         Non-polar optical phonon scattering from def. pot. (alpha stage)
    3                         Intervalley phonon scattering (alpha stage)
    4                         Polar optical phonon scattering (alpha stage)
    5                         None
    6                         None
    7                         None
    8                         None
    9                         None
    10                        None
    11                        None
    12                        Constant (energy and k-point independent)
    ========================= ====================

    Only the acoustic phonon scattering has been tested.

    Please consult the book of C. Jacoboni :cite:`jacoboni_2010_toetis` and the
    publication concerning this software by
    E. Flage-Larsen :cite:`t4me` for additional details
    on the involved parameters. Also consult the bandstructure
    configuration file for the respective constants that have to
    be set besides `select_scattering` and their units.

    .. todo:: Add more extensive documentation for the different scattering
              mechanisms.

    .. warning:: The scattering models based on density of states
                 does currently not properly involve the energy
                 shift required for the transfer energies.
                 This is quite serious, but does not influene the acoustic
                 phonon scattering (taken to be zero in the model implemented
                 here). The current approach if not using the
                 analytic parabolic models is that the transfer energies
                 for all scattering mechanisms are summed and the energy
                 of where the relaxation time is evaluated is shifted by
                 this amount during interpolation. The sum approximation
                 is not physically justified and needs additional
                 investigation, also the interpolation are sensitive and
                 can fail close to the van Hove singularities since the
                 relaxation time is propotional to the inverse of
                 the density of states. All known problems pertaining to
                 density of states are also manifested here for the
                 scattering.

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running scattering_dos.")
    logger.debug("Calculating the scattering properties based on the "
                 "models of the density of states.")

    # set temperature
    temperatures = tr.temperatures

    num_scatterings = select_scattering[0].shape[0]
    num_bands = tr.bs.energies.shape[0]
    num_energy_steps = energies.shape[0]
    temperature_steps = temperatures.shape[0]
    prefactor_scattering = np.zeros((temperature_steps, num_bands,
                                     num_scatterings))
    scattering = np.zeros((temperature_steps, num_bands, num_energy_steps,
                           num_scatterings))
    inc_in_total = np.zeros((num_bands, num_scatterings), dtype=bool)
    tr.tau_energy_trans = np.zeros((num_bands, num_scatterings))
    # prepare stuff that does not depend on temperature
    # calculate q for intervalley phonon scattering
    q_diff = tr.bs.q_energy_trans[:, 1] - tr.bs.q_energy_trans[:, 0]
    q_length = np.linalg.norm(q_diff)

    # build scaling (unit array)
    scaling = np.full((num_scatterings, 2), constants.zeroshift)
    scaling[0, 0] = 1e5 * constants.pi * \
        constants.kb / (constants.hbar * constants.jtoev)
    scaling[1, 0] = 1e2 * constants.pi / (2.0 * constants.jtoev)
    scaling[1, 1] = 1e-4 * constants.hbar
    scaling[2, 0] = 1e5 * constants.pi / (constants.hbar * constants.jtoev)
    scaling[2, 1] = 1e2 * constants.pi / (2.0 * constants.jtoev)
    scaling[3, 0] = 2 * constants.pi
    # check that we do not loop over several temperatures
    # if one or more explicit tau0 has been given, since tau0
    # often include a temperature dependent factor
    if np.any(tr.bs.explicit_prefact):
        if temperatures.shape[0] > 1:
            logger.error("Explicit tau0 have been set, but the user "
                         "also wants to use these at different "
                         "temperatures. This is often simply wrong. "
                         "Exiting.")
            sys.exit(1)

    # build prefix array
    # loop temperature
    for tempi, tempv in np.ndenumerate(temperatures):
        # loop scattering mechanisms
        for band in range(num_bands):
            emi = 0.0
            sign = 1.0
            # fetch which scattering processes set for total sum
            inc_in_total[band] = tr.bs.select_scattering[band]
            for scattering_index, _ in \
                    np.ndenumerate(select_scattering[band]):

                # This portion sets the scattering array based
                # on the more general density of states models.

                # tau=tau_0/DOS(E+E_trans),

                # where E_trans is some transfer energy.

                # However, in this routine we use the inverse values
                # (scattering rates), w (propto DOS(E))
                if scattering_index[0] == 0:
                    # Elastic acoustic phonon scattering

                    # w0[0] = pi*k*T*D_a^2/(hbar*v^2*rho)

                    # D_a = acoustic deformation potential [eV]
                    # v = speed of sound [m/s]
                    # rho = mass density [g/cm^3]

                    # and DOS in units of 1/(eV AA^3) yields w0
                    # in units of fs as for the parabolic case

                    # tau0[0] = pi*D^2/(hbar*v^2*rho)

                    # units are the same, but the kT factor is added
                    # on the fly to account for temperature
                    # dependent data arrays
                    if tr.bs.explicit_prefact[band][0]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][0]
                    else:
                        da = tr.bs.da[band]
                        rho = tr.bs.rho[band]
                        speed_sound = tr.bs.speed_sound[band]
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            scaling[0, 0] * tempv * np.power(da /
                                                             speed_sound,
                                                             2.0) / rho

                elif scattering_index[0] == 1:
                    # Nonpolar optical phonon scattering

                    # w0[1]=pi*D_o^2*(n_op+1/2-/+1/2)/(2*rho*omega_op)

                    # D_o = non-polar optical deformation potential [eV/AA]
                    # n_op = optical phonon occupation number
                    # rho = mass density [g/cm^3]
                    # omega_op = optical phonon angular frequency
                    # (often Einstein) [THz]
                    if tr.bs.explicit_prefact[band][1]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][1]
                    else:
                        d = tr.bs.do[band]
                        n = tr.bs.no[band]
                        rho = tr.bs.rho[band]
                        omega = tr.bs.omegao[band]
                        tr.tau_energy_trans[
                            band, scattering_index] = \
                            sign * scaling[1, 1] * omega

                        if tr.bs.emi[band]:
                            emi = 1.0
                        prefactor_scattering[
                            tempi, band, scattering_index] = \
                            scaling[1, 0] * np.power(d, 2.0) * \
                            (n + 1.0 * emi) / (rho * omega)

                elif scattering_index[0] == 2:
                    # Intervalley phonon scattering

                    # w0[2]=pi*D_vv'^2*(n_vv'+1/2-/+1/2)/(2*rho*omega_vv')

                    # D_vv'=sqrt(|D_a*q_vv'|^2+D_o^2)

                    # D_a = acoustic deformation potential [eV]
                    # D_o = non-polar optical deformation potential [eV/AA]
                    # n_vv' = intervalley phonon occupation number
                    # rho = mass density [g/cm^3]
                    # omega_vv' = intervalley transition phonon angular
                    #             frequency [THz]
                    if tr.bs.explicit_prefact[band][2]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][2]
                    else:
                        # calculate D_vv'
                        da = tr.bs.da[band]
                        do = tr.bs.do[band]
                        n = tr.bs.nvv[band]
                        rho = tr.bs.rho[band]
                        omega = tr.bs.omegavv[band]
                        if tr.bs.emi[band]:
                            emi = 1.0
                        dvv = np.sqrt(
                            np.power(scaling[2, 0] * da * q_length, 2.0) +
                            np.power(scaling[2, 1] * do, 2.0))
                        prefactor_scattering[
                            tempi, band, scattering_index] = \
                            scaling[2, 1] * np.power(dvv, 2.0) * \
                            (n + 1.0 * emi) / (rho * omega)

                elif scattering_index[0] == 3:
                    # Polar optical phonon scattering

                    # w0[3]=2*pi*e^2*F^2*(n_op+1/2-/+1/2/hbar

                    # F = Frohlich expression [?]
                    # n_op = optical phonon occupation number
                    if tr.bs.explicit_prefact[band][3]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][3]
                    else:
                        # GET RID OF FROHLICH AND DO THE WHOLE THING
                        f = tr.bs.f[band]
                        n = tr.bs.no[band]
                        # NO OMEGA YET, COMES FROM THE FROHLICH ++
                        tr.tau_energy_trans[
                            band,
                            scattering_index] = sign * scaling[1, 1] * omega
                        if tr.bs.emi[band]:
                            emi = 1.0
                        prefactor_scattering[
                            tempi, band, scattering_index] = \
                            scaling[3, 0] * np.power(f, 2.0) * \
                            (n + 1.0 * emi)

    # BEWARE: HERE WE DO NOT SHIFT PROPERLY FOR THE
    # TRANSFER ENERGIES, PRINT WARNING
    # UNTIL WE HAVE A SOLUTION TO THIS PROBLEM
    logger.warning("BEWARE: The scattering models based on "
                   "density of states does currently not "
                   "properly involve the energy shift required "
                   "for the transfer energies. This is quite "
                   "serious, but does not influene the acoustic "
                   "phonon scattering. The current approach if "
                   "not using the analytic parabolic models is "
                   "that the transfer energies for all scattering "
                   "mechanisms are summed and the energy of "
                   "where the relaxation time is evaluated "
                   "is shifted by this amount during "
                   "interpolation. The sum approximation is not "
                   "physically justified and needs additional "
                   "investigation, also the interpolation "
                   "are sensitive and can fail close to the van "
                   "Hove singularities since the relaxation "
                   "time is propotional to the inverse of the "
                   "density of states. Continuing.")
    # now, make sure tau_energy_trans is zero for the
    # scattering mechnisms we do not
    # want (because of the sum later runs over the whole array)
    # and values are set
    # even though select_scattering is False for these entries
    # ~ flips the bool values in the array
    tr.tau_energy_trans[~select_scattering] = 0.0

    # now squeeze in the energy dependence
    # (either directly or indirectly)
    # multiply dos_times_scattering with the prefactor
    # to obtain the scattering array
    scattering[:, :, :, 0:num_scatterings - 2] = \
        dos[np.newaxis, :, :, np.newaxis] * \
        prefactor_scattering[:, :, np.newaxis,
                             0:num_scatterings - 2]
    # now add the constant scattering part given in fs units
    # (so invert)
    prefactor_scattering[:, :, num_scatterings -
                         1] = 1.0 / tr.bs.tau0c[np.newaxis, :]
    scattering[:, :, :, num_scatterings - 1] = 1.0 / \
        tr.bs.tau0c[np.newaxis, :, np.newaxis]
    # set up array to force zeros into the sum array if one
    # only wants certain scattering mechanisms in the total sum
    iit = inc_in_total[:, np.newaxis, :] * \
        np.ones(energies.shape[0], dtype=int)[np.newaxis, :, np.newaxis]
    # now calculate the total scattering rate
    scattering_total = (np.nan_to_num(scattering) * iit).sum(-1)
    # and then, since up til now we have calculated the scattering rate
    # we invert all values in order to get the "tau" in fs units
    # also set true zero to a small value before inverting
    scattering = np.nan_to_num(scattering)
    scattering[scattering < constants.zero] = constants.zero
    scattering_total[scattering_total < constants.zero] = constants.zero
    scattering_inv = 1.0 / scattering
    scattering_inv = np.nan_to_num(scattering_inv)
    with np.errstate(over="ignore"):
        scattering_total_inv = 1.0 / scattering_total
    # remove nan, inf etc. and return
    scattering_inv = np.nan_to_num(scattering_inv)
    scattering_total_inv = np.nan_to_num(scattering_total_inv)
    # tau0 for use in the closed Fermi integrals
    prefactor_scattering[
        prefactor_scattering < constants.zero] = constants.zero
    scattering_tau0 = np.nan_to_num(1.0 / prefactor_scattering)
    # now make sure the non selected scattering mechanisms
    # contain very large value (so
    # that the scattering rate W=0)
    non_selected = np.array(1 - iit[:, 0, :], dtype=bool)
    scattering_tau0[:, non_selected] = constants.large
    return scattering_inv, scattering_total_inv, scattering_tau0


def scattering_parabolic(tr, energies, select_scattering, use_eonk=False):  # pylint: disable=too-many-locals # noqa: MC0001
    """
    Setup scattering mechnisms.

    Store values in the scattering arrays using parabolic band dispersions
    as an approximation.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    energies : ndarray
        | Dimension: (N)

        Array containing the energy in eV at N samplings
        where the scattering values are o be calculated.
    select_scattering : ndarray
        | Dimension: (12)

        Array containing integers. Set to 1 to select
        the scattering, 0 to exclude.
        The variables in `select_scattering` are set in
        the bandstructure configuration file, one value
        for each scattering and band. See notes below
        for the currrently available scattering mechnisms.
    use_eonk : boolean
        If set to True, generate the scattering values on the
        supplied energy for each band and on its k-points

    Returns
    -------
    scattering_inv : ndarray
        | Dimension: (T,N,M,12)

        The scattering array in fs units
        in the current `Transport()` object for T temperature steps,
        N number of bands, M number of energy steps and 12 number of
        scattering mechanisms
    scattering_total_inv : ndarray
        | Dimension: (T, N, M)

        The total (all mechanisms summed) scattering array
        in fs units) in the current `Transport()` object for T
        temperature steps, N number of bands and M number of
        energy steps
    scattering_tau0 : ndarray
        | Dimension: (T, N, 12)

        The scattering prefactor array, tau0 in units of fs,
        in the current `Transport()` object for T temperature
        steps, N number of bands and 12 number of
        scattering mechanisms.

    Notes
    -----
    Currently only the following scattering mechanisms are supported:

    ========================= ====================
    `select_scattering` index scattering mechanism
    ========================= ====================
    1                         Acoustic phonon scattering from def. pot.
    2                         Non-polar optical phonon scattering from def. pot.
    3                         Intervalley phonon scattering
    4                         Polar optical phonon scattering
    5                         Piezoelectric acoustic phonon scattering
    6                         Ionized impurity scattering, Brooks-Herring
    7                         Ionized impurity scattering, Conwell-Weisskopf
    8                         Alloy scattering
    9                         None
    10                        None
    11                        None
    12                        Constant (energy and k-point independent)
    ========================= ====================

    Please consult the book of C. Jacoboni :cite:`jacoboni_2010_toetis` and the
    publication concerning this software by
    E. Flage-Larsen :cite:`t4me` for additional details
    on the involved parameters. Also consult the bandstructure
    configuration file for the respective constants that have to
    be set besides `select_scattering` and their units.

    .. todo:: Add more extensive documentation for the different scattering
              mechanisms.

    .. rubric:: References

    .. bibliography:: references.bib
        :style: unsrt
        :filter: docname in docnames

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running scattering_parabolic.")
    logger.debug("Calculating the scattering properties based on the "
                 "analytic parabolic scattering models.")

    # set temperatures
    temperatures = tr.temperatures

    num_scatterings = select_scattering[0].shape[0]
    num_bands = tr.bs.energies.shape[0]
    if not use_eonk:
        # use supplied energy grid
        num_energy_steps = energies.shape[0]
    else:
        # use the k-point grid of the energies
        # assume same set of k-point grid for all bands
        num_energy_steps = tr.bs.energies.shape[1]
    factor = np.zeros((num_bands, num_scatterings))
    temperature_steps = temperatures.shape[0]
    prefactor_scattering = np.zeros((temperature_steps, num_bands,
                                     num_scatterings))
    energy_correction_prefactor = np.zeros((num_bands, num_scatterings))
    energy_r_correction = np.zeros((num_bands, num_scatterings))
    scattering = np.zeros((temperature_steps, num_bands, num_energy_steps,
                           num_scatterings))
    inc_in_total = np.zeros((num_bands, num_scatterings), dtype=bool)
    tr.tau_energy_trans = np.zeros((num_bands, num_scatterings))
    # prepare stuff that does not depend on temperature
    # calculate q for intervalley phonon scattering
    q_diff = tr.bs.q_energy_trans[:, 1] - tr.bs.q_energy_trans[:, 0]
    # remember to bring q to cartesian
    q_length = np.linalg.norm(tr.bs.lattice.dir_to_cart(q_diff))

    # build scaling (unit array)
    scaling = np.full((num_scatterings, 2), constants.zeroshift)
    r_factor = np.zeros(num_scatterings)
    scaling[0, 0] = 1e4 * constants.kb * \
        np.power(constants.elmass, 1.5) * np.sqrt(constants.jtoev) / \
        (np.sqrt(5) * constants.pi * np.power(constants.hbar, 4.0))
    scaling[1, 0] = np.sqrt(5) * np.power(constants.elmass, 1.5) * \
        np.sqrt(constants.jtoev) / \
        (constants.pi * np.power(constants.hbar, 3.0))
    scaling[1, 1] = 1e-4 * constants.hbar
    scaling[2, 0] = scaling[1, 0]
    scaling[2, 1] = scaling[1, 1]
    scaling[3, 0] = 0.1 * np.sqrt(constants.elmass) * \
        np.power(constants.elcharge, 2.0) * \
        constants.vacperm / \
        (4 * np.sqrt(20) * constants.pi * constants.hbar)
    scaling[3, 1] = scaling[1, 1]
    scaling[4, 0] = 1e7 * np.sqrt(constants.elmass) * \
        np.power(constants.elcharge, 2.0) * \
        constants.kb / \
        (np.sqrt(80) * constants.pi *
         np.power(constants.vacperm *
                  constants.hbar, 2.0))
    scaling[4, 1] = constants.bandunit
    scaling[5, 0] = 1e2 * np.sqrt(10) * \
        np.power(constants.elcharge, 4.0) / \
        (32 * np.sqrt(2) * constants.pi *
         np.power(constants.vacperm, 2.0) *
         np.sqrt(constants.elmass) *
         np.power(constants.jtoev, -1.5))
    scaling[5, 1] = constants.bandunit
    scaling[6, 0] = 1e2 * np.sqrt(10) * \
        np.power(constants.elcharge, 4.0) / \
        (32 * np.sqrt(2) * constants.pi *
         np.power(constants.vacperm, 2.0) *
         np.sqrt(constants.elmass) * np.power(constants.jtoev, -1.5))
    scaling[6, 1] = 1e-2 * 64 * np.power(constants.pi, 2.0) * \
        np.power(constants.vacperm, 2.0) / \
        np.power(constants.elmass * constants.jtoev, 2.0)
    scaling[7, 0] = np.power(constants.elmass, 1.5) / \
        (np.power(constants.hbar, 4.0) * np.sqrt(5))
    r_factor[0] = 0.0
    r_factor[1] = 0.0
    r_factor[2] = 0.0
    r_factor[3] = 1.0
    r_factor[4] = 1.0
    r_factor[5] = 2.0
    r_factor[6] = 2.0
    r_factor[7] = 0.0
    r_factor[num_scatterings - 1] = 0.5
    tr.scattering_r_factor = r_factor
    r_factor_includinghalf = 0.5 - r_factor
    # check that we do not loop over several temperatures
    # if one or more explicit tau0 has been given, since tau0
    # often include a temperature dependent factor
    if np.any(tr.bs.explicit_prefact):
        if temperatures.shape[0] > 1:
            logger.error("Explicit tau0 have been set, but the user "
                         "also wants to use these at different "
                         "temperatures. This is often simply wrong. "
                         "Exiting.")
            sys.exit(1)

    # build prefix array
    # loop temperature
    for tempi, tempv in np.ndenumerate(temperatures):
        # loop scattering mechanisms
        for band in range(num_bands):
            emi = 0.0
            sign = 1.0
            # check parabolic effmass for all bands
            effmass_vec = tr.bs.effmass[band]
            if not parabolic_effective_mass(effmass_vec):
                logger.error("The setup of scattering mechanisms using "
                             "parabolic models requires a parabolic "
                             "effective mass. Exiting.")
                sys.exit(1)
            effmass = effmass_vec[0]
            # make sure effective mass is positive
            effmass = abs(effmass)
            # fetch which scattering processes set for total sum
            inc_in_total[band] = tr.bs.select_scattering[band]
            for scattering_index, _ in \
                    np.ndenumerate(select_scattering[band]):
                # This portion sets the scattering array based
                # on the well established parabolic scattering models

                # tau=tau_0*E^{r-1/2}

                # However, in this routine we use the inverse values
                # (scattering rates), w and invert at the end
                if scattering_index[0] == 0:
                    # Elastic acoustic phonon scattering

                    # r = 0

                    # w0=(sqrt(2)*m^3/2*k*T*D^2)/(pi*hbar^4*rho*v^2)

                    # (with spin degen included)

                    # energy dep = E^1/2

                    # m = effective mass [kg]
                    # D = acoustic deformation potential [eV]
                    # v = speed of sound [m/s]
                    # rho = mass density [g/cm^3]

                    # BOTH ABSORPTION AND EMMISION

                    # ONLY ELASTIC, e.g. hbar q v << kT
                    if tr.bs.explicit_prefact[band][0]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][0]
                    else:
                        if tempi[0] == 0:
                            da = tr.bs.da[band]
                            speed_sound = tr.bs.speed_sound[band]
                            rho = tr.bs.rho[band]
                            factor[band, scattering_index] = \
                                scaling[0, 0] * np.power(effmass, 1.5) * \
                                np.power(da, 2.0) / \
                                (np.power(speed_sound, 2.0) * rho)
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index] * tempv

                if scattering_index[0] == 1:
                    # Nonpolar optical scattering

                    # r = 0

                    # w0=(m^3/2*D_o^2)/(sqrt(2)*pi*hbar^3*rho*omega_op)
                    #      *(n_op+1/2-/+1/2)

                    # energy dep = (E +/- hbar*omega_op)^1/2

                    # m = effective mass, units of m_e
                    # D = non-polar optical deformation potential [eV/AA]
                    # n_op = optical phonon occupation number
                    # rho = mass density [g/cm^3]
                    # omega_op = optical phonon angular frequency
                    #            (often Einstein) [THz]

                    # NO EMMISION IF E<hbar omega
                    if tr.bs.explicit_prefact[band][1]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][1]
                    else:
                        if tempi[0] == 0:
                            do = tr.bs.do[band]
                            rho = tr.bs.rho[band]
                            omega = tr.bs.omegao[band]
                            n = tr.bs.no[band]
                            if tr.bs.emi[band]:
                                emi = 1.0
                                sign = -1.0
                            temp = sign * scaling[1, 1] * omega
                            energy_r_correction[band, scattering_index] = temp
                            # set the transfer energy (used sometimes when
                            # integrating the explicit tau)
                            tr.tau_energy_trans[
                                band, scattering_index] = \
                                sign * scaling[1, 1] * omega
                            factor[band, scattering_index] = \
                                scaling[1, 0] * np.power(effmass, 1.5) * \
                                np.power(do, 2.0) * (n + 1.0 * emi) / \
                                (omega * rho)
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index]

                if scattering_index[0] == 2:
                    # Intervalley phonon scattering

                    # r = 0

                    # w0=m^3/2*D_vv'^2*Z_f)/(sqrt(2)*pi*hbar^3*rho*omega_vv')
                    #      *(n_vv'+1/2-/+1/2)

                    # energy dep = (E +/- hbar*omega_vv' - dE_vv')^1/2

                    # m = effective mass, units of m_e
                    # Z_f = numer of possible final states (final degeneracy)
                    # D_vv'=sqrt(|D_a*q_vv'|^2+D_o^2)
                    # D_a = acoustic deformation potential [eV]
                    # D_o = non-polar optical deformation potential [eV/AA]
                    # n_vv' = intervalley phonon occupation number
                    # rho = mass density [g/cm^3]
                    # omega_vv' = intervalley transition phonon angular
                    #             frequency [THz]

                    # NO ABSORPTION OR EMMISION IF E < hbar omega - etrans
                    # (i=initial, f=final)

                    # where

                    # etrans = energy difference between the bottoms of the
                    #          final and initial valley [eV]
                    if tr.bs.explicit_prefact[band][2]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][2]
                    else:
                        # calculate D_vv'
                        if tempi[0] == 0:
                            da = tr.bs.da[band]
                            do = tr.bs.do[band]
                            n = tr.bs.nvv[band]
                            rho = tr.bs.rho[band]
                            omega = tr.bs.omegavv[band]
                            zf = tr.bs.zf[band]
                            etrans = tr.bs.etrans[band]
                            if tr.bs.emi[band]:
                                emi = 1.0
                                sign = -1.0
                            # TODO: generalize D_a*q to vector form pylint: disable=fixme
                            dvv = np.sqrt(
                                np.power(da * q_length, 2.0) +
                                np.power(do, 2.0))
                            temp = sign * scaling[2, 1] * omega - etrans
                            energy_r_correction[band, scattering_index] = temp
                            # set the transfer energy (used sometimes when
                            # integrating the explicit tau)
                            tr.tau_energy_trans[band, scattering_index] = temp
                            factor[band, scattering_index] = \
                                scaling[2, 0] * np.power(effmass, 1.5) * \
                                np.power(dvv, 2.0) * zf * (n + 1.0 * emi) / \
                                (rho * omega)
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index]

                if scattering_index[0] == 3:
                    # Polar optical phonon scattering

                    # r = 1, with corrections (not a simple energy relation)

                    # w0 = (sqrt(m)*e^2*omega)*(1/eps(inf)-1/eps(0))*
                    #      (n+1/2-/+1/2)/(4*pi*hbar*sqrt(2))

                    # energy dep = E^-1/2 * ln(sqrt(E)+sqrt(E +/- hbar omega)/
                    #              |sqrt(E)-sqrt(E +/- hbar omega)|)

                    # m = effective mass, units of m_e
                    # omega = optical phonon angular frequency [THz]
                    # eps(inf) = electronic permitivity, units of vacuum
                    #            permitivity
                    # eps = ionic permitivity, units of vacuum permitivity
                    # n = optical phonon occupation number

                    # ONLY VALID IF E > hbar omega (assume elastic process)

                    # ALSO, NEGLECTED SCREENING (ONLY VALID FOR HIGHLY DOPED)
                    if tr.bs.explicit_prefact[band][3]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][3]
                    else:
                        if tempi[0] == 0:
                            omega = tr.bs.omegao[band]
                            epsi = tr.bs.epsi[band]
                            eps = tr.bs.eps[band]
                            n = tr.bs.no[band]
                            if tr.bs.emi[band]:
                                emi = 1.0
                                sign = -1.0
                            temp = sign * scaling[3, 1] * omega
                            energy_correction_prefactor[
                                band, scattering_index] = temp
                            # set the transfer energy (used sometimes when
                            # integrating the explicit tau)
                            tr.tau_energy_trans[band, scattering_index] = temp
                            factor[band, scattering_index] = \
                                scaling[3, 0] * np.sqrt(effmass) * \
                                omega * (1 / epsi - 1 / eps) * \
                                (n + 1.0 * emi)
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index]

                if scattering_index[0] == 4:
                    # Piezoelectric acoustic phonon scattering

                    # r = 1, with corrections (not a simple energy relation)

                    # w0 = ((p^2*sqrt(m)*e^2*k*T)/ \
                    #       (sqrt(8)*pi*eps^2*hbar^2*rho*v^2))

                    # energy dep = E^-1/2 * [ln(1+4E/E_0)-1/(1+E_0/4E)]

                    # where E_0 = hbar^2 isl^2 / 2 m

                    # m = effective mass, unitless units of m_e
                    # p = piezoelectric constant [C/m^2]
                    # T = temperature [K]
                    # eps = electronic dielectric constant [F/m]
                    # rho = mass density [g/cm^3]
                    # v = speed of sound [m/s]
                    # isl = inverse screening length [AA^-1]

                    # INCLUDES BOTH ABSORPTION AND EMMISION
                    if tr.bs.explicit_prefact[band][4]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][4]
                    else:
                        if tempi[0] == 0:
                            p = tr.bs.p[band]
                            eps = tr.bs.eps[band]
                            rho = tr.bs.rho[band]
                            speed_sound = tr.bs.speed_sound[band]
                            isl = tr.bs.isl[band]
                            if effmass < constants.zero:
                                effmass = constants.zero
                            energy_correction_prefactor[band, scattering_index] = \
                                scaling[4, 1] * np.power(isl, 2.0) / effmass
                            denom = np.power(eps * speed_sound, 2.0) * rho
                            if denom < constants.zero:
                                denom = constants.zero
                            factor[band, scattering_index] = scaling[4, 0] * \
                                np.sqrt(effmass) * \
                                np.power(p, 2.0) / denom
                            # ignore overflow when performing multiply
                            np.seterr(over='ignore')
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index] * tempv

                if scattering_index[0] == 5:
                    # Ionized impurity scattering (BH)

                    # r = 2, with corrections (not a simple energy relation)

                    # w0 = Z^2*e^4*n_i/(32*pi*sqrt(2*m)*eps^2)

                    # energy dep = E^-3/2 * (ln(1+gamma)-gamma/(1+gamma))

                    # m = effective mass, units of m_e
                    # Z = number of charge units (of e) of the impurity
                    # n_i = ionized impurity density [10^21 cm^-3]
                    # eps = electronic dielectric constant [in units of epsilon_0]
                    # gamma = 4E/E_0
                    # E_0 = hbar^2 isli^2/2m
                    # isli = inverse screening length in AA^-1
                    if tr.bs.explicit_prefact[band][5]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][5]
                    else:
                        if tempi[0] == 0:
                            z = tr.bs.z[band]
                            n_i = tr.bs.ni[band]
                            eps = tr.bs.eps[band]
                            isli = tr.bs.isli[band]
                            e_0 = scaling[5, 1] * np.power(isli, 2.0) / effmass
                            denom = np.power(eps, 2.0) * np.sqrt(effmass)
                            if denom < constants.zero:
                                denom = constants.zero
                                factor[band, scattering_index] = \
                                    scaling[5, 0] * \
                                    np.power(z, 2.0) * n_i / denom
                            energy_correction_prefactor[band,
                                                        scattering_index] = e_0
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index]

                if scattering_index[0] == 6:
                    # Ionized impurity scattering (CW)

                    # r = 2, with corrections (not a simple energy relation)

                    # w0 = e^4*n_i/(32*pi*sqrt(2m)*eps^2)

                    # energy dep = E^-3/2 * ln(1+gamma)

                    # m = effective mass, units of m_e
                    # n_i = ionized impurity density [10^21 cm^-3]
                    # eps = electronic dielectric constant [in units of epsilon_0]
                    # gamma = (8pi b epsilon E/e^2)^2
                    # b = (3/4pi n_i)^1/3
                    if tr.bs.explicit_prefact[band][6]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][6]
                    else:
                        if tempi[0] == 0:
                            z = tr.bs.z[band]
                            n_i = tr.bs.ni[band]
                            eps = tr.bs.eps[band]
                            b = np.power(  # pylint: disable=assignment-from-no-return
                                4 * constants.pi * n_i / 3.0, -1.0 / 3.0)
                            factor[band, scattering_index] = scaling[
                                6, 0] * n_i / np.sqrt(effmass)
                            energy_correction_prefactor[band, scattering_index] = \
                                scaling[6, 1] * np.power(eps * b, 2.0)
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index]

                if scattering_index[0] == 7:
                    # Alloy scattering

                    # r = 0

                    # w0 = sqrt(2) * V * m^3/2 V_diff^2 x(1-x) / (hbar^4)

                    # energy dep = E^-1/2

                    # m = effective mass, units of m_e
                    # V = volume of unit cell in AA^3
                    # V_diff = difference between the two species atomic
                    #          potential in eV
                    # x = the fraction (of e.g. the compound A_xB_(1-x)C)

                    # NEUTRAL MODEL, OTHERWISE USE IONIZED IMPURITY
                    if tr.bs.explicit_prefact[band][7]:
                        prefactor_scattering[tempi,
                                             band,
                                             scattering_index] = \
                            tr.bs.explicit_prefact_values[band][7]
                    else:
                        if tempi[0] == 0:
                            v = tr.lattice.volume
                            v_diff = tr.bs.vdiff[band]
                            x = tr.bs.alloyconc[band]
                            factor[band, scattering_index] = \
                                scaling[7, 0] * np.power(effmass, 1.5) * \
                                np.power(v_diff, 2.0) * v * x * (1 - x)
                        prefactor_scattering[tempi, band, scattering_index] = \
                            factor[band, scattering_index]

    # now squeeze in the energy dependence (either directly or indirectly)

    # multiply with energy term given by the energy spacing
    # check if we have negative values in energies, take absolute value
    # also add shifts of energy from phonons, exitations etc from
    # correction_prefactor
    if not use_eonk:
        # use supplied energy grid
        with np.errstate(divide='ignore'):
            power_energy = np.power(  # pylint: disable=assignment-from-no-return
                np.abs(energies[np.newaxis, :, np.newaxis] +
                       energy_r_correction[:, np.newaxis, :]),
                r_factor_includinghalf[np.newaxis, np.newaxis, :])
    else:
        # use energies that lies on the k-point grid
        with np.errstate(divide='ignore'):
            power_energy = np.power(  # pylint: disable=assignment-from-no-return
                np.abs(energies[:, :, np.newaxis] +
                       energy_r_correction[:, np.newaxis, :]),
                r_factor_includinghalf[np.newaxis, np.newaxis, :])

    # in the following we apply energy corrections
    # remove nan values (force to zero)
    power_energy = np.nan_to_num(power_energy)
    # ignore overflow
    with np.errstate(over='ignore'):
        scattering[:, :, :, 0:num_scatterings - 2] = \
            power_energy[np.newaxis, :, :, 0:num_scatterings - 2] * \
            prefactor_scattering[:, :, np.newaxis, 0:num_scatterings - 2]
    # TODO: HERE WE SHOULD HAVE A REMOVE NAN AND NUM AGAIN? pylint: disable=fixme

    # calculate correction factor for tau[3] (polar optical) and add
    e_dep = np.sqrt(
        np.abs(energies[np.newaxis, :]) +
        energy_correction_prefactor[:, np.newaxis, 3])
    # remove nan values (force to zero)
    e_dep = np.nan_to_num(e_dep)
    denom = np.abs(np.sqrt(np.abs(energies[np.newaxis, :])) - e_dep)
    denom[denom < constants.zero] = constants.zero
    with np.errstate(divide='ignore'):
        e_cor_fact = np.log(  # pylint: disable=assignment-from-no-return
            (np.sqrt(np.abs(energies[np.newaxis, :])) + e_dep) / denom)
    # remove inf values (force to the largest supported value)
    e_cor_fact = np.nan_to_num(e_cor_fact)
    scattering[:, :, :, 3] = scattering[:, :, :, 3] * e_cor_fact[np.
                                                                 newaxis, :, :]

    # calculate correction factor for tau[4] (piezoelectric) and add
    denom = energy_correction_prefactor[:, np.newaxis, 4]
    denom[denom < constants.zero] = constants.zero
    with np.errstate(over="ignore"):
        e_dep = 4 * np.abs(energies[np.newaxis, :]) / denom
    # remove inf values (force to the largest supported value)
    e_dep = np.nan_to_num(e_dep)
    e_cor_fact = np.log(1 + e_dep) - e_dep / (e_dep + 1)
    # remove nan values (force to zero)
    # ignore overflow
    with np.errstate(over="ignore"):
        scattering[:, :, :,
                   4] = scattering[:, :, :, 4] * e_cor_fact[np.newaxis, :, :]

    # calculate the correction factor for tau[5] (BH ionized imp) and add
    # (very similar to piezoelectric, but we assume we can have different
    # screening lengths)
    denom = energy_correction_prefactor[:, np.newaxis, 5]
    denom[denom < constants.zero] = constants.zero
    with np.errstate(over="ignore"):
        e_dep = 4 * np.abs(energies[np.newaxis, :]) / denom
    # remove inf values (force to the largest supported value)
    e_dep = np.nan_to_num(e_dep)
    e_cor_fact = np.log(1 + e_dep) - e_dep / (e_dep + 1)
    scattering[:, :, :, 5] = scattering[:, :, :, 5] * e_cor_fact[np.
                                                                 newaxis, :, :]

    # calculate the correction factor for tau[6] (CW ionized imp) and add
    # (very similar to BH, but different correction factor)
    e_dep = np.power(np.abs(energies[np.newaxis, :]), 2.0) * \
        energy_correction_prefactor[:, np.newaxis, 6]
    # remove inf values (force to the largest supported value)
    e_dep = np.nan_to_num(e_dep)
    e_cor_fact = np.log(1 + e_dep)  # pylint: disable=assignment-from-no-return
    scattering[:, :, :, 6] = scattering[:, :, :, 6] * e_cor_fact[np.
                                                                 newaxis, :, :]

    # now add the constant scattering part given in fs units (so invert)
    prefactor_scattering[:, :, num_scatterings -
                         1] = 1.0 / tr.bs.tau0c[np.newaxis, :]
    scattering[:, :, :, num_scatterings - 1] = 1.0 / \
        tr.bs.tau0c[np.newaxis, :, np.newaxis]
    # set up array to force zeros into the sum array if one only wants certain
    # scattering mechanisms in the total sum
    iit = inc_in_total[:, np.newaxis, :] * \
        np.ones(num_energy_steps, dtype=int)[np.newaxis, :, np.newaxis]
    # now calculate the total scattering rate
    scattering_total = (np.nan_to_num(scattering) * iit).sum(-1)
    # and then, since up til now we have calculated the scattering rate
    # we invert all values in order to get the "tau" in fs units
    # also set true zero to a small value before inverting
    scattering[scattering < constants.zero] = constants.zero
    scattering_total[scattering_total < constants.zero] = constants.zero
    scattering = np.nan_to_num(scattering)
    # now invert the scattering array (use same name
    # to save memory)
    scattering = 1.0 / scattering
    with np.errstate(over="ignore"):
        scattering_total_inv = 1.0 / scattering_total
    # remove nan, inf etc.
    scattering = np.nan_to_num(scattering)
    scattering_total_inv = np.nan_to_num(scattering_total_inv)
    # tau0 for use in the closed Fermi integrals
    prefactor_scattering[
        prefactor_scattering < constants.zero] = constants.zero
    scattering_tau0 = np.nan_to_num(1.0 / prefactor_scattering)
    # now make sure the non selected scattering mechanisms contain very
    # large value (so that the scattering rate W=0)
    non_selected = np.array(1 - iit[:, 0, :], dtype=bool)
    scattering_tau0[:, non_selected] = constants.large
    return scattering, scattering_total_inv, scattering_tau0


def find_r_for_closed(tr, band):
    """
    Analyze the input tau0 and find the associated scattering values r.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    band : integer
        The band index.

    Returns
    -------
    integer
        Two times the r value to avoid half integer values.

    Notes
    -----
    These are necessary for the analytic Fermi integrals.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running find_r_for_closed.")

    # check that only one scattering mech is entered
    if np.sum(tr.bs.select_scattering[band]) > 1:
        logging.error("Parabolic Fermi integral routines are only "
                      "defined for one type of scattering. Typically, "
                      "set one scattering mechnisms to a value (tau0) "
                      "and the other factors to zero. Exiting.")
        sys.exit(1)
    # return value and multiply by two
    return int(2 * tr.scattering_r_factor[tr.bs.select_scattering[band]])


def combined_scattering(tr, energy, tau0, energy_trans):
    r"""
    Calculates the total relaxation time.

    Parameters
    ----------
    tr : object
        A `Transport()` object
    energy : float
        The energy of the charge carrier in eV.
    tau0 : ndarray
        | Dimension: (12)

        Contains the relaxation time prefactors for the
        different scattering mechanisms in units of fs.
    energy_trans : ndarray
        | Dimension: (12)

        Contains the energy transitions in eV (that is added to the energy
        in :math:`\\tau=\\tau_0E^{r-1/2}`, typically,
        :math:`E=E+\\hbar \\omega`, where :math:`\\hbar \\omega`
        is the size of the energy transition. Set it to zero for
        the non-relevant scattering mechanisms.
    effmass : float
        The effective mass in units of the electron mass

    Returns
    -------
    float
        The combined relaxation time in fs.

    Notes
    -----
    Calculates the total relaxation time

    .. math:: \\frac{1}{\\tau}=\\sum_i \\frac{1}{\\tau_i},

    where :math:`\\tau=\\tau_0E^{r-1/2}`.
    The array `scattering_tau0_select` determines which
    scattering to include in the sum. Consult :func:`scattering_parabolic`
    for additional details. The scattering prefactors
    :math:`\\tau_0` are ordered in a sequence described there.
    The `scattering_tau0_select` follows this sequence and is
    set in the bandstructure configuration file.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running combined_scattering.")

    tau_decomp = (
        tau0 * np.power(energy + energy_trans, 0.5 - tr.scattering_r_factor)
    )[tr.scattering_tau0_select]
    tau_decomp = np.nan_to_num(tau_decomp)
    return np.nan_to_num(1.0 / np.sum(tau_decomp))


def interpolate(tr, method="linear"):  # pylint: disable=too-many-locals
    """
    Interpolates the scattering array on all available energies.

    Parameters
    ----------
    tr : object
        A `Transport()` object containing the scattering arrays
        and the energies etc.
    method : string, optional
        The interpolation method to use. Uses
        the :func:`interp1d` function of Scipy and this sets the
        parameter `kind`. Defaults to "linear".

    Returns
    -------
    None

    See Also
    --------
    scipy.interpolate.interp1d

    Notes
    -----
    Here we only perform an interpolation on the array containing the total
    relaxation time since this is used during the transport calculations.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running interpolate.")

    logger.info("Interpolating the scattering values on the energygrid "
                "of the band structure.")

    # check the validity of the method parameter
    if method not in constants.interp1d_methods:
        logger.error("The 'method' parameter passed is not a present "
                     "option for your interp1d version. Exiting")
        sys.exit(1)

    energies = tr.scattering_energies
    scattering_total_inv = tr.scattering_total_inv
    scattering_inv = tr.scattering_inv
    inter_energies = tr.bs.energies

    # We need to make sure that the requested energies are
    # inside the bounds of the original data set. We do this
    # simply by padding the smallest and largest values by a
    # small constant. These values should anyway be far outside
    # the important energy region anyway.
    max_inter_energies = np.amax(inter_energies)
    min_inter_energies = np.amin(inter_energies)
    # need to check if other values are similar to some constant
    # first, check largest
    replace = np.where(
        inter_energies > (max_inter_energies - constants.zerocut))

    inter_energies[replace] = inter_energies[replace] - \
        constants.zerocut

    # then smallest
    replace = np.where(
        inter_energies < (min_inter_energies + constants.zerocut))
    inter_energies[replace] = inter_energies[replace] + \
        constants.zerocut

    # this is really dirty, but since the number of temperatures
    # and bands is usually pretty limited it was not a top priority to
    # increase the speed of this one, could easily be done
    # in a different way...
    num_temp_steps = scattering_total_inv.shape[0]
    num_bands = scattering_total_inv.shape[1]
    num_energies = inter_energies.shape[1]
    num_scatterings = scattering_inv.shape[3]
    scattering_total_inv_inter = np.zeros((num_temp_steps, num_bands,
                                           num_energies))
    if not tr.param.onlytotalrate:
        scattering_inv_inter = np.zeros((num_temp_steps, num_bands,
                                         num_energies, num_scatterings))
    for temp in range(num_temp_steps):
        for band in range(num_bands):
            # do the interpolation of the total first
            inter_total_inv = \
                scipy.interpolate.interp1d(energies,
                                           scattering_total_inv[temp, band])
            scattering_total_inv_inter[temp, band] = inter_total_inv(
                inter_energies[band])

            if not tr.param.onlytotalrate:
                # and then the decomposed, but only for the selected mechanisms
                include_scattering = np.nonzero(
                    tr.bs.select_scattering[band])[0]
                exclude_scattering = np.nonzero(
                    1 - tr.bs.select_scattering[band])[0]
                for scattering in include_scattering:
                    inter_inv = scipy.interpolate.interp1d(
                        energies, scattering_inv[temp, band, :, scattering])
                    scattering_inv_inter[
                        temp, band, :, scattering] = \
                        inter_inv(inter_energies[band])
                for scattering in exclude_scattering:
                    scattering_inv_inter[temp, band, :, scattering] = np.full(
                        num_energies, constants.large)

    # store new dense scattering arrays, also make sure we have
    # no nans or inf
    tr.scattering_total_inv = np.nan_to_num(scattering_total_inv_inter)
    if not tr.param.onlytotalrate:
        tr.scattering_inv = np.nan_to_num(scattering_inv_inter)
    # now the scattering energies are in fact just the energies
    # for each band and kpoint, setting it like this
    # should avoid occupying more memory than necessary,
    # but could be a bugsource...beware!
    tr.scattering_energies = inter_energies


def pad_scattering_values(tr):
    """
    Pad the scattering values.

    Parameters
    ----------
    tr : object
        A `Transport()` object.

    Returns
    -------
    None

    Notes
    -----
    The padded values are stored in the `tr` object.

    We need to pad the energies where the dos is calculated
    with a larger number of samples such that we cover the whole
    energy range in the stored bandstructure due to later
    interpolation routines etc. not going out of bounds
    when such energies are passed to the interpolator etc.

    We can set it to a large number because 1 eV outside the
    chemical potential. We already know no states contribute
    at temperatures below 2000 K.

    """
    # fetch min and max of the energies stored for the
    # bandstructure
    scattering_inv = tr.scattering_inv
    scattering_total_inv = tr.scattering_total_inv
    energies = tr.scattering_energies
    emin, emax = tr.bs.fetch_min_max_energy()
    scattering_emin = energies[0]
    scattering_emax = energies[energies.shape[0] - 1]
    estep = energies[1] - energies[0]
    if emin < scattering_emin:
        # fetch missing interval below and pad with linear ramp
        ebelow = scattering_emin - emin
        numsteps_below = int(np.ceil(ebelow / estep))
        emin = scattering_emin - numsteps_below * estep
        energies = np.pad(
            energies, (numsteps_below, 0), 'linear_ramp', end_values=(emin, 0))
        # now pad scattering arrays with endvalues below
        scattering_inv = np.pad(scattering_inv, ((0, 0), (0, 0),
                                                 (numsteps_below, 0), (0, 0)),
                                'edge')
        scattering_total_inv = np.pad(scattering_total_inv,
                                      ((0, 0), (0, 0),
                                       (numsteps_below, 0)), 'edge')
    if emax > scattering_emax:
        # fetch missing interval above and pad with linear
        # tramp
        eabove = emax - scattering_emax
        numsteps_above = int(np.ceil(eabove / estep))
        energies = np.pad(
            energies, (0, numsteps_above), 'linear_ramp', end_values=(0, emax))
        # now pad scattering arrays with endvalues above
        scattering_inv = np.pad(scattering_inv, ((0, 0), (0, 0),
                                                 (0, numsteps_above), (0, 0)),
                                'edge')
        scattering_total_inv = np.pad(scattering_total_inv,
                                      ((0, 0), (0, 0),
                                       (0, numsteps_above)), 'edge')
    tr.scattering_inv = scattering_inv
    tr.scattering_total_inv = scattering_total_inv
    tr.scattering_energies = energies


def check_scattering(tr):
    """
    Checks the scattering arrays.

    Also that they are dimensionalized to the energy values stored in the current `Bandstructure()`
    object.

    Parameters
    ----------
    tr : object
        A `Transport()` object.

    Returns
    -------
    None

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running check_scattering.")

    try:
        tr.scattering_total_inv
    except AttributeError:
        logger.error("Could not find 'scattering_total_inv' in the "
                     "current 'Transport()' object. Exiting.")
        sys.exit(1)

    if not tr.scattering_total_inv.shape[1] == tr.bs.energies.shape[0]:
        logger.error("The array 'scattering_total_inv' does not contain "
                     "the same number of bands as present in the current "
                     "band structure. Exiting.")
        sys.exit(1)
    if not tr.scattering_total_inv.shape[2] == tr.bs.energies.shape[1]:
        logger.error("The array 'scattering_total_inv' does not contain "
                     "the same number of kpoints (energies) as the current "
                     "band structure. Exiting.")
        sys.exit(1)
