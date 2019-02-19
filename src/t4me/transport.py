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
"""Contains routines to set up the calculation of the charge carrier transport coefficients."""

# pylint: disable=useless-import-alias, too-many-arguments, invalid-name, too-many-statements, too-many-lines, global-statement

import sys
import logging
import numpy as np
import scipy

import t4me.scattering as scattering
import t4me.lbtecoeff as lbtecoeff
import t4me.constants as constants


class Transport():  # pylint: disable=too-many-instance-attributes, too-many-branches
    """
    Involves all transport related routines.

    Parameters
    ----------
    bs : object
       A `Band()` object.
    lattice : object
       A `Lattice()` object.
    param : object
       A `Param()` object.

    """

    def __init__(self, bs):
        self.bs = bs
        self.param = bs.param
        self.lattice = bs.lattice
        self.temperatures = self.fetch_temperatures()
        self.fetch_chempots()
        # check which bands to be included in the transport integration
        self.fetch_relevant_bands()
        # on init also setup scattering
        self.setup_scattering()
        self.sigma = None
        self.seebeck = None
        self.lorenz = None
        self.hall = None
        self.ccn = None
        self.ccp = None

    def setup_scattering(self,
                         dos=None,
                         dos_energies=None,
                         select_scattering=None):
        """
        Selects which how to set up the carrier scattering.

        Parameters
        ----------
        dos : ndarray
            | Dimension: (N,M)

            Array containing the partial density of states in units
            of 1/eV/AA^3, where N is the band index and M is the
            energy index.
        dos_energies : ndarray
            | Dimension: (M)

            Array containing the energy in eV at M samplings where the
            density of states is calculated.
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
        None

        See Also
        --------
        scattering.scattering_dos
        scattering.scattering_parabolic

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running setup_scattering.")

        if select_scattering is None:
            select_scattering = self.bs.select_scattering

        # numerical scattering models
        if not self.param.transport_use_analytic_scattering:
            if self.param.transport_method == "closed":
                logger.error("You cannot use density of based scattering "
                             "together with the closed analytical expression "
                             "for the transport coefficients. Exiting.")
                sys.exit(1)
            logger.info("Setting up the scattering models "
                        "based on the density of states.")
            if dos is None:
                if self.bs.dos_partial is None:
                    # need to calculate the density of states
                    logger.info("No density of states present to set up "
                                "the scattering properties. Calculating it.")
                    self.bs.calc_density_of_states(transport=True)
                else:
                    # check that the energy range is at least
                    # transport_energycutband outside the min
                    # and max of the chemical potential
                    if not self.bs.check_dos_energy_range(
                            self.param.transport_chempot_min -
                            self.param.transport_energycutband,
                            self.param.transport_chempot_max +
                            self.param.transport_energycutband):
                        logger.info("Recalculating the density of states due "
                                    "to a too narrow energy range in order "
                                    "to set up the scattering data in a "
                                    "sufficiently wide energy range.")
                        self.bs.calc_density_of_states(transport=True)

            dos = self.bs.dos_partial
            if dos_energies is None:
                energies = self.bs.dos_energies
            else:
                energies = dos_energies
            scattering_inv, scattering_total_inv, scattering_tau0 = \
                scattering.scattering_dos(self,
                                          dos,
                                          energies,
                                          select_scattering)

            # check if we want Cubature/GeometricTools that the ontfly flag
            # is enabled for the scattering
            if self.param.transport_integration_method == "cubature":
                if not self.param.transport_use_scattering_ontfly:
                    logger.info("User wants to use Cubature/GeometricTools "
                                "and with non-analytical scattering. The flag "
                                "'transport_use_scattering_ontfly' is now "
                                "forced to True.")
                    self.param.transport_use_scattering_ontfly = True

            # store commons
            self.scattering_inv = scattering_inv
            self.scattering_total_inv = scattering_total_inv
            self.scattering_tau0 = scattering_tau0
            self.scattering_energies = energies

            # now pad the values such that requests from e.g. an
            # interpolation routine does not get out of bounds
            # (dos energy samples sometimes are inside the energy
            # (range of the eigenvalues)
            scattering.pad_scattering_values(self)

            # now that this is done we have to interpolate the scattering
            # such that we have a value for each k-point entry...doing this
            # makes it possible to perform the static k-point integrals later with
            # appropriate scattering included
            # do not do this if transport_use_scattering_ontfly is set to True
            # as we then extract these values during integration
            if not self.param.transport_use_scattering_ontfly and \
               self.param.transport_method == "numerick":
                scattering.interpolate(self)

        # analytic scattering models
        else:
            logger.info("Setting up the scattering models "
                        "based on the parabolic band models.")
            if self.param.transport_integration_method == "cubature":
                # for on the fly interpolation we need a linear
                # grid for the scattering and that is it
                use_eonk = False
                # check that the energy range is at least
                # transport_energycutband outside the min
                # and max of the chemical potential
                emin = (self.param.transport_chempot_min -
                        self.param.transport_energycutband)
                emax = (self.param.transport_chempot_max +
                        self.param.transport_energycutband)
                # does the analytic dos actually exist?
                if not self.bs.check_dos_energy_range(emin, emax):
                    logger.info("Adjusting the energy range locally in "
                                "order to generate the scattering data "
                                "in a sufficiently wide energy range. "
                                "Still using 'dos_num_samples' for the "
                                "spacing.")
                else:
                    emin = self.param.dos_e_min
                    emax = self.param.dos_e_max
                # generate grid of energies, still using
                # the samples requested for the density of states
                energies = np.linspace(emin, emax, self.param.dos_num_samples)
            else:
                use_eonk = True
                energies = self.bs.energies
            scattering_inv, scattering_total_inv, scattering_tau0 = \
                scattering.scattering_parabolic(self,
                                                energies,
                                                select_scattering,
                                                use_eonk=use_eonk)

            # store commons
            self.scattering_inv = scattering_inv
            self.scattering_total_inv = scattering_total_inv
            self.scattering_tau0 = scattering_tau0
            self.scattering_energies = energies

            if self.param.transport_integration_method == "cubature":
                # now pad the values such that requests from e.g. an
                # interpolation routine does not get out of bounds
                # (dos energy samples sometimes are inside the energy
                # (range of the eigenvalues)
                scattering.pad_scattering_values(self)

    def calc_transport_tensors(  # pylint: disable=too-many-locals # noqa: MC0001
            self,
            bs=None,
            temperatures=None,
            chempots=None,
            method=None):
        r"""
        Selects which method to use when calculating the transport coefficients.

        Parameters
        ----------
        bs : A `Band()` object containing the band structure.
        temperatures : ndarray, optional
            | Dimension: (N)

            Contains N different temperatures in K. If not supplied the
            `temperature` from the active `Transport()` object is used.
        chempots : ndarray, optional
            | Dimension: (M)

            Contains M different chemical potentials in eV. If not
            supplied the `chempot` from the active `Transport()` object
            is used.
        method : {"closed", "numeric", "numerick"}
            If `method` is not supplied is defaults to "numeric" unless
            bandstructure data is read numerically or generated (all
            cases where the closed Fermi Dirac integrals cannot be
            used) when it defaults to "numerick".

            | "closed" evaluates the closed Fermi integrals where only
            | one scattering mechanism is possible per band. Only valid
            | for systems where one can strictly rely on a parametrized
            | parabolic bandstructure based on effective mass models.
            | Parameters (e.g. effective masses for each band) are set
            | in the bandstructure configuration file.
            | The driver routine is :func:`lbtecoeff.parabolic_closed`

            | "numeric" similar to "closed, but evaluates the Fermi
            | integrals in an open form (e.g. it is possible to
            | concatenate the scattering mechanisms, which is not
            | possible for the closed Fermi integrals).
            | The driver routine is :func:`lbtecoeff.parabolic_numeric`

            | "numerick" evaluates the transport integrals more generally
            | as an integral over the k-points. It is less restrictive
            | than the two other options, but also more prone to
            | convergence issues etc. However, for bandstructures
            | read from datafiles, this is the only option.
            | The driver routine is :func:`lbtecoeff.numerick`

        Returns
        -------
        sigma, seebeck, lorenz : ndarray, ndarray, ndarray
            | Dimension: (N,M,3,3), (N,M,3,3), (N,M,3,3)

            Returns the electrical condcutivity, Seebeck coefficient and
            Lorenz tensor for N temperature and M chemical potential
            steps in units of :math:`\\mathrm{S}/\\mathrm{m}`,
            :math:`\\mu \\mathrm{V}/\\mathrm{K}`,
            :math:`\\mathrm{V^2}/\\mathrm{K^2}`. These are stored in the
            current `Transport()` object.

        See Also
        --------
        lbtecoeff.parabolic_closed
        lbtecoeff.parabolic_numeric
        lbtecoeff.numerick

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running calc_transport_tensors.")

        # first check that the scattering properties if now, give user
        # error
        try:
            self.scattering_tau0
        except AttributeError:
            logger.error("The user wants to calculate the transport "
                         "tensors, but no scattering mechanisms have "
                         "been configured. Please make an object of "
                         "the Transport() class, which setups up the "
                         "scattering from the param.yml and "
                         "bandparam.yml files before caling the routine "
                         "that does the transport tensor calculations. "
                         "Exiting.")
            sys.exit(1)

        # set defaults
        if temperatures is None:
            temperatures = self.temperatures
        if chempots is None:
            chempots = self.chempots
        if method is None:
            method = self.param.transport_method

        # figure out the current configuration if VASP input, for sure,
        # we do not have analytick models
        numerick = False
        if self.param.read == "vasp" or self.param.read[:5] == "numpy" \
           or self.param.read == "w90":
            numerick = True
            # now check if the user have set transport_method to closed
            # and give warning
            if self.param.transport_method == "closed":
                logger.error("The user requests to read numerical data and "
                             "solve the transport integrals using the closed "
                             "Fermi-Dirac integrals. This is not possible. "
                             "User, please make up your mind. Exiting. ")
                sys.exit(1)
        else:
            # or if any band is generated with band type different from 0
            if np.any(self.bs.bandparams[:, 0] != 0):
                numerick = True
        # now check if user wants numerick anyway
        if method == "numerick" or self.param.transport_method == "numerick":
            numerick = True

        if method != "numerick" and numerick:
            logger.info(
                "User requested to use the method '%s' "
                "for integration, but "
                "at the same time wants to read data from "
                "VASP, NumPy, Wannier90 or have generated "
                "non-parabolic bands."
                "This is not possible and we now set the method "
                "to 'numerick'.", method)

        if bs is None:
            bs = self.bs

        # analytick expansions of energy etc.
        if not numerick:
            sigma = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
            seebeck = np.zeros((temperatures.shape[0], chempots.shape[0], 3,
                                3))
            lorenz = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
            hall = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
            ccn = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
            ccp = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
            # set transport tensor scheme for analytick in python
            # which is slow, so we can chose only to calculate certain elements
            # TODO: incorporate this for all methods pylint: disable=fixme
            # loop temperaturs
            for indext, temp in np.ndenumerate(temperatures):
                # fetch eta for each band
                etas = self.fetch_etas(chempots, temp).T
                # fetch tau0 for a given temperature
                tau0 = self.scattering_tau0[indext]
                # loop temperature and calculate closed integrals
                for indexe in range(chempots.shape[0]):
                    sigma_tensor, seebeck_tensor, lorenz_tensor, \
                        hall_tensor, cc_tensor_n, cc_tensor_p = lbtecoeff.parabolice(self, etas[indexe], temp,
                                                                                     bs, tau0, method)
                    sigma[indext, indexe] = sigma_tensor
                    seebeck[indext, indexe] = seebeck_tensor
                    lorenz[indext, indexe] = lorenz_tensor
                    hall[indext, indexe] = hall_tensor
                    ccn[indext, indexe] = cc_tensor_n
                    ccp[indext, indexe] = cc_tensor_p
        # fully numerick evaluation, for the purpose of speed, the loop
        # over temperature and chemical potential is done internally.
        # The return is (temperature,chempot,3,3) arrays
        else:
            sigma, seebeck, lorenz = lbtecoeff.numerick(
                self, chempots, temperatures, bs)
            # TODO: FIX THE HALL TENSOR ASAP (MAYBE ALSO THE NERNST) pylint: disable=fixme
            hall = sigma
            ccn = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))
            ccp = np.zeros((temperatures.shape[0], chempots.shape[0], 3, 3))

        # calculate the carrier concentration
        if numerick:
            # calculate the carrier concentration
            # check if dos exists
            if ((self.bs.dos_partial is None)
                    or (not self.param.carrier_dos_analytick)):
                self.bs.calc_density_of_states()
            # loop temperatures
            for indext, temperature in np.ndenumerate(temperatures):
                # loop chempots
                for indexe, chempot in np.ndenumerate(chempots):
                    ptype, ntype, _ = self.calc_carrier_concentration(
                        temperature, chempot)
                    ccp[indext, indexe, 0, 0] = ptype
                    ccn[indext, indexe, 0, 0] = ntype
        self.sigma = sigma
        self.seebeck = seebeck
        self.lorenz = lorenz
        self.hall = hall
        self.ccn = ccn
        self.ccp = ccp

    def fetch_relevant_bands(self, tr=None):
        """
        Locate bands that will be included in the transport integrals.

        Parameters
        ----------
        tr : object, optional
            A `Transport()` object.

        Returns
        -------
        None

        Notes
        -----
        The included bands are located by considering the input
        range of chemical potentials from `transport_chempot_min`
        and `transport_chempot_max` padded with the value
        `transport_energycutband` on
        each side (see the general configuration file).

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running fetch_relevant_bands.")

        if tr is None:
            energies = self.bs.energies
            param = self.param
        else:
            energies = tr.bs.energies
            param = tr.param

        # check if user supplied specific bands for calculation and
        # let them know if they have supplied this (easy to make mistakes)
        if param.transport_include_bands:
            # first check that we actually have all the necessary bands
            band_index = np.amax(param.transport_include_bands)
            if band_index > energies.shape[0]:
                logger.error("User requested a band that is not included in "
                             "the original dataset. Exiting.")
                sys.exit(1)
            logger.info(
                "User supplied specific bands so we are only performing "
                "transport calculation on those.")
            # shift index to zero
            transport_included_bands = [
                x - 1 for x in param.transport_include_bands
            ]

        else:
            e_min = param.transport_chempot_min - param.transport_energycutband
            e_max = param.transport_chempot_max + param.transport_energycutband
            transport_included_bands = []
            # loop bands, later add vectorize on band as well
            for band in range(energies.shape[0]):
                if energies[band][(energies[band] > e_min) &
                                  (energies[band] < e_max)].size != 0:
                    transport_included_bands.append(band)

        if tr is None:
            self.included_bands = np.array(
                transport_included_bands, dtype='intc')
        else:
            tr.included_bands = np.array(
                transport_included_bands, dtype='intc')

    def calc_carrier_concentration(  # pylint: disable=too-many-locals
            self,
            temperature,
            chempot,
            dos=None,
            dos_energies=None,
            band_decomp=False,
            defect_ionization=False):
        r"""
        Returns the charge carrier concentration.

        Parameters
        ----------
        temperature : float
            The temperature in K.
        chempot : float
            The chemical potential in eV.
        dos : ndarray, optional
            | Dimension: (N,M)

            Contains the band decomposed density of states for each
            band N and energy M. If not supplied, set to the `dos_partial`
            parameter of the current `Bandstructure()` object.
        dos_energies : ndarray, optional
            | Dimension: (M)

            The energies in eV where the density of states are sampled.
        band_decomp : boolean
            Return a band decomposed carrier concentration or not.
        defect_ionization : boolean
            Selects if defect ionization compensation should be
            included. The `donor_number`, `donor_energy`,
            `donor_degen_fact`, `acceptor_number`, `acceptor_energy`
            and `acceptor_degen_fact` need to be set in the
            general configuration file.

        Returns
        -------
        n_type : ndarray
            | Dimension: (N)

            Contains the n-type carrier concentration for each band
            index N in units of :math:`10^{21} \mathrm{cm}^{-3}`.
        p_type : ndarray
            | Dimension: (N)

            Contains the p-type carrier concentration for each band
            index N in units of :math:`10^{21} \mathrm{cm}^{-3}`.

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running calc_carrier_concentration.")

        if dos is None:
            dos = self.bs.dos_partial
        if dos_energies is None:
            dos_energies = self.bs.dos_energies

        num_bands = self.bs.bandparams.shape[0]
        n_type = np.zeros(num_bands)
        p_type = np.zeros(num_bands)
        ntype_index = np.where(
            dos_energies > self.param.carrier_conduction_energy)
        ptype_index = np.where(
            dos_energies < self.param.carrier_valence_energy)
        dos_energies_ntype = dos_energies[ntype_index]
        dos_energies_ptype = dos_energies[ptype_index]
        intrinsic = np.zeros(num_bands)
        beta = 1e5 / (constants.kb * temperature)
        for band in range(num_bands):
            if dos_energies_ntype.size > 0:
                # n-type, use only energies from carrier_conduction_energy
                # to the end of the array set in param.yml, slice
                integrand = dos[band][ntype_index] * \
                    fermi_dist(dos_energies_ntype, chempot, beta)

                n_type[band] = scipy.integrate.trapz(integrand,
                                                     dos_energies_ntype)
            # p-type, use only energies from start of array to
            # carrier_valence_energy set in param.yml, slice
            if dos_energies_ptype.size > 0:
                integrand = dos[band][ptype_index] * \
                    fermi_dist(-dos_energies_ptype, -chempot, beta)
                p_type[band] = scipy.integrate.trapz(integrand,
                                                     dos_energies_ptype)
        # make sure units of carrier concentration is 10^21 cm^-3
        n_type = 1e3 * n_type
        p_type = 1e3 * p_type
        # calculte intrinsic^2 (sum for each band first)
        intrinsic = np.multiply(n_type.sum(-1), p_type.sum(-1))
        if defect_ionization:
            donor_number = self.param.donor_number
            donor_degen_fact = self.param.donor_degen_fact
            donor_energy = self.param.donor_energy
            acceptor_number = self.param.acceptor_number
            acceptor_degen_fact = self.param.acceptor_degen_fact
            acceptor_energy = self.param.acceptor_energy
            donor_ion_number = donor_ionization(
                donor_number, donor_energy, donor_degen_fact, chempot, beta)
            acceptor_ion_number = acceptor_ionization(
                acceptor_number, acceptor_energy, acceptor_degen_fact, chempot,
                beta)
            n_type = 0.5 * (donor_ion_number - acceptor_ion_number) + \
                np.sqrt(np.power(0.5 * (donor_ion_number -
                                        acceptor_ion_number), 2.0)
                        + intrinsic)
            p_type = 0.5 * (acceptor_ion_number - donor_ion_number) + \
                np.sqrt(np.power(0.5 * (acceptor_ion_number -
                                        donor_ion_number), 2.0) +
                        intrinsic)
        if not band_decomp:
            p_type = p_type.sum(-1)
            n_type = n_type.sum(-1)
        return p_type, n_type, np.sqrt(intrinsic)

    def fetch_temperatures(self, store=True):
        """
        Set up the temperatures.

        Parameters
        ----------
        store : boolean, optional
            If given and set to True, the temperature array is in
            addition to being returned also stored in the active
            `Transport()` object.

        Returns
        -------
        temperature : (N) ndarray
             Contains N temperature linear samplings in units of K. The
             parameters `temperature_min`, `temperature_max` and
             `temperature_steps` in param.yml set the maximum and
             minimum temperature and its number of steps.

        """
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running fetch_temperatures.")

        temperature = np.linspace(self.param.temperature_min,
                                  self.param.temperature_max,
                                  self.param.temperature_steps)
        if store:
            self.temperature = temperature
            return temperature

        return temperature

    def fetch_chempots(self, store=True):
        """
        Set up the chemical potential.

        Parameters
        ----------
        store : boolean, optional
            If given and set to True, the chempot array is in addition
            to being returned also stored in the current `Transport()`
            object.

        Returns
        -------
        chempot : ndarray
             | Dimension: (N)

             Contains N chemical potential linear samplings in units of
             eV. The parameters `transport_chempot_min`,
             `transport_chempot_max` and `transport_chempot_samples` in
             param.yml set the maximum and minimum chemical potential
             and its number of samples.

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running fetch_chempots.")

        chempots = np.linspace(self.param.transport_chempot_min,
                               self.param.transport_chempot_max,
                               self.param.transport_chempot_samples)
        if store:
            self.chempots = chempots
            return chempots

        return chempots

    def fetch_etas(self, chempot, temperature):
        """
        Calculate the reduced chemical potential

        Parameters
        ----------
        chempot : ndarray
            | Dimension: (N)

            Contains N samples of the chemical potential in
            units of eV.

        temperature : float
            The temperature in K.

        Returns
        -------
        eta : ndarray
            | Dimension: (N)

            Contains N samples of the reduced chemical potential

        """

        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.debug("Running fetch_etas.")

        # convert to eta, shift and loop
        eta = np.zeros((self.bs.e0.shape[0], chempot.shape[0]))
        # valence bands: eta=e_shift-chempot
        bandtype = np.where(self.bs.status == 'v')
        eta[bandtype] = 1e5 * (self.bs.e0[bandtype, np.newaxis] -
                               chempot[np.newaxis, :]) / \
            (constants.kb * temperature)
        band = bandtype[0].shape[0]
        # conduction bands: eta=chempot-e_shift
        bandtype = np.where(self.bs.status == 'c')
        eta[bandtype] = 1e5 * (chempot[np.newaxis, :] -
                               self.bs.e0[bandtype, np.newaxis]) / \
            (constants.kb * temperature)

        band += bandtype[0].shape[0]
        # check that there are no funny bands not marked as v or c
        if band != self.bs.e0.shape[0]:
            logger.error("Some bands are not marked as a conduction or \
            valence band. Please correct input files. Exiting.")
            sys.exit(1)
        return eta


def fetch_chempot_from_etas(temperature, etas):
    r"""
    Calculate the chemical potential from eta and the temperature.

    Parameters
    ----------
    temperature : float
        The temperature in K.
    etas : ndarray
        | Dimension: N

        The unitless chemical potential, :math:`\\eta` for N
        steps.

    Returns
    -------
    chempots : ndarray
        | Dimension: N

        The chemical potentials in units of eV.

    """

    chempots = etas * constants.kb * 1e-5 * temperature

    return chempots


def donor_ionization(number, energy, degen, e_fermi, beta):
    """
    Returns the number of ionized donors.

    Parameters
    ----------
    number : float
        Number of donors.
    energy : float
        The energy in eV where the
        donor compensation is to be
        evaluated.
    degen : float
        The donor degeneration number.
    e_fermi : float
        The Fermi level in eV.
    beta : float The beta (1/kT) factor in eV.

    Returns
    -------
    float
        The donor ionization compensation.

    """

    return number / (1 + np.exp((energy - e_fermi) * beta) / degen)


def acceptor_ionization(number, energy, degen, e_fermi, beta):
    """
    Returns the number of ionized acceptors.

    Parameters
    ----------
    number : float
        Number of acceptors.
    energy : float
        The energy in eV where
        the acceptor compensation is to be evaluated.
    degen : float
        The acceptor degeneration number.
    e_fermi : float
        The Fermi level in eV.
    beta : float
        The beta (1/kT) factor in eV.

    Returns
    -------
    float
        The acceptor ionization compensation.

    """

    return number / (1 + np.exp((e_fermi - energy) * beta) / degen)


def fermi_dist(e, e_fermi, beta):
    """
    Returns the Fermi Dirac distribution function (without spin degeneracy).

    Parameters
    ----------
    e : float
        The energy in eV where the Fermi Dirac distribution is to
        be evaluated.
    e_fermi : float
        The Fermi level in eV.
    beta : float
        The beta factor (1/kT) in eV.

    Returns
    -------
    float
        The value of the Fermi function (without spin degeneracy).

    """

    return 1.0 / (1.0 + np.exp((e - e_fermi) * beta))
