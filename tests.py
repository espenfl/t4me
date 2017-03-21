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

#!/usr/local/bin/python
# python specifics
import sys
import os
import unittest
import numpy as np
import pickle
import logging
import types
import inspect

# locals
import constants
import utils
import bandstructure
import lattice
import transport
import inputoutput


def run_tests(tests="slow"):
    """
    The main driver routine for the tests.

    Parameters
    ----------
    tests : optional
        The tests to run. Defaults to "slow". One can also
        run with "fast" which excludes the time consuming tests.
        Individual tests
        can be executed. Consult the source code to find the
        specific function name that needs to be passed.

    Returns
    -------
    None

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)

    logger.info("Running test routines.")

    if tests == "slow" or tests == "fast" or tests is True:
        fast = unittest.TestLoader().loadTestsFromTestCase(
            TestTransportIntegralsFast)
        unittest.TextTestRunner(verbosity=2).run(fast)
        if tests == "slow":
            # also do slow tests
            slow = unittest.TestLoader().loadTestsFromTestCase(
                TestTransportIntegralsSlow)
            unittest.TextTestRunner(verbosity=2).run(slow)
    else:
        test = unittest.TestSuite()
        # try to find the method in all classes present
        # in this module
        g = globals().copy()
        test_found = False
        for name, object in g.items():
            if inspect.isclass(object):
                try:
                    test.addTest(object(methodName=tests))
                    test_found = True
                    break
                except:
                    pass
        if not test_found:
            logger.error("Did not find test: " + tests + " in any of "
                         "test classes. Exiting.")
            sys.exit(1)
        runner = unittest.TextTestRunner()
        runner.run(test)


class TestTransportIntegralsFast(unittest.TestCase):

    def test_closed_spherical(self):
        """
        Test the Seebeck and Lorenz coefficient, carrier
        concentration and the small-h Hall factor against
        the data published by May :cite:`rowe_2012_mpacit_may`
        for a selected set of unitless chemical potentials at
        300 K for a band with the effective mass of a free
        electron. Acoustic phonon scattering is used.

        """
        test_number = 1
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now fetch data we want to compare to
        may_data = self.may_data()

        # set up transport
        tran = transport.Transport(bs)

        # fetch chempots from etas at 300 K (may_data)
        temperature = 300
        tran.temperatures = np.array([temperature])
        tran.chempots = tran.fetch_chempot_from_etas(
            temperature, may_data[:, 0])

        # now calculate the coefficients
        tran.calc_transport_tensors(method="closed")
        seebeck = tran.seebeck.squeeze()[:, 0, 0]
        lorenz = tran.lorenz.squeeze()[:, 0, 0]
        cc = tran.ccp.squeeze()[:, 0, 0]
        hall = tran.hall.squeeze()[:, 0, 0]

        # calculate small r factor (100 for 10^21 to 10^19
        # conversion
        small_hall = 100 * constants.elcharge * cc * hall

        # stack so that the shape is similar to may_data
        this_data = np.column_stack(
            (seebeck.T, lorenz.T, cc.T, small_hall.T))

        # now calculate the relative difference
        difference = np.abs(
            (this_data - may_data[:, 1:5]) / may_data[:, 1:5])

        # print difference
        # now test (only two decimals)
        self.assertTrue(np.all(difference < 1e-2))

    def test_numeric_spherical_one_band(self):
        """
        Test the numerical solution of the energy
        dependent integral that appear when making then
        spherical approximation. These integrals are present
        at the step before selecting a scattering mechanism and
        as such can introduce the Fermi-Dirac integrals.
        Here one single valence band with the effective mass of
        an electron is tested against the closed Fermi-Dirac
        expressions for 100-700 K in 7 steps and for 10 different
        chemical potentials between -0.4 and 0.4 eV, where the band
        onset is set at 0 eV. Acoustic phonon scattering is used.

        """
        test_number = 2
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # set up transport
        tran = transport.Transport(bs)

        # now we do the closed version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = tran.sigma
        seebeck_closed = tran.seebeck
        lorenz_closed = tran.lorenz
        cc_closed = tran.ccp
        hall_closed = tran.hall
        # calculate small r factor
        small_hall_closed = 100 * constants.elcharge * \
            cc_closed * hall_closed

        # now we do the numeric version
        bs.param.transport_method = "numeric"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numeric = tran.sigma
        seebeck_numeric = tran.seebeck
        lorenz_numeric = tran.lorenz
        cc_numeric = tran.ccp
        hall_numeric = tran.hall

        # calculate small r factor
        small_hall_numeric = 100 * constants.elcharge * \
            cc_numeric * hall_numeric

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(invalid='ignore')
        difference_in_sigma = np.abs(
            (sigma_numeric - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.abs(
            (seebeck_numeric - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.abs(
            (lorenz_numeric - lorenz_closed) / lorenz_closed)
        difference_in_cc = np.abs(
            (cc_numeric - cc_closed) / cc_closed)
        difference_in_small_hall = np.abs(
            (small_hall_numeric - small_hall_closed) / small_hall_closed)

        # now test, only first element along diagonal (xx)
        # makes no sense to have something more tight than 1e-6
        # as most of the constants etc. are only given to 7 digits
        self.assertTrue(np.all(
            difference_in_sigma[:, :, 0, 0] < 1e-6))
        self.assertTrue(np.all(
            difference_in_seebeck[:, :, 0, 0] < 1e-6))
        self.assertTrue(np.all(
            difference_in_lorenz[:, :, 0, 0] < 1e-6))
        # non-degenerate limit is not so easy
        self.assertTrue(np.all(
            difference_in_cc[:, :, 0, 0] < 1e-2))
        self.assertTrue(np.all(
            difference_in_small_hall[:, :, 0, 0] < 1e-2))
        # should be better in the degenerate case
        self.assertTrue(np.all(
            difference_in_cc[:, 0:5, 0, 0] < 1e-6))
        self.assertTrue(np.all(
            difference_in_small_hall[:, 0:5, 0, 0] < 1e-6))

    def test_numeric_spherical_multi_band(self):
        """
        Same as the previous test, but here an additional
        valence band with an effective mass of 0.2 of the free
        electron mass is introduced to test multiband behavior.

        """
        test_number = 3
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the closed version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = tran.sigma
        seebeck_closed = tran.seebeck
        lorenz_closed = tran.lorenz
        cc_closed = tran.ccn
        hall_closed = tran.hall
        # calculate small r factor
        small_hall_closed = 2 * 100 * constants.elcharge * \
            cc_closed * hall_closed

        # now we do the numeric version
        bs.param.transport_method = "numeric"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numeric = tran.sigma
        seebeck_numeric = tran.seebeck
        lorenz_numeric = tran.lorenz
        cc_numeric = tran.ccn
        hall_numeric = tran.hall

        # calculate small r factor
        small_hall_numeric = 2 * 100 * constants.elcharge * \
            cc_numeric * hall_numeric

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(invalid='ignore')
        difference_in_sigma = np.abs(
            (sigma_numeric - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.abs(
            (seebeck_numeric - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.abs(
            (lorenz_numeric - lorenz_closed) / lorenz_closed)
        difference_in_cc = np.abs(
            (cc_numeric - cc_closed) / cc_closed)
        difference_in_small_hall = np.abs(
            (small_hall_numeric - small_hall_closed) / small_hall_closed)

        # now test, only first element along diagonal (xx)
        # makes no sense to have something more tight than 1e-6
        # as most of the constants etc. are only given to 7 digits
        self.assertTrue(np.all(
            difference_in_sigma[:, :, 0, 0] < 1e-6))
        self.assertTrue(np.all(
            difference_in_seebeck[:, :, 0, 0] < 1e-6))
        self.assertTrue(np.all(
            difference_in_lorenz[:, :, 0, 0] < 1e-6))
        # non-degenerate limit is not so easy
        self.assertTrue(np.all(
            difference_in_cc[:, :, 0, 0] < 1e-2))
        self.assertTrue(np.all(
            difference_in_small_hall[:, :, 0, 0] < 1e-2))
        # should be better in the degenerate case
        self.assertTrue(np.all(
            difference_in_cc[:, 5:9, 0, 0] < 1e-6))
        self.assertTrue(np.all(
            difference_in_small_hall[:, 5:9, 0, 0] < 1e-6))

    def test_trapz_spherical_one_band_fast(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 700 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 31x31x31.

        """
        test_number = 4
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the numeric version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = tran.sigma
        seebeck_closed = tran.seebeck
        lorenz_closed = tran.lorenz

        # now we do the numerick version
        bs.param.transport_method = "numerick"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numerick = tran.sigma
        seebeck_numerick = tran.seebeck
        lorenz_numerick = tran.lorenz

        # now calculate the relative difference
        # ignore invalid value and divide by zero errors
        np.seterr(invalid='ignore', divide='ignore')
        difference_in_sigma = np.abs(
            (sigma_numerick - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.abs(
            (seebeck_numerick - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.abs(
            (lorenz_numerick - lorenz_closed) / lorenz_closed)

        # now test, only first element along diagonal (xx)
        # at the supplied grid of 31x31x31 they all should
        # match to within 1e-3
        self.assertTrue(np.all(
            difference_in_sigma[:, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_seebeck[:, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_lorenz[:, :, 0, 0] < 1e-3))

    def test_trapz_preinter_spherical_one_band_fast(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 300, 600 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 11x11x11, which is
        interpolated (Akima) to 31x31x31 and then integrated by
        trapezoidal.

        """
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)

        test_number = 7
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the numeric version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = tran.sigma
        seebeck_closed = tran.seebeck
        lorenz_closed = tran.lorenz

        # now we do the numerick version
        bs.param.transport_method = "numerick"
        # interpolate
        logger.info("Pre-interpolating the bandstructure.")
        bs.interpolate(store_inter=True, ivelocities=True)
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numerick = tran.sigma
        seebeck_numerick = tran.seebeck
        lorenz_numerick = tran.lorenz

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(divide='ignore', invalid='ignore')
        difference_in_sigma = np.abs(
            (sigma_numerick - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.abs(
            (seebeck_numerick - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.abs(
            (lorenz_numerick - lorenz_closed) / lorenz_closed)

        # this data is surely not very well converged at low
        # temperatures but it still suffices as a test case for
        # preinterpolation

        # 300 K
        self.assertTrue(np.all(
            difference_in_sigma[0, :, 0, 0] < 1e-2))
        self.assertTrue(np.all(
            difference_in_seebeck[0, :, 0, 0] < 1e-1))
        self.assertTrue(np.all(
            difference_in_lorenz[0, :, 0, 0] < 1e-1))
        # 700 K
        self.assertTrue(np.all(
            difference_in_sigma[1, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_seebeck[1, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_lorenz[1, :, 0, 0] < 1e-3))

    def may_data(self):
        """
        Values of a few transport coefficients calculated with
        the closed analytic Fermi-Dirac expressions for the
        free electron mass at 300 K. The values assume acoustic
        phonon scattering.

        Parameters
        ----------
        None

        Returns
        -------
        data : ndarray
            | Dimension: (11,5)

            Contains the data from Table 11.1 in
            Ref.:cite:`rowe_2012_mpacit_may`. With 11 entries
            of different unitless chemical potentials (first index
            of the second axis). The remaining 4 indexes of the
            second axis contains the Seebeck coefficient, Lorenz
            coefficient, the charge carrier density and the unitless
            Hall coefficient. The units for the three first are

            ..math:: \\micro \\mathrm{V/K}, 10^{-8} \\mathrm{V}^2/
                 \\mathrm{K}^2, 10^{21} \\mathrm{cm}^{-3},
            respectively.

        Notes
        ----
        Used in this implementation to tests the closed analytic
        Fermi-Dirac integrals. The data given by may is converted to
        a conduction band and thus the Seebeck should be
        negative. Otherwise the values entered should be precicely the
        same as in Ref.:cite:`rowe_2012_mpacit_may`.

        """
        data = np.zeros((11, 5))

        # unitless eta values (chempot/k_bT)
        # data[:, 0] = [-3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10]
        # data[:, 0] = [-10, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3]
        data[:, 0] = [3, 2, 1, 0, -1, -2, -3, -4, -6, -8, -10]
        # Seebeck coefficient in microV/K
        data[:, 1] = [432.97, 350.24, 272.50, 204.50,
                      150.88, 112.39, 86.09, 68.22, 46.95, 35.40, 28.35]
        # Lorenz coefficient in 10^-8 V^2/K^2
        data[:, 2] = [1.49, 1.51, 1.54, 1.61,
                      1.72, 1.86, 1.99, 2.09, 2.24, 2.32, 2.36]
        # carrier concentration in 10^21 cm^-3
        data[:, 3] = [0.00123, 0.00324, 0.0082, 0.0192, 0.0395,
                      0.0709, 0.1126, 0.1634, 0.2872, 0.4355, 0.6044]
        # unitless hall coefficient
        data[:, 4] = [1.17, 1.17, 1.16, 1.13,
                      1.11, 1.08, 1.05, 1.04, 1.02, 1.01, 1.01]
        return data


class TestTransportIntegralsSlow(unittest.TestCase):

    def test_trapz_spherical_one_band(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 100, 300, 600 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 61x61x61.

        """
        test_number = 5
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the numeric version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = tran.sigma
        seebeck_closed = tran.seebeck
        lorenz_closed = tran.lorenz

        # now we do the numerick version
        bs.param.transport_method = "numerick"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numerick = tran.sigma
        seebeck_numerick = tran.seebeck
        lorenz_numerick = tran.lorenz

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(divide='ignore', invalid='ignore')
        difference_in_sigma = np.abs(
            (sigma_numerick - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.abs(
            (seebeck_numerick - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.abs(
            (lorenz_numerick - lorenz_closed) / lorenz_closed)

        # now test, only first element along diagonal (xx)
        # at the supplied grid of 61x61x61 they all should
        # match to within 1e-2 (sigma) for 100 and 300 K and
        # 1e-3 for 700 K.

        # 100 K, full range of chempot, sigma and lorenz
        self.assertTrue(np.all(
            difference_in_sigma[0, :, 0, 0] < 1e-2))
        self.assertTrue(np.all(
            difference_in_lorenz[0, :, 0, 0] < 1e-1))

        # 300 K, full range of chempot, all
        self.assertTrue(np.all(
            difference_in_sigma[1, :, 0, 0] < 1e-2))
        self.assertTrue(np.all(
            difference_in_seebeck[1, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_lorenz[1, :, 0, 0] < 1e-3))

        # 700 K, full range of  chempot, all
        self.assertTrue(np.all(
            difference_in_sigma[2, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_seebeck[2, :, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_lorenz[2, :, 0, 0] < 1e-3))

        # some additional corner cases
        # 100 K last chempot entry (non-deg limit) for seebeck
        # and lorenz
        self.assertTrue(np.all(
            difference_in_seebeck[0, 9, 0, 0] < 1e-3))
        self.assertTrue(np.all(
            difference_in_lorenz[0, 9, 0, 0] < 1e-2))

        # 300 K first chempot (degen limit) sigma
        self.assertTrue(np.all(
            difference_in_sigma[1, 0, 0, 0] < 1e-5))

        # 700 K first chempot (deg limit) for all
        self.assertTrue(np.all(
            difference_in_sigma[2, 0, 0, 0] < 1e-5))
        self.assertTrue(np.all(
            difference_in_seebeck[2, 0, 0, 0] < 1e-5))
        self.assertTrue(np.all(
            difference_in_lorenz[2, 0, 0, 0] < 1e-5))

    def test_tetra_spherical_one_band(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 100, 300, 600 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 61x61x61.

        """
        test_number = 6
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the numeric version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = tran.sigma
        seebeck_closed = tran.seebeck
        lorenz_closed = tran.lorenz

        # now we do the numerick version
        bs.param.transport_method = "numerick"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numerick = tran.sigma
        seebeck_numerick = tran.seebeck
        lorenz_numerick = tran.lorenz

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(divide='ignore', invalid='ignore')
        difference_in_sigma = np.abs(
            (sigma_numerick - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.abs(
            (seebeck_numerick - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.abs(
            (lorenz_numerick - lorenz_closed) / lorenz_closed)

        # for the tetrahedron method we need very dense grids
        # so the accuracy is not so high for the input 61x61x61 grid
        # for the current implementation

        # five first chempots for sigma
        self.assertTrue(np.all(
            difference_in_sigma[:, 0:4, 0, 0] < 1e-1))
        # five first chempots for seebeck
        self.assertTrue(np.all(
            difference_in_seebeck[:, 0:4, 0, 0] < 1e-2))
        # five first chempots for lorenz
        self.assertTrue(np.all(
            difference_in_lorenz[:, 0:4, 0, 0] < 1e-2))

    def test_trapz_preinter_spherical_one_band(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 100, 300, 600 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 11x11x11, which is
        interpolated (Akima) to 61x61x61 and then integrated by
        trapezoidal.

        """
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)

        test_number = 8
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the numeric version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_closed = np.nan_to_num(tran.sigma)
        seebeck_closed = np.nan_to_num(tran.seebeck)
        lorenz_closed = np.nan_to_num(tran.lorenz)

        # now we do the numerick version
        bs.param.transport_method = "numerick"
        # interpolate
        logger.info("Pre-interpolating the bandstructure.")
        bs.interpolate(store_inter=True, ivelocities=True)
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numerick = np.nan_to_num(tran.sigma)
        seebeck_numerick = np.nan_to_num(tran.seebeck)
        lorenz_numerick = np.nan_to_num(tran.lorenz)

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(divide='ignore', invalid='ignore')
        difference_in_sigma = np.nan_to_num(
            (sigma_numerick - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.nan_to_num(
            (seebeck_numerick - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.nan_to_num(
            (lorenz_numerick - lorenz_closed) / lorenz_closed)

        # this data should be reasonably converged, at least at
        # 300 and 600 K and close to a few percent at
        # 100 K.

        # 100 K
        self.assertTrue(np.all(np.abs(
            difference_in_sigma[0, :, 0, 0]) < 1e-2))
        # (only the last 5 chempots)
        self.assertTrue(np.all(np.abs(
            difference_in_seebeck[0, 4:9, 0, 0]) < 1e-1))
        self.assertTrue(np.all(np.abs(
            difference_in_lorenz[0, :, 0, 0]) < 1e-2))
        # 300 K
        self.assertTrue(np.all(np.abs(
            difference_in_sigma[1, :, 0, 0]) < 1e-2))
        self.assertTrue(np.all(np.abs(
            difference_in_seebeck[1, :, 0, 0]) < 1e-3))
        self.assertTrue(np.all(np.abs(
            difference_in_lorenz[1, :, 0, 0]) < 1e-3))
        # 700 K
        self.assertTrue(np.all(np.abs(
            difference_in_sigma[2, :, 0, 0]) < 1e-3))
        self.assertTrue(np.all(np.abs(
            difference_in_seebeck[2, :, 0, 0]) < 1e-3))
        self.assertTrue(np.all(np.abs(
            difference_in_lorenz[2, :, 0, 0]) < 1e-3))

    def test_cubature_spherical_one_band(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 100, 300, 600 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 11x11x11 on which
        Cubature integration with on-the-fly Akima intepolation
        is used. Isotropy is enforced and the relative error is
        set at 10^-2.

        """

        test_number = 9
        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # now we do the numeric version
        bs.param.transport_method = "closed"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        bs.param.libinfo = False
        sigma_closed = np.nan_to_num(tran.sigma)
        seebeck_closed = np.nan_to_num(tran.seebeck)
        lorenz_closed = np.nan_to_num(tran.lorenz)

        # now we do the numerick version
        bs.param.transport_method = "numerick"
        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_numerick = np.nan_to_num(tran.sigma)
        seebeck_numerick = np.nan_to_num(tran.seebeck)
        lorenz_numerick = np.nan_to_num(tran.lorenz)

        # now calculate the relative difference
        # ignore invalid value errors
        np.seterr(divide='ignore', invalid='ignore')
        difference_in_sigma = np.nan_to_num(
            (sigma_numerick - sigma_closed) / sigma_closed)
        difference_in_seebeck = np.nan_to_num(
            (seebeck_numerick - seebeck_closed) / seebeck_closed)
        difference_in_lorenz = np.nan_to_num(
            (lorenz_numerick - lorenz_closed) / lorenz_closed)

        # this data should be reasonably converged
        # 300 K
        self.assertTrue(np.all(np.abs(
            difference_in_sigma[0, :, 0, 0]) < 1e-2))
        self.assertTrue(np.all(np.abs(
            difference_in_seebeck[0, :, 0, 0]) < 1e-2))
        # only first 5 chempots
        self.assertTrue(np.all(np.abs(
            difference_in_lorenz[0, 0:4, 0, 0]) < 1e-2))
        # 700 K
        self.assertTrue(np.all(np.abs(
            difference_in_sigma[1, :, 0, 0]) < 1e-3))
        self.assertTrue(np.all(np.abs(
            difference_in_seebeck[1, :, 0, 0]) < 1e-2))
        # only first 5 chempots
        self.assertTrue(np.all(np.abs(
            difference_in_lorenz[1, 0:4, 0, 0]) < 1e-2))


def read_and_setup_bandstructure(test_number):
    """
    This routine reads the parameter files, sets up
    the lattice and calculates or reads in the bandstructure

    Parameters
    ----------
    test_number : int
        The test number, e.g. the folder under the tests
        folder.

    Returns
    -------
    bs : object
        A `Bandstructure()` object.

    """

    # unfortunately we cannot pickle with classes
    # so load data everytime this is called
    location = os.path.dirname(__file__) + \
        "/tests/" + str(test_number)
    # read parameters, generate lattice and bands
    param = inputoutput.Param(inputoutput.readparam(
        location=location))
    lat = lattice.Lattice(param, location=location)
    bs = bandstructure.Bandstructure(lat, param,
                                     location=location)
    return bs


def pickle_test_data(location):
    """
    This routine simply pickles the test data if it is
    not already present in order to speed up testing

    Parameters
    ----------
    test_number : int
        An integer describing which test to be pickled

    Returns
    -------
    None

    Notes
    -----
    Stores the pickled object `data` in the tests/test_number
    folder for later use

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)

    # read parameters, generate lattice and bands
    param = inputoutput.Param(inputoutput.readparam(
        location=location))
    lat = lattice.Lattice(param, location=location)
    bs = bandstructure.Bandstructure(lat, param,
                                     location=location)

    # now pickle data to file so that we can load it
    # later in the test set
    data_file = inputoutput.file_handler(
        filename=location + "/data", status="w")
    # pickle pickle (lat and param should now be in bs)
    pickle.dump(bs, data_file)

    inputoutput.file_handler(file_handler=data_file)


def test_skw(param):
    """
    This is a simple test of the SKW method, which is similar to the CoSb3
    test supplied with Boltztrap, for comparison reasons.

    Parameters
    ----------
    None

    Returns
    -------
    None

    Notes
    -----
    Prints the supplied energies and k-points and the interpolated
    k-point grid and its energies at the end.

    The density of the interpolated grid can be controlled with the
    tag `skw_expansion_factor` in param.yml.

    """

    # set logger
    logger = logging.getLogger(__name__)
    logger.debug("Running test_skw.")

    unitcell = np.zeros((3, 3), dtype="double")
    unitcell[0][0] = 9.0164589999999993
    unitcell[1][1] = 9.0164589999999993
    unitcell[2][2] = 9.0164589999999993

    positions = np.array([
        [0.2499994731856461, 0.2499994731856461, 0.2499994731856461],
        [0.7499995286397905, 0.7499995286397905, 0.7499995286397905],
        [0.7499995286397905, 0.7499995286397905, 0.2499994731856461],
        [0.2499994731856461, 0.2499994731856461, 0.7499995286397905],
        [0.7499995286397905, 0.2499994731856461, 0.7499995286397905],
        [0.2499994731856461, 0.7499995286397905, 0.2499994731856461],
        [0.2499994731856461, 0.7499995286397905, 0.7499995286397905],
        [0.7499995286397905, 0.2499994731856461, 0.2499994731856461],
        [0.0000000000000000, 0.3351005089692052, 0.1580498508339048],
        [0.0000000000000000, 0.6648996019390765, 0.8419502600743840],
        [0.0000000000000000, 0.6648996019390765, 0.1580498508339048],
        [0.0000000000000000, 0.3351005089692052, 0.8419502600743840],
        [0.1580498508339048, 0.0000000000000000, 0.3351005089692052],
        [0.8419502600743840, 0.0000000000000000, 0.6648996019390765],
        [0.1580498508339048, 0.0000000000000000, 0.6648996019390765],
        [0.8419502600743840, 0.0000000000000000, 0.3351005089692052],
        [0.3351005089692052, 0.1580498508339048, 0.0000000000000000],
        [0.6648996019390765, 0.8419502600743840, 0.0000000000000000],
        [0.6648996019390765, 0.1580498508339048, 0.0000000000000000],
        [0.3351005089692052, 0.8419502600743840, 0.0000000000000000],
        [0.5000000554541444, 0.8351005644233496, 0.6580499062880421],
        [0.5000000554541444, 0.1648995464849321, 0.3419502046202396],
        [0.5000000554541444, 0.1648995464849321, 0.6580499062880421],
        [0.5000000554541444, 0.8351005644233496, 0.3419502046202396],
        [0.6580499062880421, 0.5000000554541444, 0.8351005644233496],
        [0.3419502046202396, 0.5000000554541444, 0.1648995464849321],
        [0.6580499062880421, 0.5000000554541444, 0.1648995464849321],
        [0.3419502046202396, 0.5000000554541444, 0.8351005644233496],
        [0.8351005644233496, 0.6580499062880421, 0.5000000554541444],
        [0.1648995464849321, 0.3419502046202396, 0.5000000554541444],
        [0.1648995464849321, 0.6580499062880421, 0.5000000554541444],
        [0.8351005644233496, 0.3419502046202396, 0.5000000554541444]], dtype="double")

    species = np.array([27, 27, 27, 27, 27, 27, 27, 27,
                        57, 57, 57, 57, 57, 57, 57, 57,
                        57, 57, 57, 57, 57, 57, 57, 57,
                        57, 57, 57, 57, 57, 57, 57, 57], dtype="intc")

    kpoints = np.array([
        [0, 0, 0],
        [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
        [0.5, 0.5, 0.5]], dtype="double")

    energies = np.array([
        [-0.94972174,  -0.90915329,  -0.86640789],
        [-0.85421954,  -0.87713215,  -0.86640789],
        [-0.85421954,  -0.87713127,  -0.86640782],
        [-0.8117333,  -0.8418362,  -0.86640628],
        [-0.8117333,  -0.83695244,  -0.86640628],
        [-0.8117333,  -0.83695163,  -0.86640628],
        [-0.78690561,  -0.7844279,  -0.80302305],
        [-0.76824218,  -0.7742723,  -0.80302041],
        [-0.76824218,  -0.72629446,  -0.69086783],
        [-0.72479965,  -0.7262938,  -0.69086783],
        [-0.72479965,  -0.70679102,  -0.69086776],
        [-0.72479965,  -0.70679036,  -0.69086761],
        [-0.70295505,  -0.70322729,  -0.69086754],
        [-0.70295505,  -0.6978075,  -0.69086754],
        [-0.70295505,  -0.68912496,  -0.68539454],
        [-0.65385687,  -0.68912489,  -0.68539454],
        [-0.65385687,  -0.6696626,  -0.68539204],
        [-0.65385687,  -0.66966209,  -0.68539204],
        [-0.63934376,  -0.58215559,  -0.56922571],
        [-0.63934376,  -0.57712769,  -0.56922571],
        [-0.63934376,  -0.57712659,  -0.56922564],
        [-0.58730931,  -0.57437296,  -0.56922417],
        [-0.58730931,  -0.56799071,  -0.56922402],
        [-0.58730931,  -0.5679902,  -0.56922402],
        [-0.39107138,  -0.41290113,  -0.43204701],
        [-0.38480849,  -0.40850716,  -0.43204701],
        [-0.37177571,  -0.40850466,  -0.43204605],
        [-0.37177571,  -0.40739248,  -0.43204554],
        [-0.34801927,  -0.40537465,  -0.43204458],
        [-0.30569553,  -0.40537281,  -0.43204458],
        [-0.30569553,  -0.33896306,  -0.33874587],
        [-0.30569553,  -0.33896152,  -0.33874543],
        [-0.30156079,  -0.3378525,  -0.30972458],
        [-0.30156079,  -0.33160673,  -0.30972458],
        [-0.30156079,  -0.31966901,  -0.30972458],
        [-0.2934339,  -0.31966901,  -0.30972311],
        [-0.2934339,  -0.31537199,  -0.30972311],
        [-0.2934339,  -0.30710413,  -0.30972311],
        [-0.28600657,  -0.28993912,  -0.30401756],
        [-0.28600657,  -0.28993787,  -0.30401756],
        [-0.28600657,  -0.27609374,  -0.30401675],
        [-0.26674376,  -0.27609352,  -0.30401675],
        [-0.2589754,  -0.25213981,  -0.21562775],
        [-0.2589754,  -0.23631474,  -0.21562775],
        [-0.22044639,  -0.22559534,  -0.21562775],
        [-0.22044639,  -0.22559519,  -0.21562731],
        [-0.22044639,  -0.21610512,  -0.21562723],
        [-0.21893967,  -0.2161049,  -0.21562723],
        [-0.21893967,  -0.20005023,  -0.19333995],
        [-0.21893967,  -0.19947702,  -0.19333995],
        [-0.19831625,  -0.19947525,  -0.19333995],
        [-0.19831625,  -0.19495957,  -0.19333966],
        [-0.17533022,  -0.19214134,  -0.19333966],
        [-0.17417438,  -0.19053841,  -0.18735761],
        [-0.17417438,  -0.18954008,  -0.18735761],
        [-0.17222248,  -0.18953986,  -0.18735739],
        [-0.17222248,  -0.18293828,  -0.18735695],
        [-0.15854755,  -0.18102107,  -0.18735673],
        [-0.15854755,  -0.18102026,  -0.18735673],
        [-0.15305367,  -0.17789187,  -0.16730219],
        [-0.15305367,  -0.16273697,  -0.16730219],
        [-0.15305367,  -0.16273631,  -0.16730211],
        [-0.14598648,  -0.1581078,  -0.16730211],
        [-0.14598648,  -0.15810729,  -0.14824921],
        [-0.14598648,  -0.14720943,  -0.14824914],
        [-0.14481962,  -0.14720921,  -0.14824914],
        [-0.14481962,  -0.13217823,  -0.14824855],
        [-0.14481962,  -0.12921278,  -0.14824855],
        [-0.14250125,  -0.12581817,  -0.14824848],
        [-0.14250125,  -0.12581795,  -0.14287757],
        [-0.14250125,  -0.12333472,  -0.14287676],
        [-0.13534572,  -0.12286389,  -0.11372228],
        [-0.13534572,  -0.12286382,  -0.11372213],
        [-0.13534572,  -0.12081049,  -0.11372213],
        [-0.12117761,  -0.11781954,  -0.11372177],
        [-0.12117761,  -0.11781932,  -0.11372177],
        [-0.12117761,  -0.10871578,  -0.11372155],
        [-0.11731828,  -0.09701899,  -0.08945399],
        [-0.11731828,  -0.09701869,  -0.08945399],
        [-0.11646665,  -0.09593973,  -0.08945355],
        [-0.11646665,  -0.09388787,  -0.08945355],
        [-0.11646665,  -0.09388765,  -0.08753215],
        [-0.11293967,  -0.08820724,  -0.08753215],
        [-0.11293967,  -0.0876937,  -0.08753215],
        [-0.11293967,  -0.08769333,  -0.08753186],
        [-0.11083834,  -0.08256636,  -0.08753186],
        [-0.10475273,  -0.08256629,  -0.08753186],
        [-0.10475273,  -0.07325681,  -0.06237931],
        [-0.09597986,  -0.06509802,  -0.06237909],
        [-0.09597986,  -0.06208156,  -0.06237909],
        [-0.09597986,  -0.0620812,  -0.06237857],
        [-0.07047099,  -0.05647215,  -0.06237857],
        [-0.07047099,  -0.056472,  -0.06237835],
        [-0.06276569,  -0.05329774,  -0.04241186],
        [-0.00309517,  -0.04826316,  -0.04241047],
        [0.00457536,  0.03268609,  0.0443007],
        [0.00457536,  0.03268609,  0.04430151],
        [0.00457536,  0.03990666,  0.09492891],
        [0.01799753,  0.04047157,  0.09492898],
        [0.01799753,  0.05047128,  0.09492898],
        [0.01799753,  0.05047165,  0.09492942],
        [0.02620358,  0.05817996,  0.09492942],
        [0.02620358,  0.07763586,  0.09492949],
        [0.02620358,  0.10080644,  0.10527486],
        [0.03509919,  0.10080644,  0.10527486],
        [0.03509919,  0.10083415,  0.10527486],
        [0.11249162,  0.10083496,  0.10527538],
        [0.11249162,  0.11985692,  0.10527538],
        [0.11249162,  0.15708691,  0.10527538],
        [0.12245745,  0.1570875,  0.13995195],
        [0.12245745,  0.16400724,  0.13995195],
        [0.12245745,  0.16634354,  0.13995364],
        [0.13765475,  0.16634376,  0.13995364],
        [0.13765475,  0.17867625,  0.1718006],
        [0.13765475,  0.17867639,  0.1718006],
        [0.23309065,  0.18103202,  0.1718014],
        [0.25537551,  0.18489217,  0.17180148],
        [0.25537654,  0.18489254,  0.17180221],
        [0.2553855,  0.19409912,  0.17180221]], dtype="double")

    ikpoints, ienergies = skw_interface.interpolate_energies(
        energies, kpoints, unitcell, positions, species, param.skw_expansion_factor)
    print(kpoints)
    print(energies)
    print(ikpoints)
    print(ienergies)