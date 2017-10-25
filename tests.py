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

    def test_closed_parabolic(self):
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

    def test_numeric_parabolic_one_band(self):
        """
        Test the numerical solution of the energy
        dependent integral that appear when making then
        parabolic approximation. These integrals are present
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

    def test_numeric_parabolic_multi_band(self):
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

    def test_trapz_parabolic_one_band_fast(self):
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

    def test_trapz_parabolic_one_band_tri_fast(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 700 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. The supplied grid is 31x31x31.
        Here, we perform the calculation in a trigonal cell,
        which should be the same as the cubic (or the closed)
        if all metric is in ordnung.

        """
        test_number = 11
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

    def test_trapz_parabolic_one_band_numdiff_fast(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 700 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. Here, the velocities are recalculated
        by numerical difference. The supplied grid is 31x31x31.

        """
        test_number = 10
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

        # and make sure to turn on the recalculation of
        # the velocities by numerical differences
        bs.param.dispersion_velocities_numdiff = True
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

    def test_trapz_parabolic_one_band_tri_numdiff_fast(self):
        """
        Test the direct integrals in reciprocal space against
        the results of the closed Fermi-Dirac integrals
        for a single valence band with a free electron mass,
        acoustic phonon scattering, at 700 K for
        chemical potential between -0.4 and 0.4 eV
        in 10 steps. Here, the velocities are recalculated
        by numerical difference. The supplied grid is 31x31x31.
        Here, we perform the calculation in a trigonal cell,
        which should be the same as the cubic (or the closed)
        if all metric is in ordnung.

        """
        test_number = 12
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

        # and make sure to turn on the recalculation of
        # the velocities by numerical differences
        bs.param.dispersion_velocities_numdiff = True
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

    def test_trapz_preinter_parabolic_one_band_fast(self):
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

    def test_si_45_primitive(self):
        """
        Test where the primitive cell of silicon is
        used at an input grid density of 45x45x45.
        This is not fully converged, but should still
        yield a reprentative test which checks the
        integration of a first-principle input in
        a cell which is not cubic. The chemical potential
        runs between -0.4 and 1.0 eV (band gap of rougly 0.6 eV).
        Test performed at 300 K.

        """

        test_number = 13

        # fetch bandstructure
        bs = read_and_setup_bandstructure(test_number)

        # set up transport
        tran = transport.Transport(bs)
        tran.calc_transport_tensors()
        sigma_calc = np.nan_to_num(tran.sigma)
        seebeck_calc = np.nan_to_num(tran.seebeck)
        lorenz_calc = np.nan_to_num(tran.lorenz)

        # set up reference values to check against
        # only check values at -0.4, 0.0 and 1.0 eV
        sigma_ref = np.array([1.14673190e+07, 4.13673024e+05,
                              1.21445728e+07])
        seebeck_ref = np.array([2.64804235e+01, 1.76976501e+02,
                                -3.08544065e+01])
        lorenz_ref = np.array([2.16706616e+00, 1.74093914e+00,
                               2.38318180e+00])

        # now calculate the relative difference
        np.seterr(divide='ignore', invalid='ignore')
        difference_in_sigma = np.nan_to_num(
            (sigma_calc[0, [0, 5, 19], 0, 0] -
             sigma_ref) / sigma_ref)
        difference_in_seebeck = np.nan_to_num(
            (seebeck_calc[0, [0, 5, 19], 0, 0] -
             seebeck_ref) / seebeck_ref)
        difference_in_lorenz = np.nan_to_num(
            (lorenz_calc[0, [0, 5, 19], 0, 0] -
             lorenz_ref) / lorenz_ref)

        # should match down to numerical precision for all chemical
        # potentials
        self.assertTrue(np.all(np.abs(
            difference_in_sigma[:]) < 1e-7))
        self.assertTrue(np.all(np.abs(
            difference_in_seebeck[:]) < 1e-7))
        self.assertTrue(np.all(np.abs(
            difference_in_lorenz[:]) < 1e-7))

    def test_trapz_parabolic_one_band(self):
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

    def test_tetra_parabolic_one_band(self):
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

    def test_trapz_preinter_parabolic_one_band(self):
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

    def test_cubature_parabolic_one_band(self):
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
