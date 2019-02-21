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

# pylint: disable=unused-import, useless-import-alias

import pytest
import numpy as np
import logging

import t4me.constants as constants
import t4me.transport as transport
from .fixtures import may_data
from .fixtures import read_and_setup_bs


@pytest.mark.parametrize('read_and_setup_bs', ('closed_parabolic',), indirect=True)
def test_closed_parabolic(may_data, read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")

    """
    Test the Seebeck and Lorenz coefficient, carrier
    concentration and the small-h Hall factor against
    the data published by May :cite:`rowe_2012_mpacit_may`
    for a selected set of unitless chemical potentials at
    300 K for a band with the effective mass of a free
    electron. Acoustic phonon scattering is used.

    """

    bs = read_and_setup_bs

    # set up transport
    tran = transport.Transport(bs)

    # fetch chempots from etas at 300 K (may_data)
    temperature = 300
    tran.temperatures = np.array([temperature])
    tran.chempots = transport.fetch_chempot_from_etas(temperature, may_data[:, 0])

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
    this_data = np.column_stack((seebeck.T, lorenz.T, cc.T, small_hall.T))

    # now calculate the relative difference
    difference = np.abs((this_data - may_data[:, 1:5]) / may_data[:, 1:5])

    # print difference
    # now test (only two decimals)
    assert np.all(difference < 1e-2)


@pytest.mark.parametrize('read_and_setup_bs', ('numeric_parabolic_one_band',), indirect=True)
def test_numeric_parabolic_one_band(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
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

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numeric - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numeric - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numeric - lorenz_closed) / lorenz_closed)
    difference_in_cc = np.abs((cc_numeric - cc_closed) / cc_closed)
    difference_in_small_hall = np.abs((small_hall_numeric - small_hall_closed) / small_hall_closed)

    # now test, only first element along diagonal (xx)
    # makes no sense to have something more tight than 1e-6
    # as most of the constants etc. are only given to 7 digits
    assert np.all(difference_in_sigma[:, :, 0, 0] < 1e-6)
    assert np.all(difference_in_seebeck[:, :, 0, 0] < 1e-6)
    assert np.all(difference_in_lorenz[:, :, 0, 0] < 1e-6)
    # non-degenerate limit is not so easy
    assert np.all(difference_in_cc[:, :, 0, 0] < 1e-2)
    assert np.all(difference_in_small_hall[:, :, 0, 0] < 1e-2)
    # should be better in the degenerate case
    assert np.all(difference_in_cc[:, 0:5, 0, 0] < 1e-6)
    assert np.all(difference_in_small_hall[:, 0:5, 0, 0] < 1e-6)


@pytest.mark.parametrize('read_and_setup_bs', ['numeric_parabolic_multi_band'], indirect=True)
def test_numeric_parabolic_multi_band(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
    """
    Same as the previous test, but here an additional
    valence band with an effective mass of 0.2 of the free
    electron mass is introduced to test multiband behavior.

    """

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numeric - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numeric - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numeric - lorenz_closed) / lorenz_closed)
    difference_in_cc = np.abs((cc_numeric - cc_closed) / cc_closed)
    difference_in_small_hall = np.abs((small_hall_numeric - small_hall_closed) / small_hall_closed)

    # now test, only first element along diagonal (xx)
    # makes no sense to have something more tight than 1e-6
    # as most of the constants etc. are only given to 7 digits
    assert np.all(difference_in_sigma[:, :, 0, 0] < 1e-6)
    assert np.all(difference_in_seebeck[:, :, 0, 0] < 1e-6)
    assert np.all(difference_in_lorenz[:, :, 0, 0] < 1e-6)
    # non-degenerate limit is not so easy
    assert np.all(difference_in_cc[:, :, 0, 0] < 1e-2)
    assert np.all(difference_in_small_hall[:, :, 0, 0] < 1e-2)
    # should be better in the degenerate case
    assert np.all(difference_in_cc[:, 5:9, 0, 0] < 1e-6)
    assert np.all(difference_in_small_hall[:, 5:9, 0, 0] < 1e-6)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_parabolic_one_band_fast'], indirect=True)
def test_trapz_parabolic_one_band_fast(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
    """
    Test the direct integrals in reciprocal space against
    the results of the closed Fermi-Dirac integrals
    for a single valence band with a free electron mass,
    acoustic phonon scattering, at 700 K for
    chemical potential between -0.4 and 0.4 eV
    in 10 steps. The supplied grid is 31x31x31.

    """

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # now test, only first element along diagonal (xx)
    # at the supplied grid of 31x31x31 they all should
    # match to within 1e-3
    assert np.all(difference_in_sigma[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_seebeck[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[:, :, 0, 0] < 1e-3)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_parabolic_one_band_tri_fast'], indirect=True)
def test_trapz_parabolic_one_band_tri_fast(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
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

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # now test, only first element along diagonal (xx)
    # at the supplied grid of 31x31x31 they all should
    # match to within 1e-3
    assert np.all(difference_in_sigma[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_seebeck[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[:, :, 0, 0] < 1e-3)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_parabolic_one_band_numdiff_fast'], indirect=True)
def test_trapz_parabolic_one_band_numdiff_fast(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
    """
    Test the direct integrals in reciprocal space against
    the results of the closed Fermi-Dirac integrals
    for a single valence band with a free electron mass,
    acoustic phonon scattering, at 700 K for
    chemical potential between -0.4 and 0.4 eV
    in 10 steps. Here, the velocities are recalculated
    by numerical difference. The supplied grid is 31x31x31.

    """

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # now test, only first element along diagonal (xx)
    # at the supplied grid of 31x31x31 they all should
    # match to within 1e-3
    assert np.all(difference_in_sigma[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_seebeck[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[:, :, 0, 0] < 1e-3)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_parabolic_one_band_tri_numdiff_fast'], indirect=True)
def test_trapz_parabolic_one_band_tri_numdiff_fast(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
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

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # now test, only first element along diagonal (xx)
    # at the supplied grid of 31x31x31 they all should
    # match to within 1e-3
    assert np.all(difference_in_sigma[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_seebeck[:, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[:, :, 0, 0] < 1e-3)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_preinter_parabolic_one_band_fast'], indirect=True)
def test_trapz_preinter_parabolic_one_band_fast(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
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

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # this data is surely not very well converged at low
    # temperatures but it still suffices as a test case for
    # preinterpolation

    # 300 K
    assert np.all(difference_in_sigma[0, :, 0, 0] < 1e-2)
    assert np.all(difference_in_seebeck[0, :, 0, 0] < 1e-1)
    assert np.all(difference_in_lorenz[0, :, 0, 0] < 1e-1)
    # 700 K
    assert np.all(difference_in_sigma[1, :, 0, 0] < 1e-3)
    assert np.all(difference_in_seebeck[1, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[1, :, 0, 0] < 1e-3)


@pytest.mark.parametrize('read_and_setup_bs', ['si_45_primitive'], indirect=True)
def test_si_45_primitive(read_and_setup_bs):
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

    bs = read_and_setup_bs

    # set up transport
    tran = transport.Transport(bs)
    tran.calc_transport_tensors()
    sigma_calc = np.nan_to_num(tran.sigma)
    seebeck_calc = np.nan_to_num(tran.seebeck)
    lorenz_calc = np.nan_to_num(tran.lorenz)

    # set up reference values to check against
    # only check values at -0.4, 0.0 and 1.0 eV
    sigma_ref = np.array([1.14673190e+07, 4.13673024e+05, 1.21445728e+07])
    seebeck_ref = np.array([2.64804235e+01, 1.76976501e+02, -3.08544065e+01])
    lorenz_ref = np.array([2.16706616e+00, 1.74093914e+00, 2.38318180e+00])

    # now calculate the relative difference
    np.seterr(divide='ignore', invalid='ignore')
    difference_in_sigma = np.nan_to_num((sigma_calc[0, [0, 5, 19], 0, 0] - sigma_ref) / sigma_ref)
    difference_in_seebeck = np.nan_to_num((seebeck_calc[0, [0, 5, 19], 0, 0] - seebeck_ref) / seebeck_ref)
    difference_in_lorenz = np.nan_to_num((lorenz_calc[0, [0, 5, 19], 0, 0] - lorenz_ref) / lorenz_ref)

    # should match down to numerical precision for all chemical
    # potentials
    assert np.all(np.abs(difference_in_sigma[:]) < 1e-7)
    assert np.all(np.abs(difference_in_seebeck[:]) < 1e-7)
    assert np.all(np.abs(difference_in_lorenz[:]) < 1e-7)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_parabolic_one_band'], indirect=True)
def test_trapz_parabolic_one_band(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
    """
    Test the direct integrals in reciprocal space against
    the results of the closed Fermi-Dirac integrals
    for a single valence band with a free electron mass,
    acoustic phonon scattering, at 100, 300, 600 K for
    chemical potential between -0.4 and 0.4 eV
    in 10 steps. The supplied grid is 61x61x61.

    """

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # now test, only first element along diagonal (xx)
    # at the supplied grid of 61x61x61 they all should
    # match to within 1e-2 (sigma) for 100 and 300 K and
    # 1e-3 for 700 K.

    # 100 K, only last 6 chempot, sigma and lorenz
    assert np.all(difference_in_sigma[0, 4:9, 0, 0] < 1e-2)
    assert np.all(difference_in_lorenz[0, 4:9, 0, 0] < 1e-2)

    # 300 K, full range of chempot, all
    assert np.all(difference_in_sigma[1, :, 0, 0] < 1e-2)
    assert np.all(difference_in_seebeck[1, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[1, :, 0, 0] < 1e-3)

    # 700 K, full range of  chempot, all
    assert np.all(difference_in_sigma[2, :, 0, 0] < 1e-3)
    assert np.all(difference_in_seebeck[2, :, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[2, :, 0, 0] < 1e-3)

    # 100 K last chempot entry (non-deg limit) for seebeck
    # and lorenz
    assert np.all(difference_in_seebeck[0, 9, 0, 0] < 1e-3)
    assert np.all(difference_in_lorenz[0, 9, 0, 0] < 1e-2)

    # 300 K first chempot (degen limit) sigma
    assert np.all(difference_in_sigma[1, 0, 0, 0] < 1e-4)

    # 700 K first chempot (deg limit) for all
    assert np.all(difference_in_sigma[2, 0, 0, 0] < 1e-5)
    assert np.all(difference_in_seebeck[2, 0, 0, 0] < 1e-5)
    assert np.all(difference_in_lorenz[2, 0, 0, 0] < 1e-5)


@pytest.mark.parametrize('read_and_setup_bs', ['tetra_parabolic_one_band'], indirect=True)
def test_tetra_parabolic_one_band(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
    """
    Test the direct integrals in reciprocal space against
    the results of the closed Fermi-Dirac integrals
    for a single valence band with a free electron mass,
    acoustic phonon scattering, at 100, 300, 600 K for
    chemical potential between -0.4 and 0.4 eV
    in 10 steps. The supplied grid is 61x61x61.

    """

    bs = read_and_setup_bs

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
    difference_in_sigma = np.abs((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.abs((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.abs((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # for the tetrahedron method we need very dense grids
    # so the accuracy is not so high for the input 61x61x61 grid
    # for the current implementation

    # five first chempots for sigma
    assert np.all(difference_in_sigma[:, 0:4, 0, 0] < 1e-1)
    # five first chempots for seebeck
    assert np.all(difference_in_seebeck[:, 0:4, 0, 0] < 1e-2)
    # five first chempots for lorenz
    assert np.all(difference_in_lorenz[:, 0:4, 0, 0] < 1e-2)


@pytest.mark.parametrize('read_and_setup_bs', ['trapz_preinter_parabolic_one_band'], indirect=True)
def test_trapz_preinter_parabolic_one_band(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
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

    bs = read_and_setup_bs

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
    difference_in_sigma = np.nan_to_num((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.nan_to_num((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.nan_to_num((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # this data should be reasonably converged, at least at
    # 300 and 600 K and close to a few percent at
    # 100 K.

    # 100 K
    assert np.all(np.abs(difference_in_sigma[0, :, 0, 0]) < 1e-2)
    # (only the last 5 chempots)
    assert np.all(np.abs(difference_in_seebeck[0, 4:9, 0, 0]) < 1e-1)
    assert np.all(np.abs(difference_in_lorenz[0, :, 0, 0]) < 1e-2)
    # 300 K
    assert np.all(np.abs(difference_in_sigma[1, :, 0, 0]) < 1e-2)
    assert np.all(np.abs(difference_in_seebeck[1, :, 0, 0]) < 1e-3)
    assert np.all(np.abs(difference_in_lorenz[1, :, 0, 0]) < 1e-3)
    # 700 K
    assert np.all(np.abs(difference_in_sigma[2, :, 0, 0]) < 1e-3)
    assert np.all(np.abs(difference_in_seebeck[2, :, 0, 0]) < 1e-3)
    assert np.all(np.abs(difference_in_lorenz[2, :, 0, 0]) < 1e-3)


@pytest.mark.parametrize('read_and_setup_bs', ['cubature_parabolic_one_band'], indirect=True)
def test_cubature_parabolic_one_band(read_and_setup_bs):
    gsl = pytest.importorskip("t4me.gsl")
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

    bs = read_and_setup_bs

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
    difference_in_sigma = np.nan_to_num((sigma_numerick - sigma_closed) / sigma_closed)
    difference_in_seebeck = np.nan_to_num((seebeck_numerick - seebeck_closed) / seebeck_closed)
    difference_in_lorenz = np.nan_to_num((lorenz_numerick - lorenz_closed) / lorenz_closed)

    # this data should be reasonably converged
    # 300 K
    assert np.all(np.abs(difference_in_sigma[0, :, 0, 0]) < 1e-2)
    assert np.all(np.abs(difference_in_seebeck[0, :, 0, 0]) < 1e-2)
    # only first 3 chempots
    assert np.all(np.abs(difference_in_lorenz[0, 0:2, 0, 0]) < 1e-2)
    # 700 K
    assert np.all(np.abs(difference_in_sigma[1, :, 0, 0]) < 1e-3)
    assert np.all(np.abs(difference_in_seebeck[1, :, 0, 0]) < 1e-2)
    # only first 5 chempots
    assert np.all(np.abs(difference_in_lorenz[1, 0:4, 0, 0]) < 1e-2)
