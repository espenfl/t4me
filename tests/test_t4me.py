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
