#!/usr/local/bin/python

# pylint: disable=unused-import, redefined-outer-name

import pytest
from t4me.tests.fixtures import read_and_setup_lattice

@pytest.mark.parametrize('read_and_setup_lattice', ['si_45_primitive'], indirect=True)
def test_si_45_primitive(read_and_setup_lattice):
    """Test that Spglib returns the correct symmetry for the primitive cell of silicon."""

    lattice = read_and_setup_lattice
    assert lattice.spacegroup == 'Fd-3m (227)'
