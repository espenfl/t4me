import os
import pytest
import t4me.inputoutput as inputoutput
import t4me.lattice as lattice

@pytest.fixture
def read_and_setup_lattice(request):
    """
    This routine reads the parameter files, sets up
    the lattice

    Returns
    -------
    lattice : object
        A `Lattice()` object.
    """

    # unfortunately we cannot pickle with classes
    # so load data everytime this is called
    location = os.path.dirname(__file__) + \
        "/../../../tests/test_data/" + str(request.param)
    # read parameters, generate lattice and bands
    param = inputoutput.Param(inputoutput.readparam(location=location))
    lat = lattice.Lattice(param, location=location)

    return lat
