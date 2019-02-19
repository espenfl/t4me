import os
import pytest
import numpy as np
import t4me.inputoutput as inputoutput
import t4me.lattice as lattice
import t4me.bandstructure as bandstructure


@pytest.fixture
def read_and_setup_bs(request):
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
        "/test_data/" + str(request.param)
    # read parameters, generate lattice and bands
    param = inputoutput.Param(inputoutput.readparam(location=location))
    lat = lattice.Lattice(param, location=location)
    bs = bandstructure.Bandstructure(lat, param, location=location)
    return bs


@pytest.fixture
def may_data():
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
    data[:, 1] = [432.97, 350.24, 272.50, 204.50, 150.88, 112.39, 86.09, 68.22, 46.95, 35.40, 28.35]
    # Lorenz coefficient in 10^-8 V^2/K^2
    data[:, 2] = [1.49, 1.51, 1.54, 1.61, 1.72, 1.86, 1.99, 2.09, 2.24, 2.32, 2.36]
    # carrier concentration in 10^21 cm^-3
    data[:, 3] = [0.00123, 0.00324, 0.0082, 0.0192, 0.0395, 0.0709, 0.1126, 0.1634, 0.2872, 0.4355, 0.6044]
    # unitless hall coefficient
    data[:, 4] = [1.17, 1.17, 1.16, 1.13, 1.11, 1.08, 1.05, 1.04, 1.02, 1.01, 1.01]
    return data
