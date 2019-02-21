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
"""Contains various input and output routines for T4ME."""

# pylint: disable=useless-import-alias, too-many-arguments, invalid-name, too-many-statements, too-many-lines, global-statement

import os
import sys
import logging
import numpy as np
import yaml

import t4me.constants as constants

# global variables
pythtb_warning_printed = False
skw_warning_printed = False
wildmagic_warning_printed = False
cubature_warning_printed = False
einspline_warning_printed = False
alglib_warning_printed = False
mpi4py_message_printed = False

# errors that cause exits


def spglib_error():
    """An error for missing Spglib interface."""
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.error("Could not locate the Spglib interface. Please "
                 "run 'python setup.py build_ext --inplace before "
                 "continuing. Exiting.")
    sys.exit(1)


# warnings


def gsl_warning():
    """An error for missing GNU GSL interface."""
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.warning("Could not locate the GSL interface. Please "
                   "do not use this functionality. Continuing.")


def alglib_warning():
    """An error for missing Alglib."""
    global alglib_warning_printed
    if not alglib_warning_printed:
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.warning("Could not load the Alglib package. Using "
                       "the Rbf function in SciPy instead. "
                       "This is extremely memory consuming, please "
                       "install the Alglib package properly of you "
                       "want to utilize RBF interpolation. "
                       "Continuing.")
        alglib_warning_printed = True


def pythtb_warning():
    """An error for missing PythTB."""
    global pythtb_warning_printed
    if not pythtb_warning_printed:
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.warning("Could not locate PythTB. PythTB calls will "
                       "give errors. Please do not set type: 3 in "
                       "bandparams.yml to generating TB bands "
                       "(use PythTB interface). Continuing.")
        pythtb_warning_printed = True


def skw_warning():
    """An error for missing SKW."""
    global skw_warning_printed
    if not skw_warning_printed:
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.warning("Could not locate the SKW interpolation interface. "
                       "Make sure you do not call any of its functions as "
                       "this will yield errors. Continuing.")
        skw_warning_printed = True


def wildmagic_warning():
    """An error for a missing GeometricTools interface."""
    global wildmagic_warning_printed
    if not wildmagic_warning_printed:
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.warning("Could not locate the WildMagic interpolation "
                       "interface. Make sure you do not call any of its "
                       "functions as this will yield errors. Continuing.")
        wildmagic_warning_printed = True


def cubature_warning():
    """An error for a missing Cubature-GeometricTools interface."""
    global cubature_warning_printed
    if not cubature_warning_printed:
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.warning("Could not locate the Cubature-WildMagic "
                       "interface. Make sure you do not call any of its "
                       "functions as this will yield errors. Continuing.")
        cubature_warning_printed = True


def einspline_warning():
    """An error for a missing Einspline interface"""
    global einspline_warning_printed
    if not einspline_warning_printed:
        # set logger
        logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
        logger.warning("Could not locate the Einspline interpolation "
                       "interface. Make sure you do not call any of its "
                       "functions as this will yield errors. Continuing.")
        einspline_warning_printed = True


class Param():  # pylint: disable=too-few-public-methods
    """
    YAML reader for the input paramters.

    Parameters
    ----------
    data : iterable
        yaml load (typically safe_load(open(yamlfilename),"r")).


    Notes
    -----
    Read a YAML paramter file.

    """

    def __init__(self, data):
        for name, value in list(data.items()):
            setattr(self, name, self._wrap(value))

    def _wrap(self, value):
        if isinstance(value, (tuple, list, set, frozenset)):
            return type(value)([self._wrap(v) for v in value])
        return Param(value) if isinstance(value, dict) else value


def readparam(location=None, filename=None):
    """
    Load the parameters in the general configuration file.

    Parameters
    ----------
    location : string, optional
        The location of the general configuration file.
        Defaults to the "input" directory in the
        current working directory.

    filename : string, optional
        The filename for the general configuration file.
        Defaults to "param.yml".

    Returns
    -------
    iterable
        An iterable YAML object.

    Notes
    -----
    The current working directory is padded in front of
    any supplied location (or if path is given in the filename).

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running readparam.")

    if filename is None:
        filename = "param.yml"
    if location is not None:
        filename = location + "/" + filename
    else:
        filename = os.getcwd() + "/input/" + filename
    try:
        stream = open(filename, "r")
    except IOError:
        logger.error("Can not open %s. Exiting.", filename)
        sys.exit(1)
    return yaml.safe_load(stream)


def readbandparam(location=None, filename=None):
    """
    Load the parameters in the bandstructure configuration file.

    Parameters
    ----------
    location : string, optional
        The location of the bandstructure configuration file.
        Defaults to "input" directory in the current working
        directory.

    filename : string, optional
        The filename for the bandstructure configuration file.
        Defaults to "bandparam.yml".

    Returns
    -------
    iterable
        An iterable YAML object.

    Notes
    -----
    The current working directory is padded in front of
    any supplied location (or if path is given in the filename).

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running readbandparam.")

    if filename is None:
        filename = "bandparam.yml"
    if location is not None:
        filename = location + "/" + filename
    else:
        filename = os.getcwd() + "/input/" + filename
    try:
        stream = open(filename, "r")
    except IOError:
        logger.error("Can not open %s. Exiting.", filename)
        sys.exit(1)
    return yaml.safe_load(stream)


def readcellparam(location=None, filename=None):
    """
    Load the parameters in the cell configuration file.

    Parameters
    ----------
    location : string, optional
        The location of the cell configuration file.
        Defaults to the "input" directory in the
        current working directory.

    filename : string, optional
        The filename for the cell configuration file.
        Defaults to "cellparam.yml".

    Returns
    -------
    iterable
        An iterable YAML object.

    Notes
    -----
    The current working directory is padded in front of
    any supplied location (or if path is given in the filename).

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running readcellparam.")

    if filename is None:
        filename = "cellparam.yml"
    if location is not None:
        filename = location + "/" + filename
    else:
        filename = os.getcwd() + "/input/" + filename
    try:
        stream = open(filename, "r")
    except IOError:
        logger.error("Can not open %s. Exiting.", filename)
        sys.exit(1)
    return yaml.safe_load(stream)


def dump_transport_coefficients(tr, filename_tag=None):  # pylint: disable=too-many-locals
    r"""
    Writes the transport coefficients to files

    Parameters
    ----------
    tr : object
        A `Transport()` object that contains the transport coefficients
    filename_tag : string, optional
        | If `filename_tag` is not an empty string, but a string `x`, the
        | output filenames are:
        | sigma_x: contains the electrical conductivity in units of
        | :math:`\\mathrm{S}/\\mathrm{m}`
        | seebeck_x:  contains the Seebeck coefficiens in units of
        | :math:`\\mu \\mathrm{V}/\\mathrm{K}`
        | lorenz_x:  contains the Lorenz number in units of
        | :math:`10^{-8} \\mathrm{V^2}/\\mathrm{K^2}`
        | kappa_x:  contains the Seebeck coefficiens in units of
        | :math:`\\mu \\mathrm{W}/\\mathrm{mK}`
        | hall_x:  the Hall coefficient (big R) in units of
        | :math:`\\mathrm{cm^{3}}/\\mathrm{C}`
        | cc_x: the carrier concentration in units of
        | :math:`10^{21} \\mathrm{cm^{-3}}`

        The default is to write the files without the tag on the end.
        Consult header of the files for the ordering.

    Returns
    -------
    None

    Notes
    -----
    Each temperature steps have its own block.

    """
    # loaded here in order to avoid chained imported with
    # lbtecoeff/inputout for printout messages, this should be fixed
    # in the future as it is ugly
    import t4me.lbtecoeff as lbtecoeff

    if filename_tag is None:
        filename_tag = ""

    # open files
    sigma_file = file_handler(filename="output/sigma", status="w")
    seebeck_file = file_handler(filename="output/seebeck", status="w")
    lorenz_file = file_handler(filename="output/lorenz", status="w")
    kappae_file = file_handler(filename="output/kappae", status="w")
    hall_file = file_handler(filename="output/hall", status="w")
    hall_cc_file = file_handler(filename="output/hall_cc", status="w")
    cc_file = file_handler(filename="output/cc", status="w")
    hall_fact_file = file_handler(filename="output/hall_fact", status="w")

    # write headers
    # sigma
    sigma_file.write(
        "##########################################################\n")
    sigma_file.write(
        "##########################################################\n")
    sigma_file.write(
        "#                                                        #\n")
    sigma_file.write(
        "##########################################################\n")
    sigma_file.write(
        "#                                                        #\n")
    sigma_file.write(
        "#                Electrical conductivity                 #\n")
    sigma_file.write(
        "#                                                        #\n")
    sigma_file.write(
        "##########################################################\n")
    sigma_file.write(
        "# ------------------------------------------------------ #\n")
    sigma_file.write(
        "# ------------------------------------------------------ #\n")
    sigma_file.write(
        "#                                                        #\n")
    sigma_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    sigma_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    sigma_file.write(
        "#     bands)                                             #\n")
    sigma_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    sigma_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    sigma_file.write(
        "#  5. Column: electrical conductivity in S/m (xx tensor) #\n")
    sigma_file.write(
        "#  6. Column: electrical conductivity in S/m (yy tensor) #\n")
    sigma_file.write(
        "#  7. Column: electrical conductivity in S/m (zz tensor) #\n")
    sigma_file.write(
        "#  8. Column: electrical conductivity in S/m (xy tensor) #\n")
    sigma_file.write(
        "#  8. Column: electrical conductivity in S/m (xz tensor) #\n")
    sigma_file.write(
        "#  9. Column: electrical conductivity in S/m (yx tensor) #\n")
    sigma_file.write(
        "# 10. Column: electrical conductivity in S/m (yz tensor) #\n")
    sigma_file.write(
        "# 11. Column: electrical conductivity in S/m (zx tensor) #\n")
    sigma_file.write(
        "# 12. Column: electrical conductivity in S/m (zy tensor) #\n")
    sigma_file.write(
        "#                                                        #\n")
    sigma_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    sigma_file.write(
        "#                                                        #\n")
    sigma_file.write(
        "# ------------------------------------------------------ #\n")
    sigma_file.write(
        "# ------------------------------------------------------ #\n")
    sigma_file.write(
        "##########################################################\n")

    # Seebeck
    seebeck_file.write(
        "##########################################################\n")
    seebeck_file.write(
        "##########################################################\n")
    seebeck_file.write(
        "#                                                        #\n")
    seebeck_file.write(
        "##########################################################\n")
    seebeck_file.write(
        "#                                                        #\n")
    seebeck_file.write(
        "#                 Seebeck coefficient                    #\n")
    seebeck_file.write(
        "#                                                        #\n")
    seebeck_file.write(
        "##########################################################\n")
    seebeck_file.write(
        "# ------------------------------------------------------ #\n")
    seebeck_file.write(
        "# ------------------------------------------------------ #\n")
    seebeck_file.write(
        "#                                                        #\n")
    seebeck_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    seebeck_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    seebeck_file.write(
        "#     bands)                                             #\n")
    seebeck_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    seebeck_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    seebeck_file.write(
        "#  5. Column: Seebeck coefficient in muV/K (xx tensor)   #\n")
    seebeck_file.write(
        "#  6. Column: Seebeck coefficient in muV/K (yy tensor)   #\n")
    seebeck_file.write(
        "#  7. Column: Seebeck coefficient in muV/K (zz tensor)   #\n")
    seebeck_file.write(
        "#  8. Column: Seebeck coefficient in muV/K (xy tensor)   #\n")
    seebeck_file.write(
        "#  8. Column: Seebeck coefficient in muV/K (xz tensor)   #\n")
    seebeck_file.write(
        "#  9. Column: Seebeck coefficient in muV/K (yx tensor)   #\n")
    seebeck_file.write(
        "# 10. Column: Seebeck coefficient in muV/K (yz tensor)   #\n")
    seebeck_file.write(
        "# 11. Column: Seebeck coefficient in muV/K (zx tensor)   #\n")
    seebeck_file.write(
        "# 12. Column: Seebeck coefficient in muV/K (zy tensor)   #\n")
    seebeck_file.write(
        "#                                                        #\n")
    seebeck_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    seebeck_file.write(
        "#                                                        #\n")
    seebeck_file.write(
        "# ------------------------------------------------------ #\n")
    seebeck_file.write(
        "# ------------------------------------------------------ #\n")
    seebeck_file.write(
        "##########################################################\n")

    # Lorenz
    lorenz_file.write(
        "##########################################################\n")
    lorenz_file.write(
        "##########################################################\n")
    lorenz_file.write(
        "#                                                        #\n")
    lorenz_file.write(
        "##########################################################\n")
    lorenz_file.write(
        "#                                                        #\n")
    lorenz_file.write(
        "#                  Lorenz coefficient                    #\n")
    lorenz_file.write(
        "#                                                        #\n")
    lorenz_file.write(
        "##########################################################\n")
    lorenz_file.write(
        "# ------------------------------------------------------ #\n")
    lorenz_file.write(
        "# ------------------------------------------------------ #\n")
    lorenz_file.write(
        "#                                                        #\n")
    lorenz_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    lorenz_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    lorenz_file.write(
        "#     bands)                                             #\n")
    lorenz_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    lorenz_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    lorenz_file.write(
        "#  5. Column: Lorenz coeff. in 10^-8 V^2/K^2 (xx tensor) #\n")
    lorenz_file.write(
        "#  6. Column: Lorenz coeff. in 10^-8 V^2/K^2 (yy tensor) #\n")
    lorenz_file.write(
        "#  7. Column: Lorenz coeff. in 10^-8 V^2/K^2 (zz tensor) #\n")
    lorenz_file.write(
        "#  8. Column: Lorenz coeff. in 10^-8 V^2/K^2 (xy tensor) #\n")
    lorenz_file.write(
        "#  8. Column: Lorenz coeff. in 10^-8 V^2/K^2 (xz tensor) #\n")
    lorenz_file.write(
        "#  9. Column: Lorenz coeff. in 10^-8 V^2/K^2 (yx tensor) #\n")
    lorenz_file.write(
        "# 10. Column: Lorenz coeff. in 10^-8 V^2/K^2 (yz tensor) #\n")
    lorenz_file.write(
        "# 11. Column: Lorenz coeff. in 10^-8 V^2/K^2 (zx tensor) #\n")
    lorenz_file.write(
        "# 12. Column: Lorenz coeff. in 10^-8 V^2/K^2 (zy tensor) #\n")
    lorenz_file.write(
        "#                                                        #\n")
    lorenz_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    lorenz_file.write(
        "#                                                        #\n")
    lorenz_file.write(
        "# ------------------------------------------------------ #\n")
    lorenz_file.write(
        "# ------------------------------------------------------ #\n")
    lorenz_file.write(
        "##########################################################\n")

    # kappae
    kappae_file.write(
        "##########################################################\n")
    kappae_file.write(
        "##########################################################\n")
    kappae_file.write(
        "#                                                        #\n")
    kappae_file.write(
        "##########################################################\n")
    kappae_file.write(
        "#                                                        #\n")
    kappae_file.write(
        "#       Electron part of the thermal conductivity        #\n")
    kappae_file.write(
        "#                                                        #\n")
    kappae_file.write(
        "##########################################################\n")
    kappae_file.write(
        "# ------------------------------------------------------ #\n")
    kappae_file.write(
        "# ------------------------------------------------------ #\n")
    kappae_file.write(
        "#                                                        #\n")
    kappae_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    kappae_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    kappae_file.write(
        "#     bands)                                             #\n")
    kappae_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    kappae_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    kappae_file.write(
        "#  5. Column: el.therm. conductivity in W/mK (xx tensor) #\n")
    kappae_file.write(
        "#  6. Column: el.therm. conductivity in W/mK (yy tensor) #\n")
    kappae_file.write(
        "#  7. Column: el.therm. conductivity in W/mK (zz tensor) #\n")
    kappae_file.write(
        "#  8. Column: el.therm. conductivity in W/mK (xy tensor) #\n")
    kappae_file.write(
        "#  8. Column: el.therm. conductivity in W/mK (xz tensor) #\n")
    kappae_file.write(
        "#  9. Column: el.therm. conductivity in W/mK (yx tensor) #\n")
    kappae_file.write(
        "# 10. Column: el.therm. conductivity in W/mK (yz tensor) #\n")
    kappae_file.write(
        "# 11. Column: el.therm. conductivity in W/mK (zx tensor) #\n")
    kappae_file.write(
        "# 12. Column: el.therm. conductivity in W/mK (zy tensor) #\n")
    kappae_file.write(
        "#                                                        #\n")
    kappae_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    kappae_file.write(
        "#                                                        #\n")
    kappae_file.write(
        "# ------------------------------------------------------ #\n")
    kappae_file.write(
        "# ------------------------------------------------------ #\n")
    kappae_file.write(
        "##########################################################\n")

    # Hall coefficient
    hall_file.write(
        "##########################################################\n")
    hall_file.write(
        "##########################################################\n")
    hall_file.write(
        "#                                                        #\n")
    hall_file.write(
        "##########################################################\n")
    hall_file.write(
        "#                                                        #\n")
    hall_file.write(
        "#               Hall coefficient (big R)                 #\n")
    hall_file.write(
        "#                                                        #\n")
    hall_file.write(
        "##########################################################\n")
    hall_file.write(
        "# ------------------------------------------------------ #\n")
    hall_file.write(
        "# ------------------------------------------------------ #\n")
    hall_file.write(
        "#                                                        #\n")
    hall_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    hall_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    hall_file.write(
        "#     bands)                                             #\n")
    hall_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    hall_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    hall_file.write(
        "#  5. Column: Hall coefficient in cm^3/C (xx tensor)     #\n")
    hall_file.write(
        "#  6. Column: Hall coefficient in cm^3/C (yy tensor)     #\n")
    hall_file.write(
        "#  7. Column: Hall coefficient in cm^3/C (zz tensor)     #\n")
    hall_file.write(
        "#  8. Column: Hall coefficient in cm^3/C (xy tensor)     #\n")
    hall_file.write(
        "#  8. Column: Hall coefficient in cm^3/C (xz tensor)     #\n")
    hall_file.write(
        "#  9. Column: Hall coefficient in cm^3/C (yx tensor)     #\n")
    hall_file.write(
        "# 10. Column: Hall coefficient in cm^3/C (yz tensor)     #\n")
    hall_file.write(
        "# 11. Column: Hall coefficient in cm^3/C (zx tensor)     #\n")
    hall_file.write(
        "# 12. Column: Hall coefficient in cm^3/C (zy tensor)     #\n")
    hall_file.write(
        "#                                                        #\n")
    hall_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    hall_file.write(
        "#                                                        #\n")
    hall_file.write(
        "# ------------------------------------------------------ #\n")
    hall_file.write(
        "# ------------------------------------------------------ #\n")
    hall_file.write(
        "##########################################################\n")

    # Hall carrier concentration
    hall_cc_file.write(
        "##########################################################\n")
    hall_cc_file.write(
        "##########################################################\n")
    hall_cc_file.write(
        "#                                                        #\n")
    hall_cc_file.write(
        "##########################################################\n")
    hall_cc_file.write(
        "#                                                        #\n")
    hall_cc_file.write(
        "#              Hall carrier concentration                #\n")
    hall_cc_file.write(
        "#                                                        #\n")
    hall_cc_file.write(
        "##########################################################\n")
    hall_cc_file.write(
        "# ------------------------------------------------------ #\n")
    hall_cc_file.write(
        "# ------------------------------------------------------ #\n")
    hall_cc_file.write(
        "#                                                        #\n")
    hall_cc_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    hall_cc_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    hall_cc_file.write(
        "#     bands)                                             #\n")
    hall_cc_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    hall_cc_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    hall_cc_file.write(
        "#  5. Column: Hall c.conc. in 10^21 cm^-3 (xx tensor)    #\n")
    hall_cc_file.write(
        "#  6. Column: Hall c.conc. in 10^21 cm^-3 (yy tensor)    #\n")
    hall_cc_file.write(
        "#  7. Column: Hall c.conc. in 10^21 cm^-3 (zz tensor)    #\n")
    hall_cc_file.write(
        "#  8. Column: Hall c.conc. in 10^21 cm^-3 (xy tensor)    #\n")
    hall_cc_file.write(
        "#  8. Column: Hall c.conc. in 10^21 cm^-3 (xz tensor)    #\n")
    hall_cc_file.write(
        "#  9. Column: Hall c.conc. in 10^21 cm^-3 (yx tensor)    #\n")
    hall_cc_file.write(
        "# 10. Column: Hall c.conc. in 10^21 cm^-3 (yz tensor)    #\n")
    hall_cc_file.write(
        "# 11. Column: Hall c.conc. in 10^21 cm^-3 (zx tensor)    #\n")
    hall_cc_file.write(
        "# 12. Column: Hall c.conc. in 10^21 cm^-3 (zy tensor)    #\n")
    hall_cc_file.write(
        "#                                                        #\n")
    hall_cc_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    hall_cc_file.write(
        "#                                                        #\n")
    hall_cc_file.write(
        "# ------------------------------------------------------ #\n")
    hall_cc_file.write(
        "# ------------------------------------------------------ #\n")
    hall_cc_file.write(
        "##########################################################\n")

    # carrier concentration
    cc_file.write(
        "##########################################################\n")
    cc_file.write(
        "##########################################################\n")
    cc_file.write(
        "#                                                        #\n")
    cc_file.write(
        "##########################################################\n")
    cc_file.write(
        "#                                                        #\n")
    cc_file.write(
        "#                 Carrier concentration                  #\n")
    cc_file.write(
        "#                                                        #\n")
    cc_file.write(
        "##########################################################\n")
    cc_file.write(
        "# ------------------------------------------------------ #\n")
    cc_file.write(
        "# ------------------------------------------------------ #\n")
    cc_file.write(
        "#                                                        #\n")
    cc_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    cc_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    cc_file.write(
        "#     bands)                                             #\n")
    cc_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    cc_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    cc_file.write(
        "#                                                        #\n")
    cc_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    cc_file.write(
        "#                                                        #\n")
    cc_file.write(
        "# ------------------------------------------------------ #\n")
    cc_file.write(
        "# ------------------------------------------------------ #\n")
    cc_file.write(
        "##########################################################\n")

    # Hall factor
    hall_fact_file.write(
        "##########################################################\n")
    hall_fact_file.write(
        "##########################################################\n")
    hall_fact_file.write(
        "#                                                        #\n")
    hall_fact_file.write(
        "##########################################################\n")
    hall_fact_file.write(
        "#                                                        #\n")
    hall_fact_file.write(
        "#                      Hall factor                       #\n")
    hall_fact_file.write(
        "#                                                        #\n")
    hall_fact_file.write(
        "##########################################################\n")
    hall_fact_file.write(
        "# ------------------------------------------------------ #\n")
    hall_fact_file.write(
        "# ------------------------------------------------------ #\n")
    hall_fact_file.write(
        "#                                                        #\n")
    hall_fact_file.write(
        "#  1. Column: chemical potential in eV                   #\n")
    hall_fact_file.write(
        "#  2. Column: eta (not usefull for other than parabolic  #\n")
    hall_fact_file.write(
        "#     bands)                                             #\n")
    hall_fact_file.write(
        "#  3. Column: n-type carrier concen. in 10^21 cm^-3      #\n")
    hall_fact_file.write(
        "#  4. Column: p-type carrier concen. in 10^21 cm^-3      #\n")
    hall_fact_file.write(
        "#  5. Column: Hall factor (xx tensor) n-type             #\n")
    hall_fact_file.write(
        "#  6. Column: Hall factor (yy tensor) n-type             #\n")
    hall_fact_file.write(
        "#  7. Column: Hall factor (zz tensor) n-type             #\n")
    hall_fact_file.write(
        "#  8. Column: Hall factor (xy tensor) n-type             #\n")
    hall_fact_file.write(
        "#  8. Column: Hall factor (xz tensor) n-type             #\n")
    hall_fact_file.write(
        "#  9. Column: Hall factor (yx tensor) n-type             #\n")
    hall_fact_file.write(
        "# 10. Column: Hall factor (yz tensor) n-type             #\n")
    hall_fact_file.write(
        "# 11. Column: Hall factor (zx tensor) n-type             #\n")
    hall_fact_file.write(
        "# 12. Column: Hall factor (zy tensor) n-type             #\n")
    hall_fact_file.write(
        "# 13-. Column: Similar but for p-type                    #\n")
    hall_fact_file.write(
        "#                                                        #\n")
    hall_fact_file.write(
        "#  Data is blocked on temperature in K                   #\n")
    hall_fact_file.write(
        "#                                                        #\n")
    hall_fact_file.write(
        "# ------------------------------------------------------ #\n")
    hall_fact_file.write(
        "# ------------------------------------------------------ #\n")
    hall_fact_file.write(
        "##########################################################\n")

    # extract the Hall carrier concentration
    # TODO: NEED TO IMPLEMENT THIS FOR ALL METHODS pylint: disable=fixme
    if not tr.param.transport_method == "numerick":
        hall_cc = lbtecoeff.calculate_hall_carrier_concentration(tr.hall)
        # extract the Hall factor
        hall_factorn = lbtecoeff.calculate_hall_factor(tr.ccn, hall_cc)
        hall_factorp = lbtecoeff.calculate_hall_factor(tr.ccp, hall_cc)
    else:
        hall_cc = np.zeros(tr.sigma.shape)
        hall_factorn = np.zeros(tr.sigma.shape)
        hall_factorp = np.zeros(tr.sigma.shape)

    for indext, temp in np.ndenumerate(tr.temperature):
        sigma_file.write("# temperature: " + str(temp) + "\n")
        seebeck_file.write("# temperature: " + str(temp) + "\n")
        lorenz_file.write("# temperature: " + str(temp) + "\n")
        kappae_file.write("# temperature: " + str(temp) + "\n")
        hall_file.write("# temperature: " + str(temp) + "\n")
        hall_cc_file.write("# temperature: " + str(temp) + "\n")
        cc_file.write("# temperature: " + str(temp) + "\n")
        hall_fact_file.write("# temperature: " + str(temp) + "\n")
        for indexe, chempot in np.ndenumerate(tr.chempots):
            t = indext[0]
            e = indexe[0]
            # calculate kappae from sigma and Lorenz
            # also add temperature dependence and the 10^-8 factor
            # in the Lorenz
            kappae = np.dot(tr.lorenz[t, e], tr.sigma[t, e]) * temp * 1e-8
            eta = 1e5 * chempot / (constants.kb * temp)
            sigma_file.write("{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                             "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                             "{:>17.8e}{:>17.8e}{:>17.8e}\n".format(
                                 chempot, eta, tr.ccn[t, e, 0, 0],
                                 tr.ccp[t, e, 0, 0], tr.sigma[t, e, 0, 0],
                                 tr.sigma[t, e, 1, 1], tr.sigma[t, e, 2, 2],
                                 tr.sigma[t, e, 0, 1], tr.sigma[t, e, 0, 2],
                                 tr.sigma[t, e, 1, 0], tr.sigma[t, e, 1, 2],
                                 tr.sigma[t, e, 2, 0], tr.sigma[t, e, 2, 1]))
            seebeck_file.write(
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}\n".format(
                    chempot, eta, tr.ccn[t, e, 0, 0], tr.ccp[t, e, 0, 0],
                    tr.seebeck[t, e, 0, 0], tr.seebeck[t, e, 1, 1],
                    tr.seebeck[t, e, 2, 2], tr.seebeck[t, e, 0, 1],
                    tr.seebeck[t, e, 0, 2], tr.seebeck[t, e, 1, 0],
                    tr.seebeck[t, e, 1, 2], tr.seebeck[t, e, 2, 0],
                    tr.seebeck[t, e, 2, 1]))
            lorenz_file.write(
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}\n".format(
                    chempot, eta, tr.ccn[t, e, 0, 0], tr.ccp[t, e, 0, 0],
                    tr.lorenz[t, e, 0, 0], tr.lorenz[t, e, 1, 1],
                    tr.lorenz[t, e, 2, 2], tr.lorenz[t, e, 0, 1],
                    tr.lorenz[t, e, 0, 2], tr.lorenz[t, e, 1, 0],
                    tr.lorenz[t, e, 1, 2], tr.lorenz[t, e, 2, 0],
                    tr.lorenz[t, e, 2, 1]))
            kappae_file.write(
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}\n".format(chempot, eta, tr.ccn[t, e, 0, 0],
                                     tr.ccp[t, e, 0, 0], kappae[0, 0],
                                     kappae[1, 1], kappae[2, 2], kappae[0, 1],
                                     kappae[0, 2], kappae[1, 0], kappae[1, 2],
                                     kappae[2, 0], kappae[2, 1]))
            hall_file.write(
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}\n".format(chempot, eta, tr.ccn[t, e, 0, 0],
                                     tr.ccp[t, e, 0, 0], tr.hall[t, e, 0, 0],
                                     tr.hall[t, e, 1, 1], tr.hall[t, e, 2, 2],
                                     tr.hall[t, e, 0, 1], tr.hall[t, e, 0, 2],
                                     tr.hall[t, e, 1, 0], tr.hall[t, e, 1, 2],
                                     tr.hall[t, e, 2, 0], tr.hall[t, e, 2, 1]))
            hall_cc_file.write(
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}\n".format(chempot, eta, tr.ccn[t, e, 0, 0],
                                     tr.ccp[t, e, 0, 0], hall_cc[t, e, 0, 0],
                                     hall_cc[t, e, 1, 1], hall_cc[t, e, 2, 2],
                                     hall_cc[t, e, 0, 1], hall_cc[t, e, 0, 2],
                                     hall_cc[t, e, 1, 0], hall_cc[t, e, 1, 2],
                                     hall_cc[t, e, 2, 0], hall_cc[t, e, 2, 1]))
            cc_file.write("{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}\n".format(
                chempot, eta, tr.ccn[t, e, 0, 0], tr.ccp[t, e, 0, 0]))
            hall_fact_file.write(
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}{:>17.8e}{:>17.8e}"
                "{:>17.8e}{:>17.8e}\n".format(
                    chempot, eta, tr.ccn[t, e, 0, 0], tr.ccp[t, e, 0, 0],
                    hall_factorn[t, e, 0, 0], hall_factorn[t, e, 1, 1],
                    hall_factorn[t, e, 2, 2], hall_factorn[t, e, 0, 1],
                    hall_factorn[t, e, 0, 2], hall_factorn[t, e, 1, 0],
                    hall_factorn[t, e, 1, 2], hall_factorn[t, e, 2, 0],
                    hall_factorn[t, e, 2, 1], hall_factorp[t, e, 0, 0],
                    hall_factorp[t, e, 1, 1], hall_factorp[t, e, 2, 2],
                    hall_factorp[t, e, 0, 1], hall_factorp[t, e, 0, 2],
                    hall_factorp[t, e, 1, 0], hall_factorp[t, e, 1, 2],
                    hall_factorp[t, e, 2, 0], hall_factorp[t, e, 2, 1]))

        # two clear lines between blocks
        sigma_file.write("\n\n")
        seebeck_file.write("\n\n")
        lorenz_file.write("\n\n")
        kappae_file.write("\n\n")
        hall_file.write("\n\n")
        hall_cc_file.write("\n\n")
        cc_file.write("\n\n")
        hall_fact_file.write("\n\n")

    # close files
    file_handler(handler=sigma_file)
    file_handler(handler=seebeck_file)
    file_handler(handler=lorenz_file)
    file_handler(handler=kappae_file)
    file_handler(handler=hall_file)
    file_handler(handler=hall_cc_file)
    file_handler(handler=cc_file)
    file_handler(handler=hall_fact_file)


def dump_relaxation_time(tr, filename=None):  # pylint: disable=too-many-branches
    """
    Writes the relaxation time to file.

    Parameters
    ----------
    tr : object
        A `Transport()` object containing the relaxation time and
        other details related to the carrier transport.
    filename : string, optional
        The output filename, default is "scattering". The string
        "_band_N" is added to this string, where N is the band number.

    Returns
    -------
    None

    Notes
    -----
    One file per band, filename "scattering_band_x", where
    x is the band number. In each file the temperature dependence is
    blocked, while the carrier energy, total relaxation time
    and each individual relaxation times follow as columns for each block.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running dump_relaxation_time.")

    if filename is None:
        filename = "scattering"
    precalc_scattering = True
    # write one file for each band and a block for each temperature
    try:
        tr.scattering_inv
    except AttributeError:
        precalc_scattering = False
    # in order not to waste to much space the scattering_energies
    # is always static and similar for all bands unless we have
    # called the interpolate routine to lay the tau values
    # out on all E(k)
    noband_dep = True
    if len(tr.scattering_energies.shape) > 1:
        noband_dep = False
    if precalc_scattering:  # pylint: disable=too-many-nested-blocks
        for band in range(tr.scattering_total_inv[0].shape[0]):
            filename_band = filename + "_band_" + str(band + 1)
            # open
            scattering_file = file_handler(
                "output/" + filename_band, status='w')
            scattering_file.write(
                "##########################################################\n")
            scattering_file.write(
                "##########################################################\n")
            scattering_file.write(
                "#                                                        #\n")
            scattering_file.write(
                "##########################################################\n")
            scattering_file.write(
                "#                                                        #\n")
            scattering_file.write(
                "#                   Relaxation times                     #\n")
            scattering_file.write(
                "#                                                        #\n")
            scattering_file.write(
                "##########################################################\n")
            scattering_file.write(
                "# ------------------------------------------------------ #\n")
            scattering_file.write(
                "# ------------------------------------------------------ #\n")
            scattering_file.write(
                "#                                                        #\n")
            scattering_file.write(
                "#  1. Column: energy in eV                               #\n")
            scattering_file.write(
                "#  2. Column: total relaxtion time from Matthiessen's    #\n")
            scattering_file.write(
                "#     rule in fs                                         #\n")
            if not tr.param.onlytotalrate:
                scattering_file.write(
                    "#  3-14. Column: relaxation time in fs for each          #\n"
                )
                scattering_file.write(
                    "#        mechanism defined in the band parameter file    #\n"
                )
            scattering_file.write(
                "#                                                        #\n")
            scattering_file.write(
                "# ------------------------------------------------------ #\n")
            scattering_file.write(
                "# ------------------------------------------------------ #\n")
            scattering_file.write(
                "##########################################################\n")

            if noband_dep:
                scattering_energies = tr.scattering_energies
            else:
                scattering_energies = tr.scattering_energies[band]
            for temp in range(tr.scattering_total_inv.shape[0]):
                for energy in range(scattering_energies.shape[0]):
                    scattering_file.write(
                        str(scattering_energies[energy]) + " " +
                        str(tr.scattering_total_inv[temp, band, energy]))
                    if not tr.param.onlytotalrate:
                        for scatter in range(tr.scattering_inv.shape[3]):
                            scattering_file.write(" " + str(
                                tr.scattering_inv[temp, band, energy, scatter])
                                                  )
                    scattering_file.write("\n")
                # write two blank lines in order to block on temperature
                scattering_file.write("\n\n")
            # close
            file_handler(filename_band, scattering_file)
    else:
        logger.info("Did not find any precalulcated scattering and there is "
                    "no scattering array to write. However, we might "
                    "have used analytical models generated on the fly "
                    "instead.")


def dump_density_of_states(bs, dos=None, dos_energies=None, filename="dos"):
    """
    Writes the density of states to file.

    Parameters
    ----------
    bs : object
        A `Bandstructure()` object.
    dos : ndarray, optional
        | Dimension: (N,M)

        The density of states for N bands at M energy samplings
        If not supplied, set to `bs.dos`.
    dos_energies : ndarray, optional
        | Dimension: (M)

        The M energy samples used for the density of state
        If not supplied, set to `bs.dos_energies`
    filename : string, optional
        The filename used to write the density of states.
        Default is "dos".

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running dump_density_of_states.")

    if dos is None:
        dos = bs.dos
    if dos is None:
        logger.info("No entries for the density of states have been "
                    "located. If you want to dump the density of states "
                    "please enable its calculation.")
        return
    if dos_energies is None:
        dos_energies = bs.dos_energies

    # concatenate density of states for all bands, expand later and include
    # decomposed dos
    dos_file = file_handler("output/" + filename, status='w')
    dos_file.write(
        "##########################################################\n")
    dos_file.write(
        "##########################################################\n")
    dos_file.write(
        "#                                                        #\n")
    dos_file.write(
        "##########################################################\n")
    dos_file.write(
        "#                                                        #\n")
    dos_file.write(
        "#                  Density of states                     #\n")
    dos_file.write(
        "#                                                        #\n")
    dos_file.write(
        "##########################################################\n")
    dos_file.write(
        "# ------------------------------------------------------ #\n")
    dos_file.write(
        "# ------------------------------------------------------ #\n")
    dos_file.write(
        "#                                                        #\n")
    dos_file.write(
        "#  1. Column: energy in eV                               #\n")
    dos_file.write(
        "#  2. Column: density of states in 1/eV/AA^3             #\n")
    dos_file.write(
        "#                                                        #\n")
    dos_file.write(
        "# ------------------------------------------------------ #\n")
    dos_file.write(
        "# ------------------------------------------------------ #\n")
    dos_file.write(
        "##########################################################\n")
    for index, energy in np.ndenumerate(dos_energies):
        dos_file.write("{:>12.4e}".format(energy) +
                       "{:>12.4e}\n".format(dos[index]))

    file_handler(filename, dos_file)


def dump_bandstruct_line(bs,
                         kstart,
                         kend,
                         filename="band",
                         datatype="e",
                         k_direct=True,
                         itype=None,
                         itype_sub=None):
    """
    Writes the energy or velocity dispersions extracted along a line to a file.

    Parameters
    ----------
    bs : object
        A `Band()` object containing the energies and velocity dispersions.
    kstart : ndarray
        | Dimension: (3)

        The start k-point vector in cartesian coordinates.
    kend : ndarray
        | Dimension: (3)

        The end k-point vector in cartesian coordinates.
    filename : string, optional
        The filename used to write the energy or velocity dispersions.
        Defaults to "band".
    datatype : {"e","v"}
        Selects to write energy dispersions ("e") or velocity dispersions ("v").
    itype : string, optional
        | Can be any of:
        | {"linearnd", "interpn", "rbf", "einspline", "wildmagic", "skw"}

        The type of interpolate method to use. If not set, the parameter
        `dispersion_interpolate_method` in param.yml sets this.
    itype_sub : string, optional
        | Can be any of:
        | {"nearest", "linear"}, when `itype` is set to `interpn`.
        | {"multiquadric", "inverse_multiquadric", "gaussian", "linear",
        | "cubic", "quintic", "thin_plate"}, when `itype` is set to `rbf`
        | and when the Scipy variety is used (the `alglib` variable set
        | to False in the :func:`interpolate` function). If `alglib` is
        | set to True (default), then `itype_sub` does not have to be set.
        | {"natural", "flat", "periodic", "antiperiodic"}, when `itype`
        | is set to `einspline`.
        | {"trilinear, tricubic_exact, tricubic_bspline, akima"},
        | when `itype` is set to `wildmagic`.

        The subtype of the interpolation method.

    Returns
    -------
    None

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running dump_bandstruct_line.")

    # kstart and kend supplied in cartesian (next routine expects direct),
    # convert
    if not k_direct:
        ks = bs.lattice.cart_to_dir(kstart)
        ke = bs.lattice.cart_to_dir(kend)
    else:
        ks = kstart
        ke = kend

    if datatype == "e":
        # fetch energies and kpts
        e, kpts = bs.fetch_energies_along_line(
            ks, ke, itype=itype, itype_sub=itype_sub)
        kpoint_length = np.linalg.norm(kpts - kpts[0], axis=1)
        bands_file = file_handler("output/" + filename, status='w')
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "#                 Bandstructure lines                    #\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "#  1. Column: k-point distance in AA^-1                  #\n")
        bands_file.write(
            "#  2. Column: electron energy dispersion for the first   #\n")
        bands_file.write(
            "#     band in eV                                         #\n")
        bands_file.write(
            "#  3. Column: electron energy dispersion for the second  #\n")
        bands_file.write(
            "#     band in eV                                         #\n")
        bands_file.write(
            "#  n. Column: electron energy dispersion for the (n-1)   #\n")
        bands_file.write(
            "#     band in eV                                         #\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write("{:<7s}".format("Dist."))
        for band in range(e.shape[0]):
            bands_file.write("{:>12s}".format("E_{" + str(band + 1) + "}(k)"))
        bands_file.write("\n")
        for kpoint in range(kpts.shape[0]):
            bands_file.write("{:>7.4f}".format(kpoint_length[kpoint]))
            for band in range(e.shape[0]):
                bands_file.write("{:>12.4e}".format(e[band, kpoint]))
            bands_file.write("\n")
        file_handler(filename, bands_file)

    if datatype == "v":
        vel, kpts = bs.fetch_velocities_along_line(
            ks, ke, itype=itype, itype_sub=itype_sub)
        kpoint_length = np.linalg.norm(kpts - kpts[0], axis=1)
        bands_file = file_handler("output/" + filename, status='w')
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "#                 Bandstructure lines                    #\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "#  1. Column: k-point distance in AA^-1                  #\n")
        bands_file.write(
            "#  2. Column: electron group velocity dispersion for the #\n")
        bands_file.write(
            "#     first band in eVAA                                 #\n")
        bands_file.write(
            "#  3. Column: electron group velocity dispersion for the #\n")
        bands_file.write(
            "#     second band in eVAA                                #\n")
        bands_file.write(
            "#  n. Column: electron group velocity dispersion for the #\n")
        bands_file.write(
            "#     (n-1) band in eVAA                                 #\n")
        bands_file.write(
            "#                                                        #\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "# ------------------------------------------------------ #\n")
        bands_file.write(
            "##########################################################\n")
        bands_file.write("{:<7s}".format("Dist."))
        for band in range(vel.shape[0]):
            for i in range(3):
                bands_file.write("{:>16s}".format("dE_{" + str(band + 1) +
                                                  "}(k)/dk_" + str(i + 1)))
        bands_file.write("\n")
        for kpoint in range(kpts.shape[0]):
            bands_file.write("{:>7.4f}".format(kpoint_length[kpoint]))
            for band in range(vel.shape[0]):
                bands_file.write(
                    " " + "{:>16.4e}".format(vel[band, 0, kpoint]) + " " +
                    "{:>16.4e}".format(vel[band, 1, kpoint]) + " " +
                    "{:>16.4e}".format(vel[band, 2, kpoint]))
            bands_file.write("\n")
        file_handler(filename, bands_file)


def file_handler(filename="", handler=None, status=None):
    """
    Open and close files

    Parameters
    ----------
    filename : string
        Filename to be handled
    handler : object, optional
        A file object. If provided, this routine closes the file
    status : {"w", "r", "a"}
        The status, e.g. write, read, append etc.

    Returns
    -------
    file_handler : object
        A file object

    """
    # get logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access

    if status is None:
        if file_handler is None:
            logger.error("Could not close an empty file handler. Exiting.")
            sys.exit(1)
        handler.close()
        return None

    try:
        handler = open(filename, status)
        return handler
    except IOError:
        logger.error("Could not open %s. Exiting.", filename)
        sys.exit(1)


def start_message():
    """
    Prints a startup message to the log file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    # get logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access

    # dump logo
    logger.info(
        "Starting T4ME. %s\n"
        "T4ME - Transport for Materials\n"
        "Version: %s\n"
        "License: GNU GPL v3\n"
        "Documentation: https://espenfl.github.io/t4me \n"
        "Git repo: git@github.com:espenfl/t4me.git\n"
        "Developed by: Espen Flage-Larsen\n"
        "Contact: espen.flage-larsen@sintef.no\n"
        "Additional contributions or extensions of the code\n"
        "are welcome and greatly appreciated. Please contact\n"
        "the developer to coordinate.\n"
        "\n\n"
        "////   Starting calculations according to set configuration.   ////\n",
        constants.logo, constants.version)


def end_message():
    """
    Prints an end message to the log file.

    Parameters
    ----------
    None

    Returns
    -------
    None

    """

    # get logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access

    # dump logo
    logger.info(
        "Calculations finished.\n"
        "\n\n////   End of calculations                                     ////\n"
    )
