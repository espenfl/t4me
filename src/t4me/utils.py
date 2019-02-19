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
"""Containing utilitary functions for T4ME."""

#!/usr/local/bin/python
import sys
import logging
import logging.config
import os
import yaml
import numpy as np

import t4me.constants as constants  #pylint: disable=useless-import-alias


def check_file(filename):
    """
    Check if a file exists

    Parameters
    ----------
    filename : string
        The location and filename of the file to be checked.

    Returns
    -------
    None

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  #pylint: disable=protected-access
    logger.debug("Running check_file.")

    if not os.path.isfile(filename):
        logger.error("%s was not found. Exiting.", filename)
        sys.exit(1)


def create_directory(path):
    """
    Check that the directory exists

    Parameters
    ----------
    path : string
        The path to the directory to be checked.

    Returns
    -------
    None

    """

    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def check_directory(path, create=False):
    """
    Check that a directory exists.

    Parameters
    ----------
    path : string
        The path to the directory to be checked.
    create : boolean, optional
        If set to True create the directory if it
        does not exist. Defaults to False and in this
        case an error is printed if the directory is
        not found.

    Returns
    -------
    None

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  #pylint: disable=protected-access
    logger.debug("Running check_file.")

    if not os.path.isdir(path):
        if not create:
            # THIS SHOULD AT SOME POINT BE REPLACED WITH
            # THE PROPER LOCAL LOGGER. AT LEAST A CHECK IF
            # IT EXISTS.
            logger.error("The directory %s does not exist. Exiting.", path)
            sys.exit(1)
        else:
            create_directory(path)


def clean_directory(path):
    """
    Clean a directory.

    Parameters
    ----------
    path : string
        The path to the directory to be cleaned.

    Returns
    -------
    None

    """
    check_directory(path, create=True)
    os.system("rm -rf " + path + "/*")


def is_number(something):
    """
    Check if `something` is a number.

    Parameters
    ----------
    something : anything
        Something to be checked.

    Returns
    -------
    boolean
        True if `something` is a number. False otherwise.

    """
    try:
        float(something)
        return True
    except ValueError:
        return False


def is_even(number):
    """
    Check if `number` is even.

    Parameters
    ----------
    number : integer
        The integer to be checked

    Returns
    -------
    boolean
        Returns True of `number` is even, False otherwise.

    """
    return number % 2 == 0


def is_power_of_two(number):
    """
    Check that if a number is a power of two.

    Parameters
    ----------
    number : float
        The supplied number to be checked.

    Returns
    -------
    boolean
        Returns True of `number` is power of two, False
        otherwise.
    """

    test_number = 1
    while test_number < number:
        test_number <<= 1

    return test_number == number


def pull_vecs_inside_boundary(vecs, border, shift=None):
    """
    Pulls vectors into a given cubic boundary box.

    Parameters
    ----------
    vecs : ndarray
        | Dimension: (N,3)

        Contains N vectors.
    border : ndarray
        | Dimension: (6)

        Contains the entries x_min, x_max, y_min, y_max, z_min and z_max,
        respectively of the indexes to be modified, typically the border
        elements.
    shift : float, optional
        An optional shift value which brings the vectors shift more inside
        the boundary box supplied in `border`.

    Returns
    -------
    None

    """

    if shift is None:
        shift = 0.0

    for vec in vecs:
        if vec[0] < border[0]:
            vec[0] = border[0] + shift
        if vec[1] < border[2]:
            vec[1] = border[2] + shift
        if vec[2] < border[4]:
            vec[2] = border[4] + shift
        if vec[0] > border[1]:
            vec[0] = border[1] - shift
        if vec[1] > border[3]:
            vec[1] = border[3] - shift
        if vec[2] > border[5]:
            vec[2] = border[5] - shift


def invert_matrix(matrix):
    """
    Inverts a matrix and checks for ill-conditions

    Parameters
    ----------
    matrix : ndarray
        | Dimension: (N,N)

        The input matrix to be inverted.

    Returns
    -------
    inv_matrix : ndarray
        | Dimension: (N,N)

        The inverted matrix. If `matrix` is ill-conditioned, nan
        values are filled in the matrix.

    """

    np.seterr(divide="ignore")
    if np.linalg.cond(matrix) < (1.0 / np.finfo(matrix.dtype).eps):
        inv_matrix = np.linalg.inv(matrix)
        return inv_matrix
    # matrix is ill conditioned, e.g. singular and linalg.inv
    # fails, return matrix with nan instead
    inv_matrix = np.full_like(matrix, np.nan)

    return inv_matrix


def config_logger(filename="/logging.yaml", level=None):
    """
    Configure the main logger.

    Parameters
    ----------
    filename : string, optional
        The filename for the logging configuration
        file. Defaults to "logging.yaml" in the current
        working directory.
    level : object
        Sets the logging level of the Python logger.
        Defaults to `INFO` if the configuration file is
        not found.

    Returns
    -------
    None

    """

    path_to_script = os.path.dirname(os.path.realpath(__file__))
    filename = path_to_script + filename
    if os.path.exists(filename):
        logging.info("Setting logger from logging.yaml")
        with open(filename, 'rt') as file_data:
            config = yaml.load(file_data.read())
        logging.config.dictConfig(config)
    else:
        if level is None:
            level = logging.INFO
        logging.info("Setting default logger")
        logging.basicConfig(level=level)


def fetch_sorting_indexes(data, order="C"):
    """
    Fetch the sorting indexes to sort an array to either column or row order.

    Parameters
    ----------
    data : ndarray
        | Dimension: (N,M)

        The input data array to be sorted.

    order : {"C", "F"}

        The sort order.

    Returns
    -------
    sort_index : ndarray
        | Dimension: (N,M)

        The sorting indexes that can be used to order the array.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  #pylint: disable=protected-access
    logger.debug("Running fetch_sorting_indexes.")

    # before sort, round the array to predefined accuracy in
    # order to avoid scattered values, inform the user
    logger.debug("Rounding the array before sort.")
    decimals = int(str(constants.zerocut).split("-")[1].strip("0"))
    data = np.around(data, decimals=decimals)

    if order == "C":
        return np.lexsort(np.transpose(data)[::-1])
    if order == "F":
        return np.lexsort(np.transpose(data))
    logger.error("The requested sort order is not supported. "
                 "Only C or F (Fortran) order is supported. "
                 "Exiting.")
    sys.exit(1)


def pull_points_back_into_zone(points):
    """
    Pulls all points outside [-0.5,0.5] in direct coordinates back into [-0.5, 0.5].

    Parameters
    ----------
    points : ndarray
        | Dimension: (N, 3)

        The N points to be checked and thrown back into
        the zone (between [-0.5, 0.5]). Should be in direct coordinates.

    Returns
    -------
    None

    """

    # do positive side
    points[np.where(points > 0.5)] = \
        points[np.where(points > 0.5)] - 1.0

    # then negative side
    points[np.where(points < -0.5)] = \
        points[np.where(points < -0.5)] + 1.0
