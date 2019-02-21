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
"""Contains various constants used by T4ME."""

# pylint: disable=invalid-name

import math
import numpy as np

version = "1.0"

kb = 8.6173324  # eV/K *10^-5 is not defined here
hbar = 6.58211928  # eV s *10^-16 is not defined here
jtoev = 6.24150934  # eV *10^18 is not defined here
evtoj = 1.602176565  # 10^-19 is not defined here
g0 = 7.7480917346  # S, *10^-5 is not defined here
knorm = 2.0836618  # Hz/K, *10^10 is not defined here
enorm = 2.417989348  # A/J *10^14 is not defined here
magnorm = 13.99624555  # Hz/T *10^9 is not defined here
elcharge = 1.602176565  # C *10^-19 is not defined here
hbarsqcsq = 197.3269718  # (MeV fm)^2 is not defined here
elmasscsq = 0.510998928  # MeV is not defined here
elmass = 9.10938291  # kg *10^-31 is not defined here
umass = 1.660538921  # kg *10^-27 is not defined here
vacperm = 8.854187817  # F/m *10^-12 is not defined here
sqmcsq = math.sqrt(elmasscsq) / hbarsqcsq  # sqrt(eV)/MeV fm
# cm^3 * 10^-22.5 not defined here
sqmcsqinvthird = math.pow(pow(hbarsqcsq, 2.0) / (elmasscsq * kb), 1.5)
# 10^22.5/cm^3 not defined here
scmcsqthird = math.pow(elmasscsq * kb / pow(hbarsqcsq, 2.0), 1.5)
# (MeV (fm^2))^-1 not defined here
mcsqoverhbarcsq = elmasscsq / math.pow(hbarsqcsq, 2.0)
sqkb = math.sqrt(kb * 1e-5)  # sqrt(eV/K)
finestruct = 7.2973525664  # unitless
pi = math.pi
pisq = math.pow(math.pi, 2.0)  # unitless
piqu = math.pow(math.pi, 3.0)  # unitless
bandunit = 1e-4 * math.pow(hbarsqcsq, 2.0) / elmasscsq / 2
zerocut = 1e-8
zeroshift = 1e-10
large = np.finfo(dtype=float).max
largem = -np.finfo(dtype=float).max
small = np.finfo(dtype=float).min
zero = np.finfo(dtype=float).tiny  # pylint: disable=no-member
inf = np.inf

rbf_functions = [
    "multiquadric", "inverse_multiquadric", "gaussian", "linear", "cubic",
    "quintic", "thin_plate"
]

einspline_boundary_cond = ["periodic", "flat", "natural", "antiperiodic"]

wildmagic_methods = [
    "trilinear", "tricubic_exact", "tricubic_bspline", "akima"
]

interp1d_methods = [
    "linear", "nearest", "zero", "slinear", "quadratic", "cubic"
]

elements = {
    "x": 0,
    "h": 1,
    "he": 2,
    "li": 3,
    "be": 4,
    "b": 5,
    "c": 6,
    "n": 7,
    "o": 8,
    "f": 9,
    "ne": 10,
    "na": 11,
    "mg": 12,
    "al": 13,
    "si": 14,
    "p": 15,
    "s": 16,
    "cl": 17,
    "ar": 18,
    "k": 19,
    "ca": 20,
    "sc": 21,
    "ti": 22,
    "v": 23,
    "cr": 24,
    "mn": 25,
    "fe": 26,
    "co": 27,
    "ni": 28,
    "cu": 29,
    "zn": 30,
    "ga": 31,
    "ge": 32,
    "as": 33,
    "se": 34,
    "br": 35,
    "kr": 36,
    "rb": 37,
    "sr": 38,
    "y": 39,
    "zr": 40,
    "nb": 41,
    "mo": 42,
    "tc": 43,
    "ru": 44,
    "rh": 45,
    "pd": 46,
    "ag": 47,
    "cd": 48,
    "in": 49,
    "sn": 50,
    "sb": 51,
    "te": 52,
    "i": 53,
    "xe": 54,
    "cs": 55,
    "ba": 56,
    "la": 57,
    "ce": 58,
    "pr": 59,
    "nd": 60,
    "pm": 61,
    "sm": 62,
    "eu": 63,
    "gd": 64,
    "tb": 65,
    "dy": 66,
    "ho": 67,
    "er": 68,
    "tm": 69,
    "yb": 70,
    "lu": 71,
    "hf": 72,
    "ta": 73,
    "w": 74,
    "re": 75,
    "os": 76,
    "ir": 77,
    "pt": 78,
    "au": 79,
    "hg": 80,
    "tl": 81,
    "pb": 82,
    "bi": 83,
    "po": 84,
    "at": 85,
    "rn": 86,
    "fr": 87,
    "ra": 88,
    "ac": 89,
    "th": 90,
    "pa": 91,
    "u": 92,
    "np": 93,
    "pu": 94,
    "am": 95,
    "cm": 96,
    "bk": 97,
    "cf": 98,
    "es": 99,
    "fm": 100,
    "md": 101,
    "no": 102,
    "lr": 103,
    "rf": 104,
    "db": 105,
    "sg": 106,
    "bh": 107,
    "hs": 108,
    "mt": 109,
    "ds": 110,
    "rg": 111,
    "cn": 112,
    "uut": 113,
    "uuq": 114,
    "uup": 115,
    "uuh": 116,
    "uus": 117,
    "uuo": 118,
}

logo = r"""                                            \n
        _________________  ____  ___    _______________\n
       /            /   / /   / /   \__/   /          /\n
      /____    ____/   / /   / /          /   _______/ \n
          /   /   /   /_/   /_/          /   /___      \n
         /   /   /           /   /\_/   /   ____/      \n
        /   /   /_____    __/   /  /   /   /_______    \n
       /   /         /   / /   /  /   /           /    \n
      /___/         /___/ /___/  /___/___________/     \n\n"""
