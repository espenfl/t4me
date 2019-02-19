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

#!/usr/bin/python
"""Contains routines that interface T4ME, e.g. to the parameter files or to other input files."""

# pylint: disable=useless-import-alias, too-many-arguments, invalid-name, too-many-statements, too-many-lines, global-statement

import sys
import logging
import xml.etree.cElementTree as ET
import numpy as np

import t4me.utils as utils
import t4me.inputoutput as inputoutput
import t4me.constants as constants


def lattice_param_numpy(lattice, location=None, filename=None):
    r"""
    Interface used to format the elements needed to generate the `Lattice()` object.

    Used if the lattice is generated from the celldata YAML file (parameterfile and Numpy intput files).

    Parameters
    ----------
    lattice : object
        A `Lattice()` object where we can store additional parameters
        detected during setup for later access.
    location : string, optional
        The location of the YAML parameter file determining the celldata.
    filename : string, optional
        The filename of the YAML parameter file determining the celldata.

    Returns
    -------
    unitcell : ndarray
        | Dimension: (3,3)

        The unitcell in cartesian coordinates and {\AA} units.
    positions : ndarray
        | Dimension: (N,3)

        The positions of the N atoms in the unitcell in cartesian
        coordinates.
    species : ndarray
        | Dimension: (N)

        Integer atomic numbers of the atomic species in the same
        order as positions. Hydrogen starts with 1, while the element X
        is located at 0. Otherwise it follows the periodic table.
    kmesh : object
        | Dimension: (3)

        A `Kmesh()` obhect for the reciprocal mesh generation
        containment. Should include `sampling`, `mesh`, `mesh_ired`
        and other parameters needed for later processing.

    Notes
    -----
    Upon writing a custom interface, please make sure that the
    parameters in the YAML files are not overwritten.

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running lattice_param_numpy.")

    data = inputoutput.readcellparam(location, filename)
    unitcell = np.ascontiguousarray(
        np.vstack((np.array(data["a"], dtype='double'),
                   np.array(data["b"], dtype='double'),
                   np.array(data["c"], dtype='double'))),
        dtype='double')
    atomtypes = data["atomtypes"]
    species = np.array(
        [constants.elements[atom.lower()] for atom in atomtypes],
        dtype='intc',
        order='C')
    if len(species.shape) != 1:
        logger.error(
            "Please enter the atomic numbers as [n1 n2] etc. Exiting.")
        sys.exit(1)
    positions = np.array(data["pos"], dtype='double', order='C')

    lattice.unitcell = unitcell
    lattice.species = species

    if len(positions.shape) < 2:
        logger.error("Please enter the atomic positions as "
                     "[[x1 y1 z1],[x2,y2,z2]] etc. Exiting.")
        sys.exit(1)
    direct_positions = data["direct"]
    # if positions is given in direct coordinates, convert to cartesian
    # consider to change this in the future
    if direct_positions:
        positions = lattice.dir_to_cart(positions, real=True)

    lattice.positions = positions

    lattice.kdata.sampling = np.array(
        data["ksampling"], dtype='intc', order='C')

    # include borders
    lattice.kdata.borderless = False

    # no mesh or mapping relations read, generate later
    lattice.kdata.k_sort_index = None
    lattice.kdata.mapping_ibz_to_bz = None
    lattice.kdata.mapping_bz_to_ibz = None
    lattice.kdata.ibz_weights = None
    lattice.kdata.mesh = None
    lattice.kdata.mesh_ired = None

    # for the analytic models we need to generate the
    # grid and mapping tables
    # for NumPy data we assume that the user inputs
    # the full grid
    if not lattice.param.read[:5] == "numpy":
        lattice.param.work_on_full_grid = False
    else:
        lattice.param.work_on_full_grid = True


def lattice_vasp(lattice, location=None, filename=None):  # pylint: disable=too-many-locals # noqa: MC0001
    r"""
    Interface used to format the elements needed to generate the `Lattice()` object.

    Used if the lattice is generated from the VASP XML file.

    Parameters
    ----------
    lattice : object
        A `Lattice()` object where we can store additional parameters
        detected during setup for later access.
    location : string, optional
        The location of the VASP XML file determining the celldata.
    filename : string, optional
        The filename of the VASP XML file determining the celldata.

    Returns
    -------
    unitcell : ndarray
        | Dimension: (3,3)

        The unitcell in cartesian coordinates and {\AA} units.
    positions : ndarray
        | Dimension: (N,3)

        The positions of the N atoms in the unitcell in cartesian
        coordinates.
    species : ndarray
        | Dimension: (N)

        Integer atomic numbers of the atomic species in the same
        order as positions. Hydrogen starts with 1, while the
        element X is located at 0. Otherwise
        it follows the periodic table.
    kmesh : object
        | Dimension: (3)

        A `Kmesh()` object for the reciprocal mesh generation
        containment. Should include `sampling`, `mesh`, `mesh_ired`
        and other parameters needed for later processing.

    Notes
    -----
    Upon writing a custom interface, please make sure that the
    parameters in the YAML parameter files are not overwritten.

    Additional parameters pertaining VASP are stored inside the
    `Param()` object with a `vasp` preemble, i.e.
    `param.vasp_something`.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running lattice_vasp.")

    if filename is None:
        # check param as well
        vasp_filename = lattice.param.readfile
        if not vasp_filename:
            vasp_filename = "vasprun.xml"
    if location is not None:
        vasp_filename = location + "/" + vasp_filename
    else:
        vasp_filename = "input/" + vasp_filename
    utils.check_file(vasp_filename)
    tree = ET.ElementTree(file=vasp_filename)
    # fetch symprec
    lattice.param.vasp_symprec = float(
        tree.find('.//parameters/separator[@name="symmetry"]/'
                  'i[@name="SYMPREC"]').text)
    # fetch cell
    unitcell_vectors = tree.findall('.//structure[@name="finalpos"]/crystal/'
                                    'varray[@name="basis"]/v')
    unitcell = np.zeros((3, 3))
    for index, unitcell_vector in enumerate(unitcell_vectors):
        unitcell[index] = np.fromstring(unitcell_vector.text, sep=' ')
    # fetch positions
    positions_entries = tree.findall('.//structure[@name="finalpos"]/'
                                     'varray[@name="positions"]/v')
    positions = np.zeros((len(positions_entries), 3),
                         dtype='double',
                         order='C')
    species_entries = tree.findall('.//atominfo/'
                                   'array[@name="atoms"]/set/rc')
    species = np.zeros(len(positions_entries), dtype='intc', order='C')
    for index, position in enumerate(positions_entries):
        positions[index] = np.fromstring(position.text, sep=' ')
        species[index] = constants.elements[species_entries[index]
                                            [0].text.split()[0].lower()]

    lattice.unitcell = unitcell
    lattice.positions = positions
    lattice.species = species

    # fetch KINTER
    kinter = 0
    kinter_param = tree.find('.//incar/i[@name="KINTER"]')
    if kinter_param is not None:
        kinter = int(kinter_param.text)
    lattice.param.vasp_kinter = kinter

    # fetch LVEL
    lvel = False
    lvel_param = tree.find('.//incar/i[@name="LVEL"]')
    if lvel_param is not None:
        if "T" in lvel_param.text:
            lvel = True
    # now check for if velocities are obtained by Wannier
    if not lvel:
        lvel_param = tree.find('.//incar/i[@name="LINTPOL_VELOCITY"]')
        if lvel_param is not None:
            if "T" in lvel_param.text:
                lvel = True
    lattice.param.vasp_lvel = lvel

    # fetch LWANNIERINTERPOL
    lwan = False
    lwan_param = tree.find('.//incar/i[@name="LWANNIERINTERPOL"]')
    if lwan_param is not None:
        if "T" in lwan_param.text:
            lwan = True
    lattice.param.vasp_lwan = lwan

    mapping_bz_to_ibz_vasp = None

    # get divisions, IBZ kpoints and location of the eigenvalues and
    # dos
    divisions = tree.find('kpoints/generation/v[@name="divisions"]')
    lattice.kdata.sampling = np.ascontiguousarray(
        np.fromstring(divisions.text, sep=" ", dtype='intc'), dtype='intc')
    if kinter:
        lattice.kdata.sampling = lattice.kdata.sampling * abs(kinter)
    else:
        lattice.kdata.sampling = np.ascontiguousarray(
            np.fromstring(divisions.text, sep=" ", dtype='intc'), dtype='intc')
    if not lvel:
        kpoints = tree.findall('kpoints/varray[@name="kpointlist"]/v')
    else:
        kpoints = tree.findall('.//kpoints[@comment="interpolated"]/'
                               'varray[@name="kpointlist"]/v')
        kpoints_full = tree.findall(
            './/eigenvelocities[@comment="interpolated"]/kpoints/'
            'varray[@name="kpointlist"]/v')
        mapping_bz_to_ibz = tree.findall(
            './/eigenvelocities[@comment="interpolated"]/kpoints/'
            'varray[@name="ibzequiv"]/v')

    # fetch IBZ points
    kpointsvasp = np.zeros((len(kpoints), 3), dtype='double', order='C')
    for index, kpoint in enumerate(kpoints):
        kpointsvasp[index] = np.fromstring(kpoint.text, sep=' ')

    # pull k-points back into zone
    utils.pull_points_back_into_zone(kpointsvasp)
    # now sort according to k-point sort (z increasing
    # fastests) and store
    k_sort_index = utils.fetch_sorting_indexes(kpointsvasp)
    lattice.kdata.mesh_ired = np.ascontiguousarray(
        kpointsvasp[k_sort_index], dtype="double")
    lattice.kdata.mesh = None
    # then also full BZ (for e.g. velocities)
    if lvel:
        kpointsvasp_full = np.zeros((len(kpoints_full), 3),
                                    dtype='double',
                                    order='C')
        mapping_bz_to_ibz_vasp = np.zeros(len(kpoints_full), dtype='intc')
        for index, kpoint in enumerate(kpoints_full):
            kpointsvasp_full[index] = np.fromstring(kpoint.text, sep=' ')
        for index, mapping in enumerate(mapping_bz_to_ibz):
            mapping_bz_to_ibz_vasp[index] = np.fromstring(
                mapping.text, sep=' ')
        # pull k-points back into zone
        utils.pull_points_back_into_zone(kpointsvasp_full)
        k_sort_index = utils.fetch_sorting_indexes(kpointsvasp_full)
        lattice.kdata.mesh = np.ascontiguousarray(
            kpointsvasp_full[k_sort_index], dtype='double')

    # need this later to sort eigenvalues etc. (only store bz
    # if that is present, otherwise ibz)
    lattice.kdata.k_sort_index = k_sort_index

    # VASP grids are borderless
    lattice.kdata.borderless = True

    # check if bz to ibz mapping is read, if so generate ibz to bz
    # this is only done if the full grid data is present
    if mapping_bz_to_ibz_vasp is not None:
        shuffle = np.zeros(mapping_bz_to_ibz_vasp.shape[0], dtype=np.intc)
        for index, value in enumerate(mapping_bz_to_ibz_vasp):
            shuffle[index] = np.where(k_sort_index == value)[0][0]

        mapping_ibz_to_bz_vasp = shuffle[k_sort_index]
        for index, ibz_point in enumerate(np.unique(mapping_ibz_to_bz_vasp)):
            mask = np.in1d(mapping_ibz_to_bz_vasp, ibz_point)
            np.copyto(mapping_bz_to_ibz_vasp, index, where=mask)

        lattice.kdata.mapping_bz_to_ibz = mapping_bz_to_ibz_vasp
        lattice.kdata.mapping_ibz_to_bz = mapping_ibz_to_bz_vasp
    else:
        lattice.kdata.mapping_bz_to_ibz = None
        lattice.kdata.mapping_ibz_to_bz = None

    # FIX THIS AS WE NOW SIMPLY CAN USE THE MAPPING TO SET THE WEIGHTS
    lattice.kdata.ibz_weights = None

    # if LVEL have been set we already have access to the full grid
    # so set such a parameter
    lattice.param.work_on_full_grid = bool(lvel)


def lattice_w90(lattice):  # pylint: disable=too-many-locals, too-many-branches
    r"""
    Interface used to format the elements needed to generate the `Lattice()` object.

    Used if the lattice is generated from the Wannier90 win file.

    Parameters
    ----------
    lattice : object
        A `Lattice()` object where we can store additional
        parameters detected during setup for later access.

    Returns
    -------
    unitcell : ndarray
        | Dimension: (3,3)

        The unitcell in cartesian coordinates and {\AA} units.
    positions : ndarray
        | Dimension: (N,3)

        The positions of the N atoms in the unitcell in cartesian
        coordinates.
    species : ndarray
        | Dimension: (N)

        Integer atomic numbers of the atomic species in the same
        order as positions. Hydrogen starts with 1, while the element
        X is located at 0. Otherwise it follows the periodic table.
    kmesh : object
        A `Kmesh()` object for the reciprocal mesh generation containment.
        Should include `sampling`, `mesh`, `mesh_ired` and other
        parameters needed for later processing.

    Notes
    -----
    Upon writing a custom interface, please make sure that the parameters
    in the YAML parameter file is not overwritten.

    Additional parameters pertaining VASP are stored inside the `Param()`
    object with a `vasp` preemble, i.e. `param.vasp_something`.

    """
    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running lattice_w90.")

    # check filename
    if lattice.param.readfile == "":
        logger.info(
            "User did not specificy the prefix to the Wannier90 output "
            "files. Setting it to 'wannier90' and continuing.")
        wprefix = "wannier90"
    else:
        wprefix = lattice.param.readfile
    wfile = open("input/" + wprefix + ".win", "r")
    wdata = wfile.readlines()
    positions = []
    atomtypes = []
    positions_start = 0
    positions_end = 0
    positions_cart = False
    for i, line in enumerate(wdata):
        if "begin unit_cell_cart" in line.lower():
            unitcell = np.array([[float(entry) for entry in element.split()]
                                 for element in wdata[i + 1:i + 4]],
                                dtype='double',
                                order='C')
        if "mp_grid" in line.lower():
            sampling = np.array(
                [int(element) for element in wdata[i].split()[2:6]],
                dtype='intc',
                order='C')
        if "begin atoms_cart" in line.lower():
            positions_cart = True
            positions_start = i + 1
        if "begin atoms_frac" in line.lower():
            positions_start = i + 1
        if "end atoms_cart" in line.lower():
            positions_end = i
        if "end atoms_frac" in line.lower():
            positions_end = i
        if "begin kpoints" in line.lower():
            kpoints_start = i + 1
        if "end kpoints" in line.lower():
            kpoints_end = i

    numkpoints = np.prod(sampling)
    kpointswannier = np.zeros((numkpoints, 3))
    for i, line in enumerate(wdata[positions_start:positions_end]):
        splitted = line.split()
        positions.append([float(element) for element in splitted[1:4]])
        atomtypes.append(splitted[0])
    positions = np.array(positions, dtype='double', order='C')
    species = np.array(
        [constants.elements[atom.lower()] for atom in atomtypes],
        dtype='intc',
        order='C')

    for i, line in enumerate(wdata[kpoints_start:kpoints_end]):
        kpointswannier[i] = np.array(
            [float(element) for element in line.split()])

    # pull k-points back into zone
    utils.pull_points_back_into_zone(kpointswannier)

    # fetch sorting indexes
    k_sort_index = utils.fetch_sorting_indexes(kpointswannier)

    # sort k-points
    kpointswannier = np.ascontiguousarray(
        kpointswannier[k_sort_index], dtype='double')

    lattice.unitcell = unitcell
    # convert positions to direct coordinates
    if positions_cart:
        positions = lattice.cart_to_dir(positions, real=True)
    lattice.positions = positions
    lattice.species = species
    lattice.kdata.sampling = sampling
    lattice.kdata.mesh = kpointswannier
    lattice.kdata.k_sort_index = k_sort_index
    # Assume borderless, need check for this
    lattice.kdata.borderless = True
    lattice.kdata.mesh_ired = None
    lattice.kdata.mapping_bz_to_ibz = None
    lattice.kdata.mapping_ibz_to_bz = None
    lattice.kdata.ibz_weights = None


def bandstructure_param(bs, location=None, filename=None):
    """
    Sets the bandstructure from the parameters in the bandstructure configuration file (default bandparam.yml).

    Also loads and stores the parameters.

    Parameters
    ----------
    bs : object
        A `Bandstructure()` object.
    location : string, optional
        The location of the bandstructure configuration
        file. Defaults later to the "input" directory
        in the current working directory.
    filename : string, optional
        The filename of the bandstructure configuration
        file. Defaults to "bandparam.yml".

    Returns
    -------
    None

    Notes
    -----
    This interface prepares analytic and tight binding
    generation of the band structure and also loads and
    stores all bandstructure related parameters.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running bandstructure_param.")

    # read the band parameter file
    read_band_parameters(bs, 0, location, filename)

    # finally we generate the raw bandstructure
    bs.energies, bs.velocities, bs.tb_band = bs.gen_bands()
    # what about energy shifts? if the user set any of the
    # parameters relevant for the Fermi level, give warning
    # and continue
    if (bs.param.e_fermi or bs.param.e_fermi_in_gap or bs.param.e_vbm):
        logger.info("User have set 'e_fermi' and/or 'e_fermi_in_gap' "
                    "and/or 'e_vbm' to True, but this options are "
                    "not supported for the parametric band "
                    "construction. Shifting the bandstructure "
                    "by 'e_shift' and continuing.")
        # velocities are of course independent of this shift
        bs.energies = bs.energies - bs.param.e_shift

    # now calculate the density of states
    bs.dos_partial, bs.dos_energies = bs.gen_dos()
    if bs.dos_partial is None:
        bs.dos = None
    else:
        bs.dos = bs.dos_partial.sum(-2)

    # here we have the velocities
    bs.gen_velocities = False

    # these can be set manually, fix later
    bs.vbm_energy = None
    bs.vbm_band = None
    bs.vbm_kpoint = None
    bs.cbm_energy = None
    bs.cbm_band = None
    bs.cbm_kpoint = None
    bs.band_gap = None
    bs.direct = None


def bandstructure_vasp(bs, location=None, filename=None):  # pylint: disable=too-many-locals # noqa: MC0001
    """
    Sets the bandstructure from a VASP XML file.

    Loads and stores the parameters in the bandstructure configuration file (defaults to bandparam.yml).

    Parameters
    ----------
    bs : object
        A `Bandstructure()` object.
    location : string, optional
        The location of the VASP XML file. Defaults
        to the "input" directory in the current working directory.
    filename : string, optional
        The filename of the VASP XML file to be read.
        Defaults to "vasprun.xml". The bandstructure
        configuration file have to be named "bandparam.yml"
        in this case.

    Returns
    -------
    None

    Notes
    -----
    This interface read and sets up the bandstructure based
    on a VASP XML file. Currently it does not read the band
    velocities as VASP does not yet support this feature.
    However, work is in progress to enable this. The band
    velocities have to be generated by an interpolation routine
    later. Flags are automatically set for this. The bandstructure
    configuration file is still read due to the need of the
    scattering properties etc.

    This interface is enabled by setting `read` in the general
    configuration file to vasp.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running bandstructure_vasp.")

    if filename is None:
        # check param as well
        vasp_filename = bs.param.readfile
        if not vasp_filename:
            vasp_filename = "vasprun.xml"
    if location is not None:
        vasp_filename = location + "/" + vasp_filename
    else:
        vasp_filename = "input/" + vasp_filename
    utils.check_file(vasp_filename)
    tree = ET.ElementTree(file=vasp_filename)

    # fetch ispin
    ispin = int(
        tree.find('.//parameters/separator[@name="electronic"]/'
                  'separator[@name="electronic spin"]/'
                  'i[@name="ISPIN"]').text)

    # if ISPIN = 2 quit (not sufficiently tested)
    if ispin == 2:
        logger.error("The code is not yet sufficiently tested for "
                     "different spin channels (ISPIN=2). Exiting. ")

    # fetch nedos
    num_samples_dos = int(
        tree.find('.//parameters/separator[@name="dos"]/'
                  'i[@name="NEDOS"]').text)

    # check first if list generated k-points, then break
    if tree.find('.//kpoints/generation[@param="listgenerated"]'):
        logger.error("The supplied vasprun.xml contain list generated k-point "
                     "sets. Please rerun VASP using a dense converged IBZ or "
                     "with KINTER/LVEL. Exiting.")
        sys.exit(1)

    # fetch KINTER
    kinter = bs.param.vasp_kinter

    # fetch LVEL
    lvel = bs.param.vasp_lvel

    # fetch LWANNIERINTERPOL
    lwan = bs.param.vasp_lwan

    # get divisions, IBZ kpoints and location of the eigenvalues and
    # dos
    if kinter > 2:
        if lvel:
            energies_base = tree.find(
                './/eigenvelocities[@comment="interpolated"]/eigenvalues')
        else:
            energies_base = tree.find(
                './/eigenvalues[@comment="interpolated"]')
        dos_base = tree.find('.//dos[@comment="interpolated"]')
    else:
        energies_base = tree.find('.//eigenvalues')
        dos_base = tree.find('.//dos')

    # initialize arrays
    # if ISPIN=2, pad energies (num bands=2*NBANDS, set of down, then
    # set of up)
    numkpoints = bs.lattice.kmesh.shape[0]
    numbands = int(
        tree.find('.//parameters/separator[@name="electronic"]/'
                  'i[@name="NBANDS"]').text)
    if ispin > 1:
        energiesvasp = np.zeros((2 * numbands, numkpoints))
        if lvel:
            velocitiesvasp = np.zeros((2 * numbands, 3, numkpoints))
        else:
            velocitiesvasp = None
        occvasp = np.zeros((2 * numbands, numkpoints))
        dosvasp = np.zeros((2, num_samples_dos, 2))
        bs.spin_degen = np.full(2 * numbands, 1, dtype="intc")
    else:
        energiesvasp = np.zeros((numbands, numkpoints))
        if lvel:
            velocitiesvasp = np.zeros((numbands, 3, numkpoints))
        else:
            velocitiesvasp = None
        occvasp = np.zeros((numbands, numkpoints))
        dosvasp = np.zeros((num_samples_dos, 2))
        bs.spin_degen = np.full(numbands, 2, dtype="intc")
    # now read energies and/or velocities
    for spin in range(1, ispin + 1):
        energies = energies_base.findall('array/set/set[@comment="spin ' +
                                         str(spin) + '"]/set')
        dos = dos_base.findall('total/array/set/set[@comment="spin ' +
                               str(spin) + '"]/r')
        # read energies and/or velocities
        for idxk, kpoint in enumerate(energies):
            for idxe, data_per_k in enumerate(kpoint):
                data = np.fromstring(data_per_k.text, sep=' ')
                # set index to account for padded values if ISPIN=2
                energy_index = (spin - 1) * numbands + idxe
                energiesvasp[energy_index][idxk] = data[0]
                if (not lvel and kinter < 2):
                    occvasp[energy_index][idxk] = data[1]
                elif lvel:
                    # OCCUPANCIES NOT YET WRITTEN TOGETHER WITH THE
                    # VELOCITIES...PUT ALL TO ONE
                    occvasp[energy_index][idxk] = 1.0
                    velocitiesvasp[energy_index][0][idxk] = data[1]
                    velocitiesvasp[energy_index][1][idxk] = data[2]
                    velocitiesvasp[energy_index][2][idxk] = data[3]
        # read density of states
        for idx, dos_entry in enumerate(dos):
            energy_and_dos = np.fromstring(dos_entry.text, sep=' ', count=3)
            if ispin > 1:
                dosvasp[spin][idx][0] = energy_and_dos[0]
                dosvasp[spin][idx][1] = energy_and_dos[1]
            else:
                dosvasp[idx][0] = energy_and_dos[0]
                dosvasp[idx][1] = energy_and_dos[1]
    # now sort energies, velocities and occ according to k-point
    # sort (z increasing fastest), also set spin_degen
    k_sort_index = bs.lattice.k_sort_index
    for band in np.ndindex(energiesvasp.shape[0]):
        if lvel:
            # if lvel, we have the full BZ dataset, also for the
            # energies
            energiesvasp[band] = energiesvasp[band][k_sort_index]
            velocitiesvasp[band][0] = velocitiesvasp[band][0][k_sort_index]
            velocitiesvasp[band][1] = velocitiesvasp[band][1][k_sort_index]
            velocitiesvasp[band][2] = velocitiesvasp[band][2][k_sort_index]
            occvasp[band] = occvasp[band][k_sort_index]
        else:
            # else, only set energies, as we know we do not have the
            # velocities, IBZ energies are also rotated to full BZ
            num_ibz_kpoints = k_sort_index.size
            # first sort
            ibz_energies = energiesvasp[band][0:num_ibz_kpoints]
            ibz_energies = ibz_energies[k_sort_index]
            # then lay out full bz
            energiesvasp[band] = ibz_energies[bs.lattice.mapping_bz_to_ibz]
            # then occupancies
            ibz_occ = occvasp[band][0:num_ibz_kpoints]
            ibz_occ = ibz_occ[k_sort_index]
            occvasp[band] = ibz_occ[bs.lattice.mapping_bz_to_ibz]

    # fetch efermi
    if kinter < 2:
        fermi_energy = float(
            tree.find('.//calculation/dos/i[@name="efermi"]').text)
    else:
        if not lwan:
            fermi_energy = float(
                tree.find('.//calculation/dos[@comment="interpolated"]/'
                          'i[@name="efermi"]').text)
        else:
            fermi_energy = float(
                tree.find('.//calculation/dos/i[@name="efermi"]').text)

    # if we received data on the full k-point grid we do not yet
    # have access to the occupancies and it makes no real sense to
    # try to detect the gap
    if not lvel:
        bs.vbm_value, bs.vbm_band, bs.vbm_kpoint = bs.locate_vbm(
            energies=energiesvasp, occ=occvasp)
        bs.cbm_value, bs.cbm_band, bs.cbm_kpoint = bs.locate_cbm(
            energies=energiesvasp, occ=occvasp)
        bs.band_gap, bs.direct = bs.locate_bandgap(
            energies=energiesvasp, occ=occvasp)
    else:
        bs.band_gap = None

    if bs.band_gap == 0:
        if bs.direct:
            logger.info("No band gap was detected and it appears as this is "
                        "a metallic system.")
            bs.metallic = True
        else:
            logger.info("No band gap was detected and it appears as this is "
                        "a semi-metallic system.")
            bs.metallic = True

    # no occupancies, set band gap, vbm/cbm
    if bs.band_gap is None:
        logger.info("Occupancies was not available, setting band gap "
                    "to zero and vbm/cbm to the Fermi level supplied by "
                    "VASP.")
        bs.band_gap = 0.0
        bs.vbm_value = fermi_energy
        bs.cbm_value = fermi_energy
    else:
        fermi_diff = fermi_energy - bs.vbm_value
        if fermi_diff > constants.zerocut:
            logger.info(
                "The fermi_energy from VASP differs by the calculate valence "
                "band maximum by: %s", str(fermi_diff))

    e_adjust = 0.0
    # full band
    # first check that the user have only set one shift parameter
    bs.check_energyshifts()
    if ((not bs.param.transport_drop_valence)
            and (not bs.param.transport_drop_conduction)):
        if ((bs.param.e_fermi) and (not bs.param.e_fermi_in_gap)):
            logger.info("Adjusting energies to efermi in vasprun.xml.")
            e_adjust = fermi_energy
        elif ((not bs.param.e_fermi) and (bs.param.e_fermi_in_gap)):
            if not bs.metallic:
                e_adjust = bs.band_gap / 2.0 + bs.vbm_value
            else:
                logger.error(
                    "System appear to be metallic. User requested "
                    "to place Fermi energy in the middle of the band "
                    "gap, which does not exist. Please set Fermi level "
                    "manually or use default. Exiting.")
                sys.exit(1)
        else:
            logger.info("Shifting the energies by e_shift.")
            e_adjust = bs.param.e_shift
    else:
        if not bs.metallic:
            if bs.param.transport_drop_valence:
                logger.info(
                    "Dropping valence band and setting the Fermi level "
                    "to the conduction band minimum.")
                energiesvasp = np.delete(energiesvasp, bs.vbm_band, axis=0)
                e_adjust = bs.cbm_value
            else:
                logger.info("Dropping conduction band and setting the Fermi "
                            "level to the valence band maximum.")
                energiesvasp = np.delete(energiesvasp, bs.cbm_band, axis=0)
                e_adjust = bs.vbm_value
        else:
            logger.error(
                "User tries to only run calculations on the conduction "
                "or valence bands, but the system appears to be metallic"
                ", so this is not possible. Exiting.")
            sys.exit(1)
    # adjust energy
    energiesvasp = energiesvasp - e_adjust
    bs.vbm_value = bs.vbm_value - e_adjust
    bs.cbm_value = bs.cbm_value - e_adjust
    dosvasp[:, 0] = dosvasp[:, 0] - e_adjust
    # add a small shift out of zero where the energies are truly zero
    energiesvasp[np.abs(energiesvasp) < constants.zerocut] = constants.zerocut

    # bandparams contain scattering properties etc.
    data = inputoutput.readbandparam(location, filename)
    bandparams = np.zeros((numbands, 2), dtype=np.int8)
    effmass = np.zeros((numbands, 3))
    select_scattering = np.zeros((numbands, 12), dtype='intc')
    explicit_prefact = np.zeros((numbands, 11), dtype='intc')
    explicit_prefact_values = np.zeros((numbands, 11))
    q_energy_trans = np.zeros((numbands, 2, 3))
    da = np.zeros(numbands)
    do = np.zeros(numbands)
    speed_sound = np.zeros(numbands)
    no = np.zeros(numbands)
    nvv = np.zeros(numbands)
    ni = np.zeros(numbands)
    omegao = np.zeros(numbands)
    omegavv = np.zeros(numbands)
    etrans = np.zeros(numbands)
    rho = np.zeros(numbands)
    zf = np.zeros(numbands)
    f = np.zeros(numbands)
    z = np.zeros(numbands)
    isl = np.zeros(numbands)
    isli = np.zeros(numbands)
    vdiff = np.zeros(numbands)
    alloyconc = np.zeros(numbands)
    p = np.zeros(numbands)
    eps = np.zeros(numbands)
    epsi = np.zeros(numbands)
    emi = np.zeros(numbands, dtype=bool)
    tau0c = np.zeros(numbands)
    bandcount = 0
    for bands in list(data.keys()):
        banddata = data[bands]
        bandstring = bands.split()[1].split("-")
        # single band parameters
        if len(bandstring) == 1:
            band = int(bandstring[0])
            bandparams[band][0] = 10
            bandparams[band][1] = 0
            effmass[band] = np.array(banddata["effmass"])
            select_scattering[band] = np.array(banddata["select_scattering"])
            q_energy_trans[band] = np.array(banddata["q_energy_trans"])
            explicit_prefact[band] = np.array(banddata["explicit_prefact"])
            explicit_prefact_values[band] = np.array(
                banddata["explicit_prefact_values"])
            da[band] = banddata["d_a"]
            do[band] = banddata["d_o"]
            speed_sound[band] = banddata["speed_sound"]
            no[band] = banddata["n_o"]
            nvv[band] = banddata["n_vv"]
            ni[band] = banddata["n_i"]
            omegao[band] = banddata["omega_o"]
            omegavv[band] = banddata["omega_vv"]
            etrans[band] = banddata["etrans"]
            rho[band] = banddata["rho"]
            zf[band] = banddata["zf"]
            f[band] = banddata["f"]
            z[band] = banddata["z"]
            isl[band] = banddata["isl"]
            isli[band] = banddata["isl_i"]
            vdiff[band] = banddata["vdiff"]
            alloyconc[band] = banddata["alloyconc"]
            eps[band] = banddata["eps"]
            epsi[band] = banddata["epsi"]
            p[band] = banddata["p"]
            emi[band] = banddata["emission"]
            tau0c[band] = banddata["tau0_c"]
            bandcount += 1
        # multi band parameters
        else:
            band_lower = int(bandstring[0]) - 1
            # no upper limit
            if bandstring[1] == '':
                band_upper = numbands
            # range
            else:
                band_upper = int(bandstring[1])
            for band in range(band_lower, band_upper):
                bandcount += 1
                bandparams[band][0] = 10
                bandparams[band][1] = 0
                effmass[band] = np.array(banddata["effmass"])
                select_scattering[band] = np.array(
                    banddata["select_scattering"])
                q_energy_trans[band] = np.array(banddata["q_energy_trans"])
                explicit_prefact[band] = np.array(banddata["explicit_prefact"])
                explicit_prefact_values[band] = np.array(
                    banddata["explicit_prefact_values"])
                da[band] = banddata["d_a"]
                do[band] = banddata["d_o"]
                speed_sound[band] = banddata["speed_sound"]
                no[band] = banddata["n_o"]
                nvv[band] = banddata["n_vv"]
                ni[band] = banddata["n_i"]
                omegao[band] = banddata["omega_o"]
                omegavv[band] = banddata["omega_vv"]
                etrans[band] = banddata["etrans"]
                zf[band] = banddata["zf"]
                f[band] = banddata["f"]
                z[band] = banddata["z"]
                rho[band] = banddata["rho"]
                eps[band] = banddata["eps"]
                epsi[band] = banddata["epsi"]
                isl[band] = banddata["isl"]
                isli[band] = banddata["isl_i"]
                vdiff[band] = banddata["vdiff"]
                alloyconc[band] = banddata["alloyconc"]
                p[band] = banddata["p"]
                emi[band] = banddata["emission"]
                tau0c[band] = banddata["tau0_c"]
    # now check if analytic scattering is set and print warning
    # (bands from VASP is seldom parabolic)
    if bs.param.transport_use_analytic_scattering:
        logger.warning("transport_use_analytic is set, but data from "
                       "VASP is processed. These bands are seldom analytic "
                       "and the analytic scattering models can then not "
                       "be used. Hope you know what your are doing. "
                       "We recommend setting transport_use_analytic to "
                       "False in order to invoke the numeric routines to "
                       "calculate the carrier scattering. "
                       "Also, in this case "
                       "the program might use an extreme amount of memory "
                       "since the scattering arrays are stored for all "
                       "possible mechanisms, not just the enabled ones. "
                       "Continuing.")
    if bandcount != numbands:
        logger.error(
            "The number of band parameters in bandparams.yml does not "
            "correspond to the number of bands supplied in vasprun.xml. "
            "Please correct bandparams.yml and rerun. Exiting.")
        sys.exit(1)
    bs.bandparams = bandparams
    bs.effmass = effmass
    bs.q_energy_trans = q_energy_trans
    bs.da = da
    bs.do = do
    speed_sound[np.abs(speed_sound) < constants.zero] = constants.zero
    bs.speed_sound = speed_sound
    bs.no = no
    bs.nvv = nvv
    ni[np.abs(ni) < constants.zero] = constants.zero
    bs.ni = ni
    omegao[np.abs(omegao) < constants.zero] = constants.zero
    bs.omegao = omegao
    omegavv[np.abs(omegavv) < constants.zero] = constants.zero
    bs.omegavv = omegavv
    bs.etrans = etrans
    rho[np.abs(rho) < constants.zero] = constants.zero
    bs.rho = rho
    bs.emi = emi
    bs.a = None
    bs.f = f
    bs.e0 = None
    bs.tau0c = tau0c
    bs.zf = zf
    z[np.abs(z) < constants.zero] = constants.zero
    bs.z = z
    isl[np.abs(isl) < constants.zero] = constants.zero
    bs.isl = isl
    bs.isli = isli
    bs.vdiff = vdiff
    bs.alloyconc = alloyconc
    eps[np.abs(eps) < constants.zero] = constants.zero
    bs.eps = eps
    epsi[np.abs(epsi) < constants.zero] = constants.zero
    bs.epsi = epsi
    bs.p = p
    bs.status = None
    bs.kshift = None
    bs.select_scattering = (select_scattering == 1)
    bs.explicit_prefact = explicit_prefact
    bs.explicit_prefact_values = explicit_prefact_values
    # vasp includes a volume factor, remove this
    # so that the units are 1/eVAA^3 (needed for scattering
    # routines etc).
    volume = np.linalg.det(bs.lattice.unitcell)
    bs.dos = dosvasp[:, 1] / volume
    bs.dos_partial = None
    bs.dos_energies = dosvasp[:, 0]
    bs.energies = energiesvasp
    bs.velocities = velocitiesvasp
    bs.occ = occvasp
    bs.tight_hop = None
    bs.tight_orb = None
    bs.tight_onsite = None
    bs.tight_adj_onsite = None
    # set flag to generate velocities later
    if bs.velocities is None:
        bs.gen_velocities = True
        # also generate velocity arrays
        bs.velocities = np.zeros((numbands, 3, bs.lattice.kmesh.shape[0]),
                                 dtype=np.double)
    else:
        bs.gen_velocities = False


def bandstructure_numpy(bs, filename, location=None):  # pylint: disable=too-many-locals # noqa: MC0001
    """
    Sets the bandstructure from a NumPy datafile file.

    Loads and stores the parameters in the bandstructure configuration file (defaults to bandparam.yml).

    Parameters
    ----------
    bs : object
        A `Bandstructure()` object.
    filename : string
        The filename of the NumPy data file to be read.
        The bandstructure
        configuration file have to be named "bandparam.yml"
        in this case.
    location : string, optional
        The location of the NumPy data
        file. Defaults to the "input" directory in
        the current working directory.

    Returns
    -------
    None

    Notes
    -----
    This routine read NumPy datafiles containing the
    electron energy dispersions and optionally the band
    velocities.

    | The datastructure of the supplied numpy array
    | should be on the following format:
    | [
    | [kx], [ky], [kz], [e_1], [v_x_1], [v_y_1], [v_z_1],
    | [e_2], [v_x_2], [v_y_2], [v_z_2], ... ,
    | [e_n], [v_x_n], [v_y_n], [v_z_n]
    | ]

    If the band velocities are not supplied they are simply not
    present. Each column of data has the length of the number of
    k-point in the full BZ.

    The bandstructure configuration file is still read due to
    the need of the scattering properties etc.

    This interface is enabled by setting `read` in the general
    configuration file to "numpy" (datafile with only
    electron energy dispersions) or "numpyv" (datafile with
    electron energy and group velocity dispersion)

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running bandstructure_numpy.")

    if filename is None:
        # check param file
        if bs.param.readfile == "":
            logger.error("User have to specificy readfile in order to read "
                         "Numpy data files. Exiting.")
            sys.exit(1)
        else:
            numpy_filename = "input/" + bs.param.readfile
    else:
        if location is not None:
            numpy_filename = location + "/" + filename
        else:
            numpy_filename = "input/" + filename
    data = np.load(numpy_filename)
    # overwrite the ksampling and unit vectors that are calculated based
    # on the YAML parameter file for celldata
    ksampling = np.zeros(3, dtype="intc")
    for i in range(3):
        ksampling[i] = data.shape[i + 1]
    bs.lattice.ksampling = ksampling
    for i in range(3):
        bs.lattice.ksampling[i] = data.shape[i + 1]

    numbands = data.shape[0] - 3
    # set flag to generate velocities later
    if bs.param.read != "numpyv":
        bs.gen_velocities = True
        # also generate velocity arrays
        bs.velocities = np.zeros((numbands, 3, bs.lattice.kmesh.shape[0]),
                                 dtype=np.double)
        numbands = data.shape[0] - 3
        # silly but this is going to be reworked anyway, should use
        # broadcasting
        for band in range(numbands):
            bs.energies[band] = data[3 + band].flatten()
    else:
        bs.gen_velocities = False
        numbands = numbands / 4
    numkpoints = np.prod(bs.lattice.ksampling)
    bs.energies = np.zeros((numbands, numkpoints))
    bs.velocities = np.zeros((numbands, 3, numkpoints))
    bs.lattice.kmesh = np.zeros((numkpoints, 3))
    bs.lattice.kmesh[:, 0] = data[0].flatten()
    bs.lattice.kmesh[:, 1] = data[1].flatten()
    bs.lattice.kmesh[:, 2] = data[2].flatten()
    # silly but this is going to be reworked anyway, should use
    # broadcasting
    for band in range(numbands):
        bs.energies[band] = data[3 + band].flatten()
        if bs.param.read == "numpyv":
            for x in range(3):
                bs.velocities[band][x] = data[3 + band +
                                              (x + 1) * numbands].flatten()

    # what about energy shifts? if the user set any of the parameters
    # relevant for the Fermi level, give warning and continue
    if bs.param.e_fermi or bs.param.e_fermi_in_gap or bs.param.e_vbm:
        logger.info("User have set 'e_fermi' and/or 'e_fermi_in_gap' "
                    "and/or 'e_vbm' to True, but this options are not "
                    "supported for the reading of Numpy datafiles. "
                    "Shifting the bandstructure by 'e_shift' and "
                    "continuing.")
        # velocities are of course independent of this shift
        bs.energies = bs.energies - bs.param.e_shift

    spin_degen = np.zeros(numbands, dtype="intc")
    # bandparams contain scattering properties etc.
    data = inputoutput.readbandparam(location, filename="bandparam.yml")
    bandparams = np.zeros((numbands, 2), dtype=np.int8)
    effmass = np.zeros((numbands, 3))
    select_scattering = np.zeros((numbands, 12), dtype='intc')
    explicit_prefact = np.zeros((numbands, 11), dtype='intc')
    explicit_prefact_values = np.zeros((numbands, 11))
    q_energy_trans = np.zeros((numbands, 2, 3))
    da = np.zeros(numbands)
    do = np.zeros(numbands)
    speed_sound = np.zeros(numbands)
    no = np.zeros(numbands)
    nvv = np.zeros(numbands)
    ni = np.zeros(numbands)
    omegao = np.zeros(numbands)
    omegavv = np.zeros(numbands)
    etrans = np.zeros(numbands)
    rho = np.zeros(numbands)
    zf = np.zeros(numbands)
    f = np.zeros(numbands)
    z = np.zeros(numbands)
    isl = np.zeros(numbands)
    isli = np.zeros(numbands)
    vdiff = np.zeros(numbands)
    alloyconc = np.zeros(numbands)
    p = np.zeros(numbands)
    eps = np.zeros(numbands)
    epsi = np.zeros(numbands)
    emi = np.zeros(numbands, dtype=bool)
    tau0c = np.zeros(numbands)
    bandcount = 0
    for bands in list(data.keys()):
        banddata = data[bands]
        bandstring = bands.split()[1].split("-")
        # single band parameters
        if len(bandstring) == 1:
            band = int(bandstring[0])
            bandparams[band][0] = 10
            bandparams[band][1] = 0
            effmass[band] = np.array(banddata["effmass"])
            select_scattering[band] = np.array(banddata["select_scattering"])
            q_energy_trans[band] = np.array(banddata["q_energy_trans"])
            explicit_prefact[band] = np.array(banddata["explicit_prefact"])
            explicit_prefact_values[band] = np.array(
                banddata["explicit_prefact_values"])
            da[band] = banddata["d_a"]
            do[band] = banddata["d_o"]
            speed_sound[band] = banddata["speed_sound"]
            no[band] = banddata["n_o"]
            nvv[band] = banddata["n_vv"]
            ni[band] = banddata["n_i"]
            omegao[band] = banddata["omega_o"]
            omegavv[band] = banddata["omega_vv"]
            etrans[band] = banddata["etrans"]
            rho[band] = banddata["rho"]
            zf[band] = banddata["zf"]
            f[band] = banddata["f"]
            z[band] = banddata["z"]
            isl[band] = banddata["isl"]
            isli[band] = banddata["isl_i"]
            vdiff[band] = banddata["vdiff"]
            alloyconc[band] = banddata["alloyconc"]
            eps[band] = banddata["eps"]
            epsi[band] = banddata["epsi"]
            p[band] = banddata["p"]
            emi[band] = banddata["emission"]
            spin_degen[band] = banddata["spin_degen"]
            tau0c[band] = banddata["tau0_c"]
            bandcount += 1
        # multi band parameters
        else:
            band_lower = int(bandstring[0]) - 1
            # no upper limit
            if bandstring[1] == '':
                band_upper = numbands
                # range
            else:
                band_upper = int(bandstring[1])
            for band in range(band_lower, band_upper):
                bandcount += 1
                bandparams[band][0] = 10
                bandparams[band][1] = 0
                effmass[band] = np.array(banddata["effmass"])
                select_scattering[band] = np.array(
                    banddata["select_scattering"])
                q_energy_trans[band] = np.array(banddata["q_energy_trans"])
                explicit_prefact[band] = np.array(banddata["explicit_prefact"])
                explicit_prefact_values[band] = np.array(
                    banddata["explicit_prefact_values"])
                da[band] = banddata["d_a"]
                do[band] = banddata["d_o"]
                speed_sound[band] = banddata["speed_sound"]
                no[band] = banddata["n_o"]
                nvv[band] = banddata["n_vv"]
                ni[band] = banddata["n_i"]
                omegao[band] = banddata["omega_o"]
                omegavv[band] = banddata["omega_vv"]
                etrans[band] = banddata["etrans"]
                zf[band] = banddata["zf"]
                f[band] = banddata["f"]
                z[band] = banddata["z"]
                rho[band] = banddata["rho"]
                eps[band] = banddata["eps"]
                epsi[band] = banddata["epsi"]
                isl[band] = banddata["isl"]
                isli[band] = banddata["isl_i"]
                vdiff[band] = banddata["vdiff"]
                alloyconc[band] = banddata["alloyconc"]
                p[band] = banddata["p"]
                spin_degen[band] = banddata["spin_degen"]
                emi[band] = banddata["emission"]
                tau0c[band] = banddata["tau0_c"]
        # now check if analytic scattering is set and print warning
        # (bands from numpy is seldom parabolic)
    if bs.param.transport_use_analytic_scattering:
        logger.warning("transport_use_analytic is set, but data from "
                       "Numpy is processed. These bands are seldom analytic "
                       "and the analytic scattering models can then not "
                       "be used. Hope you know what your are doing. We "
                       "recommend setting transport_use_analytic to False "
                       "in order to invoke the numeric routines to "
                       "calculate the carrier scattering. Continuing.")
    if bandcount != numbands:
        logger.error(
            "The number of band parameters in bandparams.yml does not "
            "correspond to the number of bands supplied in the Numpy "
            "file. Please correct bandparams.yml and rerun. Exiting.")
        sys.exit(1)
    bs.bandparams = bandparams
    bs.effmass = effmass
    bs.q_energy_trans = q_energy_trans
    bs.da = da
    bs.do = do
    speed_sound[np.abs(speed_sound) < constants.zero] = constants.zero
    bs.speed_sound = speed_sound
    bs.no = no
    bs.nvv = nvv
    ni[np.abs(ni) < constants.zero] = constants.zero
    bs.ni = ni
    omegao[np.abs(omegao) < constants.zero] = constants.zero
    bs.omegao = omegao
    omegavv[np.abs(omegavv) < constants.zero] = constants.zero
    bs.omegavv = omegavv
    bs.etrans = etrans
    rho[np.abs(rho) < constants.zero] = constants.zero
    bs.rho = rho
    bs.emi = emi
    bs.a = None
    bs.f = f
    bs.e0 = None
    bs.tau0c = tau0c
    bs.zf = zf
    z[np.abs(z) < constants.zero] = constants.zero
    bs.z = z
    isl[np.abs(isl) < constants.zero] = constants.zero
    bs.isl = isl
    bs.isli = isli
    bs.vdiff = vdiff
    bs.alloyconc = alloyconc
    eps[np.abs(eps) < constants.zero] = constants.zero
    bs.eps = eps
    epsi[np.abs(epsi) < constants.zero] = constants.zero
    bs.epsi = epsi
    bs.p = p
    bs.status = None
    bs.kshift = None
    bs.spin_degen = spin_degen
    bs.select_scattering = (select_scattering == 1)
    bs.explicit_prefact = explicit_prefact
    bs.explicit_prefact_values = explicit_prefact_values
    bs.dos_partial = None
    bs.tight_hop = None
    bs.tight_orb = None
    bs.tight_onsite = None
    bs.tight_adj_onsite = None
    # upon reading Numpy datafiles, we have
    # no access to occupancies, so set the following
    # to None
    bs.vbm_energy = None
    bs.vbm_band = None
    bs.vbm_kpoint = None
    bs.cbm_energy = None
    bs.cbm_band = None
    bs.cbm_kpoint = None
    bs.band_gap = None
    bs.direct = None


def bandstructure_w90(bs, location=None, filename=None):  # pylint: disable=too-many-locals, too-many-branches
    """
    Sets the bandstructure from a Wannier90 output file.

    This is fed into PythTB to reconstruct and (extrapolate
    and interpolate the bandstructure on as dense grid as
    necessary). It also loads and stores the parameters in
    the bandstructure configuration file (defaults to bandparam.yml).

    Parameters
    ----------
    bs : object
        A `Bandstructure()` object.
    location : string, optional
        The location of the Wannier90 input and output files.
        Defaults to the "input" directory in the current
        working directory.
    filename : string, optional
        The prefix of the Wannier90 input and output files
        to be read. Defaults to "wannier90" (default for VASP output).
        The bandstructure configuration file have to be named
        "bandparam.yml" in this case.

    Returns
    -------
    None

    Notes
    -----
    This interface load input and output files from a Wannier90
    calculation (which again can be done post-first-principles).
    This is then fed into PythTB (if present) in order to reconstruct
    the electronic structure and interpolate the electron energy
    dispersions on as dense grid as needed. Presently there is no
    way to generate the band velocities so these has to be extracted
    by an interpolation routine later. Flags are automatically set
    for this. Please consult the manual of
    `PythTB <http://www.physics.rutgers.edu/pythtb/>`_ for how to
    prepare the correct Wannier90 input files and produce the
    necessary output files.

    This interface is enabled by setting `read` in the general
    configuration file to "w90".

    .. warning:: Please be aware of the importance of knowing what
                 the Wannier orbitals are at the end of the Wannier90
                 run.

    """
    # lazy import of PythTB (optional)
    import pythtb  # pylint: disable=import-error

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running bandstructure_w90.")

    if filename is None:
        # check param file
        if bs.param.readfile == "":
            wprefix = "wannier90"
        else:
            wprefix = bs.param.readfile
    else:
        wprefix = filename

    wannierdata = pythtb.w90("input", wprefix)
    tb = wannierdata.model(
        zero_energy=bs.param.dispersion_w90_tb_zero_energy,
        min_hopping_norm=bs.param.dispersion_w90_tb_min_hopping_norm,
        max_distance=bs.param.dispersion_w90_tb_max_distance)

    # do we want the energies on the original grid or a denser one?
    if (bs.param.dispersion_interpolate
            and (bs.param.dispersion_interpolate_method == "tb")):
        logger.info("Detected that the user want extract the energies "
                    "at a denser grid than what was supplied when "
                    "constructing the Wannier model. Switching the "
                    "grid to the target grid.")
        # fetch denser grid sampling
        iksampling = bs.lattice.fetch_iksampling()
        # regenerate grid and store
        bs.lattice.create_kmesh(iksampling, borderless=True)

    kmesh = bs.lattice.kmesh + 0.5
    energies = tb.solve_all(kmesh)
    # now, we do not have the velocities, but let us
    # not care about that for now (similar to VASP)
    numbands = energies.shape[0]
    numkpoints = np.prod(bs.lattice.ksampling)
    bs.energies = energies
    bs.velocities = np.zeros((numbands, 3, numkpoints))
    bs.gen_velocities = True
    # what about energy shifts? if the user set any of the parameters
    # relevant for the Fermi level, give warning and continue
    if bs.param.e_fermi or bs.param.e_fermi_in_gap or bs.param.e_vbm:
        logger.info("User have set 'e_fermi' and/or 'e_fermi_in_gap' "
                    "and/or 'e_vbm' to True, but this options are not "
                    "supported for the reading of Wannier90 datafiles. "
                    "Shifting the bandstructure by 'e_shift' and "
                    "continuing.")
        # velocities are of course independent of this shift
        bs.energies = bs.energies - bs.param.e_shift

    # now read the bandparam file
    spin_degen = np.zeros(numbands, dtype="intc")
    # bandparams contain scattering properties etc.
    data = inputoutput.readbandparam(location, filename="bandparam.yml")
    bandparams = np.zeros((numbands, 2), dtype=np.int8)
    effmass = np.zeros((numbands, 3))
    select_scattering = np.zeros((numbands, 12), dtype='intc')
    explicit_prefact = np.zeros((numbands, 11), dtype='intc')
    explicit_prefact_values = np.zeros((numbands, 11))
    q_energy_trans = np.zeros((numbands, 2, 3))
    da = np.zeros(numbands)
    do = np.zeros(numbands)
    speed_sound = np.zeros(numbands)
    no = np.zeros(numbands)
    nvv = np.zeros(numbands)
    ni = np.zeros(numbands)
    omegao = np.zeros(numbands)
    omegavv = np.zeros(numbands)
    etrans = np.zeros(numbands)
    rho = np.zeros(numbands)
    zf = np.zeros(numbands)
    f = np.zeros(numbands)
    z = np.zeros(numbands)
    isl = np.zeros(numbands)
    isli = np.zeros(numbands)
    vdiff = np.zeros(numbands)
    alloyconc = np.zeros(numbands)
    p = np.zeros(numbands)
    eps = np.zeros(numbands)
    epsi = np.zeros(numbands)
    emi = np.zeros(numbands, dtype=bool)
    tau0c = np.zeros(numbands)
    bandcount = 0
    for bands in list(data.keys()):
        banddata = data[bands]
        bandstring = bands.split()[1].split("-")
        # single band parameters
        if len(bandstring) == 1:
            band = int(bandstring[0])
            bandparams[band][0] = 10
            bandparams[band][1] = 0
            effmass[band] = np.array(banddata["effmass"])
            select_scattering[band] = np.array(banddata["select_scattering"])
            q_energy_trans[band] = np.array(banddata["q_energy_trans"])
            explicit_prefact[band] = np.array(banddata["explicit_prefact"])
            explicit_prefact_values[band] = np.array(
                banddata["explicit_prefact_values"])
            da[band] = banddata["d_a"]
            do[band] = banddata["d_o"]
            speed_sound[band] = banddata["speed_sound"]
            no[band] = banddata["n_o"]
            nvv[band] = banddata["n_vv"]
            ni[band] = banddata["n_i"]
            omegao[band] = banddata["omega_o"]
            omegavv[band] = banddata["omega_vv"]
            etrans[band] = banddata["etrans"]
            rho[band] = banddata["rho"]
            zf[band] = banddata["zf"]
            f[band] = banddata["f"]
            z[band] = banddata["z"]
            isl[band] = banddata["isl"]
            isli[band] = banddata["isl_i"]
            vdiff[band] = banddata["vdiff"]
            alloyconc[band] = banddata["alloyconc"]
            eps[band] = banddata["eps"]
            epsi[band] = banddata["epsi"]
            p[band] = banddata["p"]
            emi[band] = banddata["emission"]
            spin_degen[band] = banddata["spin_degen"]
            tau0c[band] = banddata["tau0_c"]
            bandcount += 1
            # multi band parameters
        else:
            band_lower = int(bandstring[0]) - 1
            # no upper limit
            if bandstring[1] == '':
                band_upper = numbands
            # range
            else:
                band_upper = int(bandstring[1])
            for band in range(band_lower, band_upper):
                bandcount += 1
                bandparams[band][0] = 10
                bandparams[band][1] = 0
                effmass[band] = np.array(banddata["effmass"])
                select_scattering[band] = np.array(
                    banddata["select_scattering"])
                q_energy_trans[band] = np.array(banddata["q_energy_trans"])
                explicit_prefact[band] = np.array(banddata["explicit_prefact"])
                explicit_prefact_values[band] = np.array(
                    banddata["explicit_prefact_values"])
                da[band] = banddata["d_a"]
                do[band] = banddata["d_o"]
                speed_sound[band] = banddata["speed_sound"]
                no[band] = banddata["n_o"]
                nvv[band] = banddata["n_vv"]
                ni[band] = banddata["n_i"]
                omegao[band] = banddata["omega_o"]
                omegavv[band] = banddata["omega_vv"]
                etrans[band] = banddata["etrans"]
                zf[band] = banddata["zf"]
                f[band] = banddata["f"]
                z[band] = banddata["z"]
                rho[band] = banddata["rho"]
                eps[band] = banddata["eps"]
                epsi[band] = banddata["epsi"]
                isl[band] = banddata["isl"]
                isli[band] = banddata["isl_i"]
                vdiff[band] = banddata["vdiff"]
                alloyconc[band] = banddata["alloyconc"]
                p[band] = banddata["p"]
                spin_degen[band] = banddata["spin_degen"]
                emi[band] = banddata["emission"]
                tau0c[band] = banddata["tau0_c"]
    # now check if analytic scattering is set and print warning
    # (bands from w90 is seldom parabolic)
    if bs.param.transport_use_analytic_scattering:
        logger.warning("transport_use_analytic is set, but data from "
                       "Numpy is processed. These bands are seldom "
                       "analytic and the analytic scattering models "
                       "can then not be used. Hope you know what your "
                       "are doing. We recommend setting "
                       "transport_use_analytic to False in order to "
                       "invoke the numeric routines to calculate the "
                       "carrier scattering. Continuing.")
    if bandcount != numbands:
        logger.error(
            "The number of band parameters in bandparams.yml does not "
            "correspond to the number of bands supplied in the "
            "Wannier90 .win file. Please correct bandparams.yml and "
            "rerun. Exiting.")
        sys.exit(1)
    bs.bandparams = bandparams
    bs.effmass = effmass
    bs.q_energy_trans = q_energy_trans
    bs.da = da
    bs.do = do
    speed_sound[np.abs(speed_sound) < constants.zero] = constants.zero
    bs.speed_sound = speed_sound
    bs.no = no
    bs.nvv = nvv
    ni[np.abs(ni) < constants.zero] = constants.zero
    bs.ni = ni
    omegao[np.abs(omegao) < constants.zero] = constants.zero
    bs.omegao = omegao
    omegavv[np.abs(omegavv) < constants.zero] = constants.zero
    bs.omegavv = omegavv
    bs.etrans = etrans
    rho[np.abs(rho) < constants.zero] = constants.zero
    bs.rho = rho
    bs.emi = emi
    bs.a = None
    bs.f = f
    bs.e0 = None
    bs.tau0c = tau0c
    bs.zf = zf
    z[np.abs(z) < constants.zero] = constants.zero
    bs.z = z
    isl[np.abs(isl) < constants.zero] = constants.zero
    bs.isl = isl
    bs.isli = isli
    bs.vdiff = vdiff
    bs.alloyconc = alloyconc
    eps[np.abs(eps) < constants.zero] = constants.zero
    bs.eps = eps
    epsi[np.abs(epsi) < constants.zero] = constants.zero
    bs.epsi = epsi
    bs.p = p
    bs.status = None
    bs.kshift = None
    bs.spin_degen = spin_degen
    bs.select_scattering = (select_scattering == 1)
    bs.explicit_prefact = explicit_prefact
    bs.explicit_prefact_values = explicit_prefact_values
    bs.dos_partial = None
    bs.tight_hop = None
    bs.tight_orb = None
    bs.tight_onsite = None
    bs.tight_adj_onsite = None
    # upon reading Wannier90 datafiles, we have
    # no access to occupancies, so set the following
    # to None
    bs.vbm_energy = None
    bs.vbm_band = None
    bs.vbm_kpoint = None
    bs.cbm_energy = None
    bs.cbm_band = None
    bs.cbm_kpoint = None
    bs.band_gap = None
    bs.direct = None


def read_band_parameters(bs, numbands, location=None, filename=None):  # pylint: disable=too-many-locals
    """
    Reads and stores the information in the band parameters configuration file (bandparam.yml).

    Parameters
    ----------
    bs : object
        The active `Bandstructure()` object.
    numbands : int
        The number of bands
    location : string, optional
        The folder in which the band configuration
        file is placed. Defaults to the relative
        folder "input".
    filename : string, optional
        The filename of the band configuration file.
        Defaults to bandparam.yml.

    Returns
    -------
    None

    Notes
    -----
    Reads and stores the values in the band configuration file.
    When writing custom interfaces it is sufficient to call this
    routine in order for the setup of the individual band
    parameters to be consistent.

    """

    # set logger
    logger = logging.getLogger(sys._getframe().f_code.co_name)  # pylint: disable=protected-access
    logger.debug("Running read_band_parameters.")

    # read band parameter file
    data = inputoutput.readbandparam(location, filename)
    numbandparams = len(list(data.values()))
    numbands = 0
    # now check if we might run into tight binding generation
    # (this can squeeze more bands into the range, so wait
    # with the array declarations etc.)
    # the number of tight binding orbitals is set always
    # len(torb), which defaults to the number of positions
    # in the Lattice() object.
    # this could probably be done a bit smarter
    bandparam_for_band = []
    for bandparam in range(numbandparams):
        try:
            bandparamdata = data["Band " + str(bandparam + 1)]
        except KeyError:
            logger.error("The Band X segments in bandparam.yml is not to spec "
                         "or you have not specified the right read flag in "
                         "param.yml. Maybe you are switching to/from a read "
                         "from external data and param reading and forgot to "
                         " also modify bandparams.yml? Exiting.")
            sys.exit(1)
        # check for tight binding entries
        if bandparamdata["type"] == 3:
            numtbbands = len(bandparamdata["torb"])
            if numtbbands == 0:
                numtbbands = bs.lattice.positions.shape[0]
            for _ in range(numtbbands):
                bandparam_for_band.append(bandparam)
            numbands = numbands + numtbbands
        else:
            numbands = numbands + 1
            bandparam_for_band.append(bandparam)

    # declare the arrays
    bandparams = np.zeros((numbands, 2), dtype=np.int8)
    effmass = np.zeros((numbands, 3))
    a = np.zeros((numbands, 3))
    ascale = np.zeros(numbands)
    e0 = np.zeros(numbands)
    status = np.empty(numbands, dtype='str')
    kshift = np.zeros((numbands, 3))
    select_scattering = np.zeros((numbands, 12), dtype='intc')
    explicit_prefact = np.zeros((numbands, 11), dtype='intc')
    explicit_prefact_values = np.zeros((numbands, 11))
    q_energy_trans = np.zeros((numbands, 2, 3))
    da = np.zeros(numbands)
    do = np.zeros(numbands)
    speed_sound = np.zeros(numbands)
    no = np.zeros(numbands)
    nvv = np.zeros(numbands)
    ni = np.zeros(numbands)
    omegao = np.zeros(numbands)
    omegavv = np.zeros(numbands)
    etrans = np.zeros(numbands)
    rho = np.zeros(numbands)
    zf = np.zeros(numbands)
    z = np.zeros(numbands)
    eps = np.zeros(numbands)
    f = np.zeros(numbands)
    epsi = np.zeros(numbands)
    p = np.zeros(numbands)
    isl = np.zeros(numbands)
    isli = np.zeros(numbands)
    tau0c = np.zeros(numbands)
    vdiff = np.zeros(numbands)
    alloyconc = np.zeros(numbands)
    emi = np.zeros(numbands, dtype=bool)
    spin_degen = np.zeros(numbands, dtype="intc")
    tight_hop = []
    tight_orb = []
    tight_onsite = []
    tight_adj_onsite = []
    for band in range(numbands):
        try:
            banddata = data["Band " + str(bandparam_for_band[band] + 1)]
        except KeyError:
            logger.error("The Band X segments in bandparam.yml is not to spec "
                         "or you have not specified the right read flag in "
                         "param.yml. Maybe you are switching to/from a read "
                         "from external data and param reading and forgot to "
                         " also modify bandparams.yml? Exiting.")
            sys.exit(1)
        bandparams[band] = np.array([banddata["type"], banddata["folding"]])
        effmass[band] = np.array(banddata["effmass"])
        a[band] = np.array(banddata["a"])
        ascale[band] = banddata["ascale"]
        e0[band] = banddata["e0"]
        status[band] = banddata["status"]
        kshift[band] = np.array(banddata["kshift"])
        select_scattering[band] = np.array(banddata["select_scattering"])
        explicit_prefact[band] = np.array(banddata["explicit_prefact"])
        explicit_prefact_values[band] = np.array(
            banddata["explicit_prefact_values"])
        q_energy_trans[band] = np.array(banddata["q_energy_trans"])
        da[band] = banddata["d_a"]
        do[band] = banddata["d_o"]
        speed_sound[band] = banddata["speed_sound"]
        no[band] = banddata["n_o"]
        nvv[band] = banddata["n_vv"]
        ni[band] = banddata["n_i"]
        omegao[band] = banddata["omega_o"]
        omegavv[band] = banddata["omega_vv"]
        etrans[band] = banddata["etrans"]
        rho[band] = banddata["rho"]
        zf[band] = banddata["zf"]
        f[band] = banddata["f"]
        z[band] = banddata["z"]
        eps[band] = banddata["eps"]
        epsi[band] = banddata["epsi"]
        p[band] = banddata["p"]
        isl[band] = banddata["isl"]
        isli[band] = banddata["isl_i"]
        vdiff[band] = banddata["vdiff"]
        alloyconc[band] = banddata["alloyconc"]
        tau0c[band] = banddata["tau0_c"]
        emi[band] = banddata["emission"]
        spin_degen[band] = banddata["spin_degen"]
        tight_hop.append(banddata["thop"])
        tight_orb.append(banddata["torb"])
        tight_onsite.append(banddata["tonsite"])
        # if band is non-parabolic, print warning
        # if transport_use_analytic_scattering is also set
        if (bandparams[band][0] != 0
                and bs.param.transport_use_analytic_scattering):
            logger.warning("The supplied band is non-parabolic, while "
                           "transport_use_analytic_scattering is set to True. "
                           "The analytic scattering models are currently only "
                           "valid for parabolic bands. Continuing.")
    bs.bandparams = bandparams
    bs.effmass = effmass
    bs.q_energy_trans = q_energy_trans
    bs.da = da
    bs.do = do
    speed_sound[np.abs(speed_sound) < constants.zero] = constants.zero
    bs.speed_sound = speed_sound
    bs.no = no
    bs.nvv = nvv
    ni[np.abs(ni) < constants.zero] = constants.zero
    bs.ni = ni
    omegao[np.abs(omegao) < constants.zero] = constants.zero
    bs.omegao = omegao
    omegavv[np.abs(omegavv) < constants.zero] = constants.zero
    bs.omegavv = omegavv
    rho[np.abs(rho) < constants.zero] = constants.zero
    bs.rho = rho
    bs.tau0c = tau0c
    bs.emi = emi
    bs.f = f
    bs.a = a
    bs.ascale = ascale
    bs.e0 = e0
    bs.status = status
    bs.kshift = kshift
    bs.select_scattering = (select_scattering == 1)
    bs.explicit_prefact = explicit_prefact
    bs.explicit_prefact_values = explicit_prefact_values
    bs.zf = zf
    z[np.abs(z) < constants.zero] = constants.zero
    bs.z = z
    bs.etrans = etrans
    isl[np.abs(isl) < constants.zero] = constants.zero
    bs.isl = isl
    bs.isli = isli
    bs.vdiff = vdiff
    bs.alloyconc = alloyconc
    eps[np.abs(eps) < constants.zero] = constants.zero
    bs.eps = eps
    epsi[np.abs(epsi) < constants.zero] = constants.zero
    bs.epsi = epsi
    bs.p = p
    bs.spin_degen = spin_degen
    bs.tight_hop = tight_hop
    bs.tight_orb = np.array(tight_orb)
    bs.tight_onsite = np.array(tight_onsite)
    bs.tight_adj_onsite = np.array(tight_adj_onsite)
