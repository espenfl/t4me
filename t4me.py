#!/usr/bin/python
# python specifics
import sys
import logging
import numpy as np

# locals
import utils
import constants


def main():
    """
    Main routine.

    Sets up the calculations and calls necessary sub-routines.

    Notes
    -----
    This routine can be modified for more advanced usage. For
    day-day operations, the most used configurations are available
    by setting the parameters in the general configuration file (
    defaults to param.yml).

    """

    # check output folder (do this before logger config
    # in case user wants the log for the this job in the output
    # folder) and clean it
    utils.clean_directory("output")

    # configure logger
    utils.config_logger()
    logger = logging.getLogger(sys._getframe().f_code.co_name)

    # these modules are loaded after we configure the logger
    # in case there are problems with third party libraries etc.
    import bandstructure
    import lattice
    import transport
    import inputoutput

    # this is not so very clean, but in order to avoid
    # code running at the top of our modules (sphinx complains e.g.)
    # we here check if all modules can be imported, if not, prints
    # a warning and later we just load the module in a subroutine
    # and assume the user knows that is going on

    # check GSL and load
    try:
        import gsl
    except ImportError:
        inputoutput.gsl_error()

    # check spglib and load
    try:
        import spglib_interface
    except ImportError:
        inputoutput.spglib_error()

    # check PythTB and load
    try:
        import pythtb
    except ImportError:
        inputoutput.pythtb_warning()

    # check skw_interface and load
    try:
        import skw_interface
    except ImportError:
        inputoutput.skw_warning()

    # check wildmagic and load
    try:
        import wildmagic
    except ImportError:
        inputoutput.wildmagic_warning()

    # check einspline and load
    try:
        import einspline
    except ImportError:
        inputoutput.einspline_warning()

    # check cubature_wildmagic and load
    try:
        import cubature_wildmagic
    except ImportError:
        inputoutput.cubature_warning()

    # read param file
    param = inputoutput.Param(inputoutput.readparam())

    # optional, run tests for functionality
    if param.run_tests:
        import tests
        tests.run_tests(param.run_tests)
        sys.exit()

    # generate lattice from param and input file
    lat = lattice.Lattice(param)

    # generate or read bandstructure data
    bs = bandstructure.Bandstructure(lat, param)

    # generate the velocities if they are not present
    # and if the user does not want to use the extraction
    # of the interpolation routines
    # makes more sense to have this in lbtecoeff.py or
    # similar, but sometimes one would want to check the
    # velocities and we thus need them before checking if
    # user wants to dump dispersion relation (imcludes velocities
    # if present).
    if bs.gen_velocities and param.dispersion_velocities_numdiff \
       and param.transport_calc:
        bs.calc_velocities()

    # dump the dispersions?
    if param.dispersion_write_preinter:
        inputoutput.dump_bandstruct_line(bs, param.dispersion_write_start,
                                         param.dispersion_write_end,
                                         datatype="e", filename="bands",
                                         k_direct=True,
                                         itype="linearnd",
                                         itype_sub="linear")

        # same as above just for the velocities
        # check if the velocities exists, if not print a message and
        # continue
        if bs.gen_velocities:
            logger.info(
                "Data for the band velocities does not exist. "
                "Skipping writing the band velocities to file "
                "and continuing.")
        else:
            inputoutput.dump_bandstruct_line(bs, param.dispersion_write_start,
                                             param.dispersion_write_end,
                                             datatype="v", filename="velocities",
                                             k_direct=True,
                                             itype="linearnd",
                                             itype_sub="linear")

    # maybe the user wants to pre-interpolate?
    if param.dispersion_interpolate:
        logger.info("Pre-interpolating the dispersion data.")
        if bs.gen_velocities and param.transport_calc:
            # here we need the velocities
            bs.interpolate(store_inter=True, ivelocities=True)
            # if velocities did not exists, we now have them, so set
            bs.gen_velocities = False
        else:
            # no velocities needed here
            bs.interpolate(store_inter=True, ivelocities=False)

        # dump the dispersions after interpolations?
        if param.dispersion_write_postinter:
            inputoutput.dump_bandstruct_line(bs, param.dispersion_write_start,
                                             param.dispersion_write_end,
                                             datatype="e", filename="bands_inter",
                                             k_direct=True,
                                             itype="linearnd",
                                             itype_sub="linear")

            # same as above just for the velocities
            if bs.gen_velocities:
                logger.info(
                    "Data for the band velocities does not exist. "
                    "Skipping writing the band velocities to file "
                    "and continuing.")
            else:
                inputoutput.dump_bandstruct_line(bs, param.dispersion_write_start,
                                                 param.dispersion_write_end,
                                                 datatype="v", filename="velocities_inter",
                                                 k_direct=True,
                                                 itype="linearnd",
                                                 itype_sub="linear")

    sys.exit(1)

    # calculation of the density of states?
    if param.dos_calc:
        # calculate
        bs.calc_density_of_states(spin_degen=True)
        # write dos
        inputoutput.dump_density_of_states(bs, filename="dos")

    # transport calcs?
    if param.transport_calc:
        # initialize parameters, scattering etc.
        tran = transport.Transport(bs)
        # write relaxation times
        if not param.transport_method == "closed":
            inputoutput.dump_relaxation_time(tran)
        # perform the actual calculation of the transport coefficients
        tran.calc_transport_tensors()
        # write transport coefficients
        inputoutput.dump_transport_coefficients(tran)
        # dump dos as well
        if not param.transport_use_analytic_scattering:
            inputoutput.dump_density_of_states(bs, filename="dos")

if __name__ == '__main__':

    main()
