from __future__ import division

import numpy as np
from scipy.interpolate import interp1d


def _import_emis_file(rcp):
    if rcp in ['rcp3pd', 'rcp26']:
        from ..RCPs.rcp26 import Emissions as rcp_emis
    elif rcp=='rcp45':
        from ..RCPs.rcp45 import Emissions as rcp_emis
    elif rcp in ['rcp6', 'rcp60']:
        from ..RCPs.rcp60 import Emissions as rcp_emis
    elif rcp=='rcp85':
        from ..RCPs.rcp85 import Emissions as rcp_emis
    else:
        raise ValueError('rcp must be rcp26, rcp45, rcp60 or rcp85')
    return rcp_emis


def scen_open(filename,
              include_cfcs='rcp45',
              startyear=1765,
              harmonise=None):

    """
    Opens a MAGICC6 .SCEN file and extracts the data. Interpolates linearly
    between non-consecutive years in the SCEN file. Fills in chlorinated gases
    from a specified RCP scenario or from custom emissions.

    Inputs:
        filename: the .SCEN file to open

    Keywords:
        include_cfcs: string, False, or nt x 16 numpy array
            MAGICC files do not come loaded with CFCs (indices 24-39).
            They are given in the harmonised files at
            http://www.pik-potsdam.de/~mmalte/rcps/.
            - Specify 'rcp3pd', 'rcp45', 'rcp6' or 'rcp85' to use these RCPs.
            - Use False to ignore and create a 24-species emission file.
            - Provide an array to tack your own chlorinated gases onto the SCEN
        startyear: First year of output file. If before first year of the SCEN
            file, use RCP4.5 to fill
        harmonise: None, or year
            Linearly interpolate between 2000 in the RCP file and the specified
            year. If None, do not harmonise

    Returns:
        nt x 40 numpy emissions array

    nt is defined as <last year of SCEN file> -
        <earlier of startyear and first year of SCEN file> + 1

    It is assumed that the .SCEN files follow the formatting convention on
    the MAGICC wiki at
    http://wiki.magicc.org/index.php?title=Creating_MAGICC_Scenario_Files.
    """

    with open(filename) as f:
        # First line is the number of time steps in the SCEN file
        str_nt = f.readline().strip()
        nt = int(str_nt)

        # Next 6 lines are unused by FaIR.
        for i in range(6):
            f.readline()

        # Eighth line is the column headers. When object orientation is
        # implemented this will be important.
        headers = f.readline().split()

        # Ninth line is the units. Again this will be important in OO-FaIR
        units = f.readline().split()

    # Now the data!
    scen_emissions = np.genfromtxt(filename, skip_header=9, max_rows=nt)
    emissions = np.copy(scen_emissions)
    scen_years = scen_emissions[:,0]

    # Interpolate between non-consecutive years in SCEN file
    f = interp1d(scen_years, emissions[:,1:], kind='linear', axis=0,
        assume_sorted=True, bounds_error=True)
    full_years = np.arange(scen_years[0], scen_years[-1]+1)
    emissions_filled = np.hstack((full_years[:,None], f(full_years)))

    # Add CFCs if requested
    if type(include_cfcs) is np.ndarray:
        if include_cfcs.shape != (len(full_years), 16):
            raise ValueError("If chlorinated gas emissions are provided by " +
            "the user they should be of size nt x 16 where nt is the number "+
            "of consecutive years. Size is %s." % str(include_cfcs.shape))
        emissions_filled = np.append(emissions_filled, include_cfcs, axis=1)
    elif include_cfcs==False:
        pass
    elif include_cfcs.lower()[:3]=='rcp':
        rcp_emis = _import_emis_file(include_cfcs.lower()).emissions
        # Need to ensure only years present in the SCEN file are taken
        # Cheers: https://stackoverflow.com/questions/3522946/
        # using-numpy-arrays-as-lookup-tables
        if int(scen_years[0])<1765:
            raise ValueError("CFCs can only be infilled from RCPs as far "+
            "back as 1765 at present")
        rcp_years = np.arange(scen_years[0], scen_years[-1]+1)
        mapping = dict(zip(rcp_emis[:,0], range(rcp_emis.shape[0])))
        rcp_cfcs = np.array([rcp_emis[mapping[key],24:] for key in rcp_years])
        emissions_filled = np.append(emissions_filled, rcp_cfcs, axis=1)
    else:
        raise ValueError("include_cfcs should be an nt x 16 numpy array, a " +
        "string (rcp3pd, rcp45, rcp6 or rcp85) or False.")

    # Fill in any pre-SCEN years from RCP4.5. All pathways are identical <2000
    if scen_years[0]>startyear:
        if scen_years[0]>2000:
            raise ValueError("Can only fill in history unambiguously if "    +
            "first year in SCEN file is 2000 or earlier. You have requested "+
            "startyear=%d and the first year in SCEN file is %d"
            % (startyear, scen_years[0]))
        else:
            # tack RCP45 on to beginning
            rcp_emis = _import_emis_file('rcp45').emissions
            rcp_years = np.arange(startyear, scen_years[0])
            mapping = dict(zip(rcp_emis[:,0], range(rcp_emis.shape[0])))
            if include_cfcs==False:
                rcp_all = np.array([rcp_emis[mapping[key],:24] for key in rcp_years])
            else:
                rcp_all = np.array([rcp_emis[mapping[key],:] for key in rcp_years])
            emissions_filled = np.insert(emissions_filled, 0, rcp_all, axis=0)

        # harmonise?
        if harmonise is not None:
            harmonise = int(harmonise)
            if harmonise < 2000:
                raise ValueError("Cannot harmonise before 2000.")
            elif harmonise > scen_years[-1]:
                 raise ValueError("Cannot harmonise after last year of "      +
                 "input dataset")
            rcp_emis_2000 = rcp_emis[2000-startyear,:]
            for j in range(1,emissions_filled.shape[1]):
                f = interp1d((2000,harmonise), (rcp_emis_2000[j],
                    emissions_filled[harmonise-startyear,j]))
                emissions_filled[2000-startyear:harmonise-startyear,j] = f(
                    np.arange(2000,harmonise))

    return emissions_filled
