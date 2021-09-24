__all__ = [
    'VolcanoAllocator', 'msvolso2l4', 'rc2nc',
    'runcfg', 'runmsvolso2l4', 'runrc', 'version'
]


import pandas as pd
import warnings
import configparser
from ._core import VolcanoAllocator, msvolso2l4, rc2nc, __version__

__doc__ = """
Code for making CMAQ-ready emission files. Files can be produced from either
GEOS-Chem's *.rc files (rc2cn) or OMI's eruption database (msvolso2l4).

Outputs for degassing have been checked for consistency with Itahashi et al.[1]
and are consistent for mass, vertical allocation, and general characterisitics.

For eruptions, the vertical allocation follows the recommendation in the GEOS-
Chem rc file comments.

[1] Itahashi et al., Incorporation of volcanic SO2 emissions in the Hemispheric
 CMAQ (H-CMAQ) version 5.2 modeling system and assessing their impacts on
 sulfate aerosol over the Northern Hemisphere, GMD 2021.
 https://gmd.copernicus.org/articles/14/5751/2021/
"""

version = __version__


def runcfg(cfgobjs, cfgtype='path', warningfilter='ignore'):
    """
    Arguments
    ---------
    cfgobjs : list
        List of paths, files, or dictionaries
    cfgtype : str
        'path', 'file', or 'dict'
    warningfilter : str
        str accepted by warnings.simplefilter
    """
    warnings.simplefilter(warningfilter)
    config = configparser.ConfigParser(
        default_section='common',
        interpolation=configparser.ExtendedInterpolation()
    )
    config.read_dict({'common': {'overwrite': False, 'verbose': 0}})

    if cfgtype == 'path':
        config.read(cfgobjs)
    elif cfgtype == 'file':
        for cfgf in cfgobjs:
            config.read_file(cfgf)
    elif cfgtype == 'dict':
        for cfgd in cfgobjs:
            config.read_dict(cfgd)
    else:
        raise ValueError(
            'cfgtype must be either path, file or dict; got {cfgtype}'
        )

    overwrite = config.getboolean('common', 'overwrite')
    verbose = config.getint('common', 'verbose')
    start_date = config.get('common', 'start_date')
    end_date = config.get('common', 'end_date')

    m3dtmpl = config.get('mcip', 'm3dtmpl')
    g2dtmpl = config.get('mcip', 'g2dtmpl')

    dates = pd.date_range(start_date, end_date, freq='D')

    sections = config.sections()
    if 'msvolso2l4' in sections:
        l4inpath = config.get('msvolso2l4', 'inpath')
        l4outtmpl = config.get('msvolso2l4', 'outtmpl')

    if 'rc' in sections:
        rcintmpl = config.get('rc', 'intmpl')
        rcouttmpl = config.get('rc', 'outtmpl')

    if 'msvolso2l4' in sections:
        runmsvolso2l4(dates, l4inpath, l4outtmpl, m3dtmpl, g2dtmpl, overwrite, verbose)

    if 'rc' in sections:
        runrc(dates, rcintmpl, rcouttmpl, m3dtmpl, g2dtmpl, overwrite, verbose)


def runmsvolso2l4(dates, l4inpath, l4outtmpl, m3dtmpl, g2dtmpl, overwrite, verbose):
    print('Starting MSVOLSO2L4 processing')

    msvol = msvolso2l4(
        inpath=l4inpath, m3dtmpl=m3dtmpl, g2dtmpl=g2dtmpl, verbose=verbose
    )

    processed = 'None'
    nodata = 'None'
    existed = 'None'

    for date in dates:
        datestr = date.strftime('%F')
        outpath = msvol.to_netcdf(date, l4outtmpl, overwrite=overwrite)
        if outpath is None:
            existed = datestr
        elif outpath == 0:
            nodata = datestr
        else:
            processed = datestr
        print(
            f'\rLast {datestr}, done {processed}, no data {nodata} or exists {existed}',
            end='', flush=True
        )

    print('\nOMI L4 Done')


def runrc(dates, rcintmpl, rcouttmpl, m3dtmpl, g2dtmpl, overwrite, verbose):
    print('Starting rc processing')
    existed = 'None'
    processed = 'None'
    for date in dates:
        datestr = date.strftime('%F')
        outpath = rc2nc(
            rcpaths=[date.strftime(rcintmpl)],
            m3dtmpl=m3dtmpl, g2dtmpl=g2dtmpl
        ).to_netcdf(
            outtmpl=rcouttmpl, verbose=0, overwrite=overwrite
        )
        if outpath is None:
            existed = datestr
        else:
            processed = datestr

        print(
            f'\rLast done {processed}; last existed {existed}', end='', flush=True
        )

    print('\nDone')
