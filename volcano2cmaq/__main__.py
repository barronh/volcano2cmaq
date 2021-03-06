import argparse

if __name__ == '__main__':
    from . import runcfg
    parser = argparse.ArgumentParser(
        prog='python -m volcano2cmaq',
        description=(
            'Make volcano emission files in IOAPI-like format for CMAQ.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        'config', nargs='*', default=['run.cfg'],
        help=(
            'Configuraton file or files. Parsed in order according to the'
            + ' configparser approach using extended interpolation'
        )
    )
    parser.epilog = """
Example:
    $ cat run.cfg
    [common]
    # Start and end dates (inclusive) to make daily files
    start_date = 2018-08-01
    end_date = 2018-08-31

    # Used in other sections, but not by the program
    dom = 108NHEMI2

    # Overwrite existing outputs; defaults to no
    # overwrite = yes
    # How much information to provide; defaults to 0
    # verbose = 3

    [mcip]
    # Must have ZF or ZH
    m3dtmpl = mcip/${dom}/METCRO3D.${dom}.44L.16%m
    # Must have HT
    g2dtmpl = mcip/${dom}/GRIDCRO2D.${dom}.44L.16%m

    [msvolso2l4]
    # This section provides configuration for the OMI MSVOLSO2 L4 processor
    # Downloaded from https://so2.gsfc.nasa.gov/measures.html
    inpath = input/MSVOLSO2L4_20210621.txt
    # Creates output where % are replaced with date components from strftime
    outtmpl = output/${dom}/MSVOLSO2L4_%Y-%m-%d.${dom}.nc

    [rc]
    # This section provides configuration for the HEMCO rc file processor
    # Path convenience variables; not used in program
    rcroot = ftp.as.harvard.edu/gcgrid/data/ExtData/HEMCO/VOLCANO/v2019-08

    # Inpaths where % are replaced with date components from strftime
    intmpl = ${rcroot}/%Y/%m/so2_volcanic_emissions_Carns.%Y%m%d.rc
    # Outpaths where % are replaced with date components from strftime
    outtmpl = output/${dom}/so2_volcanic_emissions_Carns.%Y-%m-%d.${dom}.nc

    $ python -m volcano2cmaq run.cfg
"""

    args = parser.parse_args()
    runcfg(args.config)
