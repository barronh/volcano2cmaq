import configparser
import argparse
import pandas as pd
import os
import warnings
import volcano2ioapi

warnings.simplefilter('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('config', nargs='*', default=['run.cfg'])

args = parser.parse_args()

config = configparser.ConfigParser(
    default_section='common',
    interpolation=configparser.ExtendedInterpolation()
)
config.read(args.config)

overwrite = config.getboolean('common', 'overwrite')
dom = config.get('common', 'dom')
start_date = config.get('common', 'start_date')
end_date = config.get('common', 'end_date')

m3dtmpl = config.get('mcip', 'm3dtmpl')
g2dtmpl = config.get('mcip', 'g2dtmpl')

dates = pd.date_range(start_date, end_date, freq='D')

l4inpath = config.get('msvolso2l4', 'inpath')
l4outtmpl = config.get('msvolso2l4', 'outtmpl')

rcintmpl = config.get('rc', 'intmpl')
rcouttmpl = config.get('rc', 'outtmpl')

msvol = volcano2ioapi.msvolso2l4(
    inpath=l4inpath, m3dtmpl=m3dtmpl, g2dtmpl=g2dtmpl, verbose=0
)

processed = 'None'
skipped = 'None'
existed = 'None'

for date in dates:
    datestr = date.strftime('%F')
    outpath = date.strftime(l4outtmpl)
    outf = msvol.allocate(date, outpath, overwrite=overwrite)
    if outf is not None:
        processed = datestr
    elif os.path.exists(outpath):
        skipped = datestr
    else:
        existed = datestr
    print(
        f'\rLast done {processed}, skipped {skipped}, existed {existed}',
        end='', flush=True
    )
print('\nDone')

existed = 'None'

for date in dates:
    outpath = date.strftime(rcouttmpl)
    print(
        f'\rWorking on {outpath}; last skipped {existed}', end='', flush=True
    )
    volcano2ioapi.rc2nc(
        rcpaths=[date.strftime(rcintmpl)], m3dtmpl=m3dtmpl, g2dtmpl=g2dtmpl
    ).allocate(
        outpath=outpath, verbose=0, overwrite=overwrite
    )

print('\nDone')
