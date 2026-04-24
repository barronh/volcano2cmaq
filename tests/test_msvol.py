import os
import pandas as pd
import volcano2cmaq

outpath = 'outputs/MSVOLSO2L4_20260129_1188NHEMI2_202007.nc'
inpath = 'inputs/MSVOLSO2L4_20260129.txt'
g2dtmpl = 'inputs/mcip/GRIDCRO2D_1188NHEMI2'
m3dtmpl = 'inputs/mcip/METCRO3D_1188NHEMI2'

if os.path.exists(outpath):
    os.remove(outpath)

dates = pd.date_range('2020-06-30', '2020-08-01')

vp = volcano2cmaq.msvolso2l4(inpath, m3dtmpl, g2dtmpl, verbose=9)
vp.to_netcdf(dates, outpath)