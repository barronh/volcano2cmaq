import os
import pandas as pd
import volcano2cmaq

outpath = 'outputs/so2_volcanic_emissions_Carns_1188NHEMI2_202007.nc'
g2dtmpl = 'inputs/mcip/GRIDCRO2D_1188NHEMI2'
m3dtmpl = 'inputs/mcip/METCRO3D_1188NHEMI2'
if os.path.exists(outpath):
    os.remove(outpath)
dates = pd.date_range('2020-06-30', '2020-08-01')
rcpaths = dates.strftime(f'HEMCO/VOLCANO/v2021-09/%Y/%m/so2_volcanic_emissions_Carns.%Y%m%d.rc')

vp = volcano2cmaq.rc2nc(rcpaths, m3dtmpl, g2dtmpl, verbose=9)
vp.to_netcdf(outpath)