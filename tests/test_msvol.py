import os
import pandas as pd
import volcano2cmaq
import pyrsig
import pycno

outpath = 'outputs/MSVOLSO2L4_20260129_1188NHEMI2_202007.nc'
inpath = 'inputs/MSVOLSO2L4_20260129.txt'
g2dtmpl = 'inputs/mcip/GRIDCRO2D_1188NHEMI2'
m3dtmpl = 'inputs/mcip/METCRO3D_1188NHEMI2'

if os.path.exists(outpath):
    os.remove(outpath)

dates = pd.date_range('2020-06-30', '2020-08-01')

vp = volcano2cmaq.msvolso2l4(inpath, m3dtmpl, g2dtmpl, verbose=9)
vp.to_netcdf(dates, outpath)
f = pyrsig.open_ioapi(outpath)
cno = pycno.cno(f.crs_proj4)
os.makedirs('figs', exist_ok=True)
for key in ['SO2_ERUPT', 'SO2_DEGAS']:
    figstem = os.path.join('figs', os.path.basename(outpath))
    if key in f:
        qm = f[key].sum(('ROW', 'COL')).where(lambda x:x>0).T.plot(figsize=(8, 4))
        qm.axes.set(ylim=(1, 0), title=outpath)
        qm.figure.savefig(f'{figstem}_{key}_timebylayer.png')
        qm = f[key].sum(('TSTEP', 'LAY')).where(lambda x:x>0).plot(figsize=(6, 4))
        cno.draw(ax=qm.axes)
        qm.figure.savefig(f'{figstem}_{key}_map.png')
