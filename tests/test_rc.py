import os
import pandas as pd
import volcano2cmaq
import pyrsig
import pycno

for sfx in ['', '_ZHONLY']:
    outpath = f'outputs/so2_volcanic_emissions_Carns_1188NHEMI2_202007{sfx}.nc'
    g2dtmpl = 'inputs/mcip/GRIDCRO2D_1188NHEMI2'
    m3dtmpl = f'inputs/mcip/METCRO3D_1188NHEMI2{sfx}'
    if os.path.exists(outpath):
        os.remove(outpath)
    dates = pd.date_range('2020-06-30', '2020-08-01')
    rcpaths = dates.strftime(f'inputs/HEMCO/VOLCANO/v2024-04/%Y/%m/so2_volcanic_emissions_Carns.%Y%m%d.rc')
    vp = volcano2cmaq.rc2nc(rcpaths, m3dtmpl, g2dtmpl, verbose=9)
    vp.to_netcdf(outpath)
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
            