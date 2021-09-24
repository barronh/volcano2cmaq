import os
import PseudoNetCDF as pnc
import pandas as pd
import numpy as np
from datetime import datetime
from warnings import warn

__version__ = '0.1.1'
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


class VolcanoAllocator:
    def __init__(self, m3dtmpl, g2dtmpl, cache=True, verbose=0):
        """
        Arguments
        ---------
        m3dtmpl : str
            strftime template to construct a path by day for METCRO3D file or
            any file that has ZF or ZH in meters agl
        g2dtmpl : str
            strftime template to construct a path by day for GRIDCRO2D file or
            any file that has HT in meters asl.
        cache : bool
            If True, store met files for reuse. This is useful when met files
            are climatology.
        """
        self._zfs = {}
        self._cache = cache
        self._cache = cache
        self._g2dtmpl = g2dtmpl
        self._m3dtmpl = m3dtmpl
        self.verbose = verbose

    def _prepoutfile(self, dateobj, erupt=True, degas=True):
        keys = []
        if erupt:
            keys.append('SO2_ERUPT')
        if degas:
            keys.append('SO2_DEGAS')
        rawzf = pnc.pncopen(dateobj.strftime(self._m3dtmpl), format='ioapi')
        if 'ZF' in rawzf.variables:
            zkey = 'ZF'
        else:
            zkey = 'ZH'
        outf = rawzf.subset([zkey]).apply(TSTEP='mean').copy(
        ).renameVariables(**{zkey: keys[0]})
        outf.SDATE = int(dateobj.strftime('%Y%j'))
        outf.STIME = 0
        outf.TSTEP = 240000
        outv = outf.variables['SO2_ERUPT']
        outv[:] = 0
        outv.units = 'moles/s'
        outv.long_name = keys[0].ljust(16)
        var_descs = {
            'SO2_ERUPT': 'Sulfur as SO2 from eruptions',
            'SO2_DEGAS': 'Sulfur as SO2 from degassing'
        }
        outv.var_desc = var_descs[keys[0]].ljust(80)
        if len(keys) > 1:
            outv = outf.copyVariable(outv, key=keys[1])
            outv.long_name = keys[1].ljust(16)
            outv.var_desc = var_descs[keys[1]].ljust(80)

        if isinstance(self, rc2nc):
            filedesc = f"""
AeroCom volcanic emissions of SO2 for 1978-2019 obtained from Christoph Keller
(NASA/GMAO). Translated from GEOS-Chem rc files to IOAPI netCDF files by Barron
H. Henderson for use in CMAQ.

Consider speciation. Perhaps assume a ratio of SO2 to PM and then use profile
433012.5 and 4330110 to speciate PM.

volcano2ioapi.rc2nc v{__version__}
"""
        elif isinstance(self, msvolso2l4):
            filedesc = f"""
Eruption sulfur mass[1] emitted into top 1/3 of plume cloud height using met
files[2,3]. Original mass in kt(so2). Converted to g/s (kt * 1e9 / 24 / 3600)
and then converted to moles by divideing by 64.066.

volcano2ioapi.msvolso2l4 v{__version__}

[1] https://so2.gsfc.nasa.gov/eruptions/MSVOLSO2L4_20210621.txt
[2] {self._g2dtmpl}
[3] {self._m3dtmpl}
"""
        else:
            warn('Type unknown')
            filedesc = ''
        outf.FILEDESC = filedesc
        outf.HISTORY = ''
        delattr(outf, 'VAR-LIST')
        outf.updatemeta(sortmeta=True)
        outf.CDATE = outf.WDATE
        outf.CTIME = outf.WTIME
        outf.updatetflag(overwrite=True)
        outf.updatemeta(sortmeta=True)
        outf.dimensions.move_to_end('TSTEP', last=False)
        return outf

    # def ll2ij(self, lon, lat):
    #     gf = None
    #     return gf.ll2ij(lon, lat)

    # def get_zedges(self, i, j):
    #     pass

    def layerfrac(self, dateobj, i, j, elevation, cloud_column_height):
        """
        Arguments
        ---------
        dateobj : datetime.datetime
            date must implement strftime
        elevation : float
            meters above sea level of the volcano surface
        cloud_column_height : float
            meters above sea level of the plume height
        """
        verbose = self.verbose
        zfile = self.get_zedgefile(dateobj)

        yp = [0, 1]
        if 'ZF_ASL' in zfile.variables:
            ze = zfile.variables['ZF_ASL'][:].mean(0)[:, j, i]
            heights = ze[:].copy()
            # If the surface height is below the volcano, reset minimum
            # This should be rare for coarse domains.
            heights[0] = np.minimum(heights[0], elevation)
            if elevation == cloud_column_height:
                # Allocate to surface layer
                xp = np.array([
                    elevation, cloud_column_height + 0.001
                ])
            else:
                xp = np.array([
                    elevation, cloud_column_height
                ])

            layercdf = np.interp(heights, xp, yp, left=0, right=1)
            layerfrac = np.diff(layercdf)
        elif 'ZH_ASL' in zfile.variables:
            xp = ze = heights = zfile.variables['ZH_ASL'][0, 1:, j, i]
            yp = np.arange(heights.size)
            if elevation == cloud_column_height:
                k = int(np.interp(elevation, heights, yp, left=0).round(0))
                layerfrac = np.zeros_like(heights)
                layerfrac[k] = 1
            else:
                bottom = (
                    cloud_column_height - (cloud_column_height - elevation) / 3
                )
                top = cloud_column_height
                intopthird = (heights > bottom) & (heights < top)
                layerfrac = intopthird / intopthird.sum()
                assert(layerfrac.sum().round(5) == 1)
        else:
            raise KeyError('Must have either ZF_ASL or ZH_ASL')

        if verbose > 1:
            print(xp, yp)
            print(ze[:])
            print(heights)
            print(layerfrac)

        return layerfrac

    def get_gridfile(self, dateobj):
        """
        Returns a file with the ll2ij interface

        Arguments
        ---------
        dateobj : datetime.datetime
            date must implement strftime
        """
        g2f = pnc.pncopen(
            dateobj.strftime(self._g2dtmpl), format='ioapi'
        )
        return g2f

    def get_zedgefile(self, dateobj):
        """
        From a date and file strftime templates, create a file with ZF and ZH
        in meters above sea level (ASL) with a first layer height

        Arguments
        ---------
        dateobj : datetime.datetime
            date must implement strftime
        """
        if dateobj in self._zfs:
            return self._zfs[dateobj]
        if self.verbose > 0:
            print(f'Opening {dateobj.strftime("%F")}')
        m3f = pnc.pncopen(
            dateobj.strftime(self._m3dtmpl), format='ioapi'
        )
        keepkeys = []
        for zkey in ['ZF', 'ZH']:
            if zkey in m3f.variables:
                keepkeys.append(zkey)
        kidx = np.append(0, np.arange(m3f.NLAYS))
        zfile = m3f.subset(keepkeys).apply(TSTEP='mean')
        zefile = zfile.slice(LAY=kidx)
        g2f = pnc.pncopen(
            dateobj.strftime(self._g2dtmpl), format='ioapi'
        )

        zfile.copyVariable(g2f.variables['HT'], key='HT')
        if 'ZH' in keepkeys:
            zvar = zefile.copyVariable(
                zfile.variables['ZH'], key='ZH_ASL', withdata=False
            )
            zvar[:] = np.concatenate(
                [
                    zfile.variables['HT'][:, [0]],
                    zfile.variables['ZH'][:] + zfile.variables['HT'][:]
                ], axis=1
            )

        if 'ZF' in keepkeys:
            zvar = zefile.copyVariable(
                zfile.variables['ZF'], key='ZF_ASL', withdata=False
            )
            zvar[:] = np.concatenate(
                [
                    zfile.variables['HT'][:, [0]],
                    zfile.variables['ZF'][:] + zfile.variables['HT'][:]
                ], axis=1
            )
        if self._cache:
            self._zfs[dateobj] = zefile

        return zefile


class msvolso2l4(VolcanoAllocator):
    def __init__(
        self, inpath, cache=True, m3dtmpl=None, g2dtmpl=None, verbose=0
    ):
        """
        OMI Volcano eruption text file processor

        Arguments
        ---------
        inpath : str
            Tab delimited file with the following columns volcano, lat, lon,
            v_alt, yyyy, mm, dd, type, vei, p_alt_obs, p_alt_est, and so2(kt).
            Text file is available from https://so2.gsfc.nasa.gov/measures.html
            e.g., https://so2.gsfc.nasa.gov/eruptions/MSVOLSO2L4_20210621.txt
            date likely to change

        Notes
        -----
        see VolcanoAllocator for other keywords
        """
        # Open data and prep meta variables
        data = self._data = pd.read_csv(
            inpath, delimiter=r'\s+', usecols=list(range(12)),
            na_values=[-999]
        )
        self._data['DateObj'] = pd.to_datetime(
            (data.yyyy * 10000 + data.mm * 100 + data.dd).astype('str')
        )
        self._data['DOY'] = data['DateObj'].dt.dayofyear

        # Initialize met processing capability
        VolcanoAllocator.__init__(
            self, g2dtmpl=g2dtmpl, m3dtmpl=m3dtmpl, cache=True, verbose=verbose
        )

    def allocate(self, dateobj, outpath=None, overwrite=False, **save_kw):
        """
        Allocate volcanic emissions from date to a IOAPI file.

        Arguments
        ---------
        dateobj : datetime.datetime
            Date to allocate
        """
        verbose = self.verbose
        if outpath is not None:
            if os.path.exists(outpath) and not overwrite:
                if verbose > 0:
                    print(f'{outpath} exists')
                return None

        if verbose > 0:
            print('Querying data')

        # Query data for same day data
        alldata = self._data.query(
            f'DateObj == "{dateobj.strftime("%Y-%m-%d")}"'
        ).copy()

        if alldata.shape[0] == 0:
            if verbose > 0:
                print(f'No data found for {dateobj.strftime("%Y-%m-%d")}')
            return None

        if verbose > 0:
            print('Calculating plume top and I/J')

        # Set p_alt to prefer observations, but settle for estimation
        alldata['p_alt'] = alldata['p_alt_obs']
        alldata.loc[alldata['p_alt'].isna(), 'p_alt'] = (
            alldata.loc[alldata['p_alt'].isna(), 'p_alt_est']
        )

        gfile = self.get_gridfile(dateobj)
        alldata['I'], alldata['J'] = gfile.ll2ij(
            alldata['lon'], alldata['lat']
        )

        # Subset data for within the domain
        data = alldata.query(
            f'I >= 0 and I < {gfile.NCOLS} and J >= 0 and J < {gfile.NROWS}'
        )
        if data.shape[0] == 0:
            if verbose > 0:
                print(
                    'No data found in domain for '
                    + dateobj.strftime("%Y-%m-%d")
                )
            return None

        if verbose > 0:
            print('Grouping and aggregating')

        gbdata = data.groupby(
            ['I', 'J', 'type', 'v_alt', 'p_alt'], as_index=False
        )
        adata = gbdata.agg(SO2=('so2(kt)', 'sum'))

        if verbose > 0:
            print('Calculating plume range')

        adata['stop'] = adata['p_alt'] * 1000
        # Calculate top 1/3 bottom
        adata['start'] = (
            adata['stop'] - (adata['stop'] - adata['v_alt'] * 1000) / 3
        )

        if verbose > 0:
            print('Preparing output file')

        outf = self._prepoutfile(dateobj, degas=False)
        outv = outf.variables['SO2_ERUPT']

        if verbose > 0:
            print('Allocating data to grid and layers')

        for idx, row in adata.iterrows():
            j, i = row['J'], row['I']
            so2molerate = row['SO2'] * 1000000000 / (24 * 3600) / 64.066
            layerfrac = self.layerfrac(
                dateobj, i, j, adata['start'].values[0],
                adata['stop'].values[0]
            )

            outvals = layerfrac * so2molerate
            if self.verbose > 0:
                print(outvals)
            outv[0, :, j, i] += outvals
        if outpath is not None:
            save_kw.setdefault('complevel', 1)
            diskf = outf.save(outpath, verbose=verbose, **save_kw)
            return diskf
        else:
            return outf


class rc2nc(VolcanoAllocator):
    def __init__(
        self, rcpaths, m3dtmpl, g2dtmpl, process=True, verbose=0,
    ):
        """
        rc2nc reads text files from GEOS-Chem and HEMCO's volcano extension
        and can create IOAPI-like input files

        Arguments
        ---------
        rcpaths : list
            One or more paths. See addrcfile for more details.
        outpath : str
            If None, defaults to output/Volcano.STARTDATE-ENDDATE.nc
        process : bool
            If true, run default process
        verbose : int
            increasing levels of verbosity

        Notes
        -----
        see VolcanoAllocator for other keywords
        """
        self._rcpaths = rcpaths
        sdatestr = rcpaths[0].split('.')[-2]
        self.sdateobj = datetime.strptime(sdatestr, '%Y%m%d')
        VolcanoAllocator.__init__(
            self, cache=True, m3dtmpl=m3dtmpl, g2dtmpl=g2dtmpl, verbose=verbose
        )

    def addrcfile(self, dateobj, rcpath, efile):
        """
        Add mass from rcpath to emission file (efile) using the heights
        from the height file (zfile).

        Arguments
        ---------
        dateobj : datetime.datetime
            Date to allocate
        rcpath : str
            Path to a text file with comments preceded by # followed
            by the world "volcano::". All subsequent lines are space delimited
            with columns lat (-90, 90), lon (-180, 180), massrate (kg S/s),
            elevation (m), and cloud_column_height (m). These files are
            GEOS-Chem inputs
            https://ftp.as.harvard.edu/gcgrid/data/ExtData/HEMCO/VOLCANO/
        efile : netcdf-like
            File for output that has IOAPI properties, implements ll2ij, and
            has SO2_ERUPT and SO2_DEGAS variables
        """
        verbose = self.verbose
        # lat_degrees lon_degrees sulfure_kgSps elevation_m
        # cloud_column_height_m
        vdata = pd.read_csv(
            rcpath, delimiter=r'\s+', skiprows=4, comment=':',
            names=['lat', 'lon', 'massrate', 'elev', 'cch']
        )
        i, j = efile.ll2ij(vdata.lon.values, vdata.lat.values)
        vdata['I'] = i
        vdata['J'] = j
        if verbose > 0:
            print('Total Sulfur')
            print(rcpath, vdata.massrate.sum())
            print('All degas (True); All erupt (False); Both (True, False)')
            print(np.unique(vdata.elev == vdata.cch))
        udata = vdata.query(
          f'(I >= 0) and (I < {efile.NCOLS})'
          + f' and (J >= 0) and (J < {efile.NROWS})'
        )
        ijecs = udata.groupby(['I', 'J', 'elev', 'cch']).sum()
        levels = np.zeros((44,), dtype='f')
        for (i, j, elevation, cloud_column_height), row in ijecs.iterrows():
            levels[:] = 0
            # Raw unit is kgS/s
            # Converting to moles/s
            molerate = row['massrate'] * 1000 / 32

            layerfrac = self.layerfrac(
                dateobj, i, j, elevation, cloud_column_height
            )

            if elevation == cloud_column_height:
                outv = efile.variables['SO2_DEGAS']
            else:
                outv = efile.variables['SO2_ERUPT']

            outv[0, :, j, i] += molerate * layerfrac

    def allocate(self, verbose=0, outpath=None, overwrite=False, **save_kw):
        """
        Use all rcpaths supplied at intialization to create an output.

        Arguments
        ---------
        verbose : int
            increasing levels of verbosity
        """
        verbose = self.verbose
        if outpath is not None:
            if os.path.exists(outpath) and not overwrite:
                if verbose > 0:
                    print(f'{outpath} exists')
                return None

        rcpaths = self._rcpaths
        sdatestr = rcpaths[0].split('.')[-2]
        sdateobj = datetime.strptime(sdatestr, '%Y%m%d')
        efile = self._prepoutfile(sdateobj)
        efiles = [efile.copy() for i in rcpaths]
        # edatestr = rcpaths[-1].split('.')[-2]
        for efile, rcpath in zip(efiles, rcpaths):
            datestr = rcpath.split('.')[-2]
            dateobj = datetime.strptime(datestr, '%Y%m%d')
            efile.SDATE = int(dateobj.strftime('%Y%j'))
            efile.STIME = 0
            efile.TSTEP = 240000
            self.addrcfile(dateobj, rcpath, efile)

        sdate = datetime.strptime(sdatestr, '%Y%m%d')
        # edate = datetime.strptime(edatestr, '%Y%m%d')
        outf = efiles[0].stack(efiles[1:], 'TSTEP')
        outf.SDATE = int(sdate.strftime('%Y%j'))
        outf.STIME = 0
        outf.TSTEP = 240000
        outf.updatetflag(overwrite=True)
        outf.variables.move_to_end('TFLAG', last=False)
        outf.dimensions.move_to_end('TSTEP', last=False)
        if outpath is not None:
            save_kw.setdefault('complevel', 1)
            diskf = outf.save(outpath, verbose=verbose, **save_kw)
            return diskf
        else:
            return outf
