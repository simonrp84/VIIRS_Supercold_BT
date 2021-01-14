# !/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2020 Simon R Proud
#
# This file is part of the code for Supercold_BT.
#
# Supercold_BT is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Supercold_BT is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# Supercold_BT.  If not, see <http://www.gnu.org/licenses/>.

"""Find minimum brightness temperatures in a scene and save to file."""
from datetime import datetime
from os.path import exists
from satpy import Scene
from glob import glob

import numpy as np
import warnings
import pathlib
import netCDF4
import sys

warnings.filterwarnings("ignore")


def make_ncdf(out_f, nc_data, arr_size):
    fid = netCDF4.Dataset(out_f, 'w')
    fid.createDimension('lat', 180 * arr_size)
    fid.createDimension('lon', 360 * arr_size)
    later = fid.createVariable('lat', np.float32, ('lat',))
    loner = fid.createVariable('lon', np.float32, ('lon',))
    spf = fid.createVariable('AQUA_MINTEMP',
                             np.float32,
                             ('lat', 'lon',),
                             fill_value=9999.)

    later[:] = np.arange(-90, 90, 1.0 / sizer)
    loner[:] = np.arange(-180, 180, 1.0 / sizer)

    spf[:] = nc_data[:, :]

    fid.close()


def compute_map(in_arr, in_lats, in_lons, data_t, in_size):
    """Compute per-pixel mintemps"""

    in_lats = np.array(np.round((in_lats + 90) * in_size), dtype=np.int)
    in_lons = np.array(np.round((in_lons + 180) * in_size), dtype=np.int)

    gd_pts = (in_lats >= 0).nonzero()
    in_lats = in_lats[gd_pts]
    in_lons = in_lons[gd_pts]
    data_t = data_t[gd_pts]
    gd_pts = (in_lats < sizer*180).nonzero()
    in_lats = in_lats[gd_pts]
    in_lons = in_lons[gd_pts]
    data_t = data_t[gd_pts]
    gd_pts = (in_lons > 0).nonzero()
    in_lats = in_lats[gd_pts]
    in_lons = in_lons[gd_pts]
    data_t = data_t[gd_pts]
    gd_pts = (in_lons < sizer*360).nonzero()
    in_lats = in_lats[gd_pts]
    in_lons = in_lons[gd_pts]
    data_t = data_t[gd_pts]

    for i in range(0, len(data_t)):
        if data_t[i] < in_arr[in_lats[i], in_lons[i]]:
            in_arr[in_lats[i], in_lons[i]] = data_t[i]

    return in_arr


def preproc_data(infile):
    """Read MODIS data and find good pixels."""
    # Load relevant bands
    scn = Scene([infile], reader='modis_l1b')
    scn.load(['31'], resolution=1000)

    # Get lat / lon grid to ensure we focus on tropics
    in_lons, in_lats = scn['31'].area.get_lonlats()
    in_lats = in_lats.compute()
    in_lons = in_lons.compute()

    in_lats = in_lats.ravel()
    in_lons = in_lons.ravel()

    lats2 = np.abs(in_lats)

    if np.nanmin(lats2) < 33.5:
        scn_data = scn['31'].values
        scn_data = scn_data.ravel()
        # Only look at points in tropics
        gd_pts = (lats2 < 33.5).nonzero()
        scn_data = scn_data[gd_pts]
        in_lons = in_lons[gd_pts]
        in_lats = in_lats[gd_pts]
        lats2 = lats2[gd_pts]

        # Sometimes there's bad data, <155K BT
        gd_pts = (scn_data > 155).nonzero()
        scn_data = scn_data[gd_pts]
        in_lons = in_lons[gd_pts]
        in_lats = in_lats[gd_pts]
        lats2 = lats2[gd_pts]

        # Now find points colder than 205K for stats
        gd_pts = (scn_data < 205).nonzero()
        scn_data = scn_data[gd_pts]
        in_lons = in_lons[gd_pts]
        in_lats = in_lats[gd_pts]
        lats2 = lats2[gd_pts]
        if len(gd_pts[0]) < 1:
            return None, None, None
        if np.nanmin(lats2) > 33.5:
            return None, None, None
    else:
        # Every everything is high latitude, return nothing
        print(" -   All data-points outside latitude range, skipping.")
        return None, None, None
    # Return data
    return in_lats, in_lons, scn_data


# Define satellite we're processing
satproc = 'MYD'

# Define input directories, prime and backup (CEMS has lots of missing data)
indir_m2 = '/neodc/modis/data/' + satproc + '021KM/collection61/'
indir_m2_bk = '/gws/nopw/j04/nerc_avmet/MODIS/MOD_REPL' + satproc + '021KM/'


# Define output directory (personal GWS as above)
odir = '/gws/nopw/j04/nerc_avmet/MODIS/MINT/'

# Find date to process from command line
dtstr = sys.argv[1]
dater = datetime.strptime(dtstr, "%Y%m%d")
yrstr = dater.strftime("%Y")
mnstr = dater.strftime("%m")
dystr = dater.strftime("%d")
doystr = dater.strftime("%j")


# Set up output directory
# This is for extreme values in MODIS
odir_ex = odir + '/EXTREME/' + dater.strftime("%Y/%m") + '/'
# This is for daily stats
odir_dl = odir + '/DAILY/' + dater.strftime("%Y/%m") + '/'
# This is for global map of coldest per-pix temperatures
odir_ar = odir + '/MAP/' + dater.strftime("%Y/%m") + '/'

# Make directories
pathlib.Path(odir_ex).mkdir(parents=True, exist_ok=True)
pathlib.Path(odir_dl).mkdir(parents=True, exist_ok=True)
pathlib.Path(odir_ar).mkdir(parents=True, exist_ok=True)

# Find L1b radiances, have to check both directories due to missing
# data on CEDA
m2_files = glob(indir_m2 + yrstr + '/' + mnstr +
                '/' + dystr + '/*' + yrstr + doystr + '*.hdf')
m2_files_bk = glob(indir_m2_bk + yrstr + '/' + mnstr +
                   '/' + dystr + '/*' + yrstr + doystr + '*.hdf')

m2_files = m2_files + m2_files_bk
m2_files.sort()

# Output filenames
outf_ex = odir_ex + satproc + '_EX_' + yrstr + mnstr + dystr + '.csv'
outf_dl = odir_dl + satproc + '_DAILY_' + yrstr + mnstr + dystr + '.csv'
outf_ar = odir_ar + satproc + '_MINT_' + yrstr + mnstr + dystr + '.nc'

# Set up array for geographic mintemp
sizer = 4
mainarr = np.zeros((180*sizer, 360*sizer))
mainarr[:, :] = 9999.

if exists(outf_ex) and exists(outf_dl) and exists(outf_ar):
    print("Already processed:", yrstr + mnstr + dystr)
    quit()

odata = np.array([])

n_files = len(m2_files)
count = 1

t_thr = 178.15

fid_ex = open(outf_ex, 'w')
fid_dl = open(outf_dl, 'w')

# Now loop over all files for this day
for inf in m2_files:
    print(inf)
    try:
        # Load data
        lats, lons, data = preproc_data(inf)
        if data is None:
            continue
        mainarr = compute_map(mainarr, lats, lons, np.copy(moddata), sizer)
        # Find if any points are in the extreme cold category
        if np.nanmin(data) < t_thr:
            pts = (data < t_thr).nonzero()
            # If there are, write to file
            for pt in pts[0]:
                outstr = inf + ',' + str(lats[pt]) + ','
                outstr = outstr + str(lons[pt]) + ',' + str(data[pt]) + '\n'
                fid_ex.write(outstr)
        # Only append if there's data
        if data.shape[0] > 0:
            odata = np.concatenate([odata, data])
    except Exception as e:
        print("Problem with: " + inf + '. Error:', e)
        pass
    count += 1

# Get lists of points for BTs below various thresholds
pts165 = (odata < 165).nonzero()
pts170 = (odata < 170).nonzero()
pts175 = (odata < 175).nonzero()
pts180 = (odata < 180).nonzero()
pts185 = (odata < 185).nonzero()
pts190 = (odata < 190).nonzero()
pts195 = (odata < 195).nonzero()

# Add some stats to the output string
outstr = yrstr + '-' + mnstr + '-' + dystr + ','
outstr = outstr + str(np.nanmin(odata)) + ','
outstr = outstr + str(np.nanpercentile(odata, 1)) + ','
outstr = outstr + str(np.nanpercentile(odata, 2)) + ','
outstr = outstr + str(np.nanpercentile(odata, 5)) + ','
outstr = outstr + str(np.nanpercentile(odata, 10)) + ','
outstr = outstr + str(np.nanpercentile(odata, 20)) + ','

# This tedious section adds, if they exist, data about any particularly cold BTs.
# Otherwise, it skips that temperature.
try:
    outstr = outstr + str(len(pts165[0])) + ','
except:
    outstr = outstr + ','
try:
    outstr = outstr + str(len(pts170[0])) + ','
except:
    outstr = outstr + ','
try:
    outstr = outstr + str(len(pts175[0])) + ','
except:
    outstr = outstr + ','
try:
    outstr = outstr + str(len(pts180[0])) + ','
except:
    outstr = outstr + ','
try:
    outstr = outstr + str(len(pts185[0])) + ','
except:
    outstr = outstr + ','
try:
    outstr = outstr + str(len(pts190[0])) + ','
except:
    outstr = outstr + ','
try:
    outstr = outstr + str(len(pts195[0])) + ','
except:
    outstr = outstr + ','

print(outstr)
outstr = outstr + '\n'

pts = odata

make_ncdf(outf_ar, mainarr, sizer)

fid_dl.write(outstr)
fid_dl.close()
