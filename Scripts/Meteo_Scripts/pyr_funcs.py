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

"""Helper functions for computing the data required for pyrttov."""

from .consts import get_avec, get_bvec
from netCDF4 import Dataset
import numpy as np
try:
    import pyrttov
except:
    pyrttov = None

def expand2nprofiles(n, nprof):
    """Transform 1D array to a [nprof, nlevels] array."""
    outp = np.empty((nprof, len(n)), dtype=n.dtype)
    for i in range(nprof):
        outp[i, :] = n[:]
    return outp


def get_sfc(inf_s, latp, lonp):
    """Get surface variables from ecmwf."""

    fid = Dataset(inf_s, 'r')

    # First get the lats + lons to compute closest point in array
    lats = fid['latitude'][:]
    lons = fid['longitude'][:]
    dist_la = np.abs(latp - lats)
    dist_lo = np.abs(lonp - lons)
    la_pt = (np.nanmin(dist_la) == dist_la).nonzero()[0]
    lo_pt = (np.nanmin(dist_lo) == dist_lo).nonzero()[0]

    # Now use this point to retrieve other values
    skint = np.squeeze(fid['skt'][:])[la_pt, lo_pt][0]
    sp = np.squeeze(fid['sp'][:])[la_pt, lo_pt][0]
    u10 = np.squeeze(fid['u10'][:])[la_pt, lo_pt][0]
    v10 = np.squeeze(fid['v10'][:])[la_pt, lo_pt][0]
    fid.close()
    
    return la_pt, lo_pt, skint, sp, u10, v10
  
    
def get_ml(inf_m, la_pt, lo_pt):
    """Get model level fields."""
    
    fid = Dataset(inf_m, 'r')
    spechum = np.squeeze(np.squeeze(fid['q'][:])[:, la_pt, lo_pt])
    temp = np.squeeze(np.squeeze(fid['t'][:])[:, la_pt, lo_pt])
    geop_m = np.squeeze(np.squeeze(fid['z'][:])[:, la_pt, lo_pt])
    ozone = np.squeeze(np.squeeze(fid['o3'][:])[:, la_pt, lo_pt])
    fid.close()

    # Number of levels is defined by the ECMWF file, usually 137
    nlevs = len(ozone)
    
    return spechum, temp, geop_m, ozone, nlevs


def correct_angs(satangs):
    """Ensure azimuths in right range."""
    if satangs[1] < 0:
        satangs[1] = 360 + satangs[1]
    if satangs[3] < 0:
        satangs[3] = 360 + satangs[3]
    return satangs


def get_pres(sp, nlevs):
    pres = np.zeros((nlevs))
    avec = get_avec()
    bvec = get_bvec()
    for i in range(1, nlevs+1):
        pres[i-1] = (avec[i] + bvec[i] * sp) / 100.

    return pres

def setup_atlas(ir_atlas_dir, brdf_atlas_dir, mn, inst_inst, nchan, nprofiles):
    """Set up the emissivity and BRDF atlases."""
    irAtlas = pyrttov.Atlas()
    irAtlas.AtlasPath = ir_atlas_dir
    irAtlas.loadIrEmisAtlas(mn, inst_inst, ang_corr=True)

    brdfAtlas = pyrttov.Atlas()
    brdfAtlas.AtlasPath = brdf_atlas_dir
    brdfAtlas.loadBrdfAtlas(mn, inst_inst) # Supply RTTOV object to enable single-instrument initialisation
    brdfAtlas.IncSea = False  
    surfemisrefl = np.zeros((4, nprofiles, nchan), dtype=np.float64)
    surfemisrefl[:,:,:] = -1.
    surfemisrefl[0,:,:] = irAtlas.getEmisBrdf(inst_inst)
    surfemisrefl[1,:,:] = brdfAtlas.getEmisBrdf(inst_inst)
    return surfemisrefl
    
