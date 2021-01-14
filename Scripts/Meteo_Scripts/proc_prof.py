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

from .pyr_funcs import *
import numpy as np
try:
    import pyrttov
except:
    pyrttov = None


def make_prof_temp(trop_temp, trop_alt, ot_temp, lapse_rate):
    """Find a convective overshoot altitude."""
    del_h = (ot_temp - trop_temp) / lapse_rate
    return trop_alt + del_h


def compute_geop(q, t, sp, z, nlevs):
    from .consts import get_avec, get_bvec
    """Find geopotential for every gridpoint and model level.
    
    Adapted from: https://confluence.ecmwf.int/display/CKB/ERA5%3A+compute+geopotential+on+model+levels"""
    
    # Get the a and b coefficients for 137 level ECMWF
    avec = get_avec()
    bvec = get_bvec()
    
    rd = 287.06
    g_const = 9.81
    
    z_h = np.copy(q)
    z_h[:] = 0
    tmp_zh = z
    
    for curlev in list(range(nlevs - 1, -1, -1)):
        moist_t = t[curlev] * (1 + q[curlev])
        philev = avec[curlev - 1] + (bvec[curlev - 1] * sp)
        philevp1 = avec[curlev] + (bvec[curlev] * sp)
        if curlev == 1:
            dlogp = np.log(philevp1 / 0.1)
        else:
            dlogp = np.log(philevp1 / philev)
        moist_t = moist_t * rd
        tmp_zh = tmp_zh + (moist_t * dlogp)
        z_h[curlev] = tmp_zh / g_const

    return z_h


def set_rttov_dirs():
    """Set up the RTTOV main and data dirs."""
    
    # Note, you'll need to change this for your local PC
    rttov_top = '/network/aopp/apres/users/proud/Progs/rttov130/'
    ir_atlas_dir = rttov_top + '/emis_data/'
    brdf_atlas_dir = rttov_top + '/brdf_data/'
    top_rttov_coef = rttov_top + '/rtcoef_rttov13/'
    
    return rttov_top, ir_atlas_dir, brdf_atlas_dir, top_rttov_coef


def get_viirs_itemps():
    """Retrieve VIIRS I-band temps from file."""
    # These should not be relative paths, but this is generic
    i4 = np.load('../../data/Meteo/bt_section_i4.npy')
    i5 = np.load('../../data/Meteo/bt_section_i5.npy')

    return np.array(i4), np.array(i5)


def get_viirs_mtemps():
    """Retrieve VIIRS M-band temps from file."""
    # These should not be relative paths, but this is generic
    m14 = np.load('../../data/bt_section_m14.npy')
    m15 = np.load('../../data/bt_section_m15.npy')
    m16 = np.load('../../data/bt_section_m16.npy')

    return np.array(m14), np.array(m15), np.array(m16)


def init_atlas(ir_at_dir, brdf_at_dir, month):
    """Initiate the IR and BRDF atlases in pyRTTOV."""
    iratlas = pyrttov.Atlas()
    iratlas.AtlasPath = ir_at_dir
    iratlas.loadIrEmisAtlas(month, ang_corr=True)

    brdf_atlas = pyrttov.Atlas()
    brdf_atlas.AtlasPath = brdf_at_dir
    brdf_atlas.loadBrdfAtlas(month)
    brdf_atlas.IncSea = True

    return iratlas, brdf_atlas


def setup_surf_r_e(rttov_inst, n_profs, n_chans, i_rat, b_rat):
    """Set up the surface reflectance."""
    surf_er = np.zeros((4, n_profs, n_chans), dtype=np.float64)
    rttov_inst.SurfEmisRefl = surf_er
    surf_er[:, :, :] = -1.
    surf_er[0, :, :] = i_rat.getEmisBrdf(rttov_inst)
    surf_er[1, :, :] = b_rat.getEmisBrdf(rttov_inst)


def setup_rttov(chan_list, chanfile, cldfile, in_prof, n_threads):
    """Main RTTOV setup functions."""
    my_rttov = pyrttov.Rttov()
    my_rttov.FileCoef = chanfile
    my_rttov.FileSccld = cldfile
    my_rttov.Options.AddInterp = True
    my_rttov.Options.AddSolar = False
    my_rttov.Options.VerboseWrapper = True
    my_rttov.Options.AddClouds = True
    my_rttov.Options.CheckOpts = True
    my_rttov.Options.DoCheckinput = True
    my_rttov.Options.Nthreads = n_threads
    my_rttov.Options.NprofsPerCall = 5
    my_rttov.loadInst(chan_list)

    my_rttov.Profiles = in_prof
    in_prof.Icede[:, :] = 0

    return my_rttov, len(chan_list)


def make_profile(n_profiles, n_levels, ice_scheme, ice_paramer,
                 pres, temp, spechum, ozone, sat_angs, s2m,
                 skin, surf_type, surf_geom, datetimes):
    """Create a set of RTTOV profiles."""

    my_profiles = pyrttov.Profiles(n_profiles, n_levels)

    # We want clouds in g/m3
    my_profiles.MmrCldAer = 0

    # Set ice cloud types
    icecloud = np.array([ice_scheme, ice_paramer], dtype=np.int32).transpose()

    # We want kg/kg for O3 and Q
    my_profiles.GasUnits = 1

    # Set up the actual profiles, expanding for nprofs
    # This essentially copies the RTTOV example code
    my_profiles.P = expand2nprofiles(pres, n_profiles)
    my_profiles.T = expand2nprofiles(temp, n_profiles)
    my_profiles.Q = expand2nprofiles(spechum, n_profiles)
    my_profiles.O3 = expand2nprofiles(ozone, n_profiles)

    my_profiles.Angles = expand2nprofiles(sat_angs, n_profiles)
    my_profiles.S2m = expand2nprofiles(s2m, n_profiles)
    my_profiles.Skin = expand2nprofiles(skin, n_profiles)
    my_profiles.SurfType = expand2nprofiles(surf_type, n_profiles)
    my_profiles.SurfGeom = expand2nprofiles(surf_geom, n_profiles)

    my_profiles.DateTimes = expand2nprofiles(datetimes, n_profiles)

    my_profiles.IceCloud = expand2nprofiles(icecloud, n_profiles)
    my_profiles.Icede = expand2nprofiles(np.zeros(n_levels), n_profiles)
    my_profiles.Cirr = expand2nprofiles(np.zeros(n_levels), n_profiles)
    my_profiles.Cfrac = expand2nprofiles(np.zeros(n_levels), n_profiles)

    return my_profiles


def make_generic_profile(inf_m, inf_s, lat_storm, lon_storm, sat_angs,
                         icesch, iceparm, n_profiles, in_date):
    """Set up an RTTOV profile from ML and SFC files."""
    # Get surface and model level data
    la_pt, lo_pt, skint, sp, u10, v10 = get_sfc(inf_s, lat_storm, lon_storm)
    spechum, temp, geop_m, ozone, n_levels = get_ml(inf_m, la_pt, lo_pt)

    # Get pressure on model levels
    pres = get_pres(sp, n_levels)
    # Set up near surface data. Estimate some params
    s2m = np.array([sp / 100., temp[136], spechum[136], u10, v10, 100000.], dtype=np.float64).transpose()
    # Set up skin parameters
    skin = np.array([skint, 35., 0., 0., 3.0, 5.0, 15.0, 0.1, 0.3], dtype=np.float64).transpose()

    # Set up surface
    surftype = np.array([1, 1], dtype=np.int32).transpose()
    surfgeom = np.array([lat_storm, lon_storm, 0], dtype=np.float64).transpose()

    # Set up processing date
    date_in = [in_date.year, in_date.month, in_date.day,
               in_date.hour, in_date.minute, 0]
    datetimes = np.array(date_in, dtype=np.int32).transpose()

    # This makes the actual RTTOV profiles, which then need to be populated
    # with the cloud properties.
    my_profiles = make_profile(n_profiles, n_levels, icesch, iceparm,
                               pres, temp, spechum, ozone, sat_angs, s2m,
                               skin, surftype, surfgeom, datetimes)
    return my_profiles, n_levels


def get_data(inf_m, inf_s, lat_storm, lon_storm, sat_angs,
             icesch, iceparm, cirrange,
             cfrrange, upbin, in_date):
    """Retrieve data from file and use to construct profile."""

    # Number of profiles is defined by us, how many clouds to simulate
    n_profiles = cirrange[2] * cfrrange[2] * 6

    my_profiles, n_levels = make_generic_profile(inf_m, inf_s, lat_storm, lon_storm, sat_angs,
                                                 icesch, iceparm, n_profiles, in_date)

    # Counter for number of profiles
    npr = 0
    
    # This section of code sets the cloud properties that we use to
    # simulate the VIIRS scene in RTTOV.
    # upbin is the presumed location of the cloud. Because this is a manual
    # case study I just extrapolated a rough value from the temperature profile
    # coupled with the VIIRS BT and an estimate of OT lapse rate (-7.34K/K).
    
    # First, set everything below the cloud to be optically thick cloud
    my_profiles.Cirr[0:, upbin + 1:] = 3.
    my_profiles.Cfrac[0:, upbin + 1:] = 1.
    
    # Simulate some possible cloud tops. This is a poor way of doing it, but as we
    # are interested in an illustrative spectra + BTs rather than an actual cloud top 
    # properties retrieval it's acceptable here.
    
    for cirr in np.linspace(cirrange[0], cirrange[1], cirrange[2]):
        for cfrac in np.linspace(cfrrange[0], cfrrange[1], cfrrange[2]):
            # Set cloud top bin to specified cirrus IWC and fraction
            my_profiles.Cirr[npr, upbin] = cirr
            my_profiles.Cfrac[npr, upbin] = cfrac
            npr += 1
            
            # Now set some profiles as if the cloud is in two bins
            my_profiles.Cirr[npr, upbin] = cirr
            my_profiles.Cfrac[npr, upbin] = cfrac    
            my_profiles.Cirr[npr, upbin + 1] = cirr
            my_profiles.Cfrac[npr, upbin + 1] = cfrac  
            npr += 1 
            
            # Now a test for if there's a small above anvil cirrus plume
            my_profiles.Cirr[npr, upbin:upbin + 2] = 0
            my_profiles.Cfrac[npr, upbin:upbin + 2] = 0
            my_profiles.Cirr[npr, upbin + 2] = cirr
            my_profiles.Cfrac[npr, upbin + 2] = cfrac
            npr += 1
            
            # Now in case we got the cloud height bin incorrect
            my_profiles.Cirr[npr, upbin:upbin + 2] = 3
            my_profiles.Cfrac[npr, upbin:upbin + 2] = 1
            my_profiles.Cirr[npr, upbin + 3] = cirr
            my_profiles.Cfrac[npr, upbin + 3] = cfrac
            npr += 1

            # Cloud one layer lower
            my_profiles.Cirr[npr, upbin - 1] = cirr
            my_profiles.Cfrac[npr, upbin - 1] = cfrac
            npr += 1
            
            # Cloud takes up muliple layers
            my_profiles.Cirr[npr, upbin - 1:upbin + 2] = cirr
            my_profiles.Cfrac[npr, upbin - 1:upbin + 2] = cfrac
            npr += 1

    return my_profiles, n_levels, n_profiles
