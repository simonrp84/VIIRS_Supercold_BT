"""
Create a plot showing atmospheric temperature alongside different satellite bands.

This script runs RTTOV for an IASI simulation of the supercold BT. The measured BTs
from a variety of other satellites are then over-plotted on the IASI spectra for comparison.

"""
from pyspectral.rsr_reader import RelativeSpectralResponse
from pyspectral.utils import get_wave_range
import Meteo_Scripts.proc_prof as pp
import Meteo_Scripts.load_sat as ls
import matplotlib.pyplot as plt
from datetime import datetime
import xarray as xr
import numpy as np
import warnings
import pytz

warnings.filterwarnings('ignore')

np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

# Set up some RTTOV locations, these are hardcoded and should be changed
ir_atlas_dir = '/network/aopp/apres/users/proud/Progs/rttov130/emis_data/'
brdf_atlas_dir = '/network/aopp/apres/users/proud/Progs/rttov130/brdf_data/'
top_rttov_coef = '/network/aopp/apres/users/proud/Progs/rttov130/rtcoef_rttov13/'

iasi_rtc = top_rttov_coef + 'rttov9pred101L/rtcoef_metop_2_iasi_allgas.H5'
iasi_cld = top_rttov_coef + 'cldaer_visir/sccldcoef_metop_2_iasi.H5'

# Storm position and solar angles, which aren't really necessary (this is IR only) but
# best to include the real angles in case VIS simulations would be useful at some point.
sza = -62.64299
saa = 165.6918

storm_lat = -3.2609048
storm_lon = 163.2608


storm_time = datetime(2018, 12, 29, 13, 30, tzinfo=pytz.utc)

# Below are pre-defined BT values for the storm coordinates and time given above.
# Also included are the wavelength ranges required for the plotting later on.

# This threshold is used to define min/max wavelength
rsr_thresh = 0.15

# NOAA-20 / VIIRS first
rsr = RelativeSpectralResponse('NOAA-20', 'viirs')
viirsi_data = {'I04': [208.000] + get_wave_range(rsr.rsr['I4']['det-1'], rsr_thresh),
               'I05': [161.955] + get_wave_range(rsr.rsr['I5']['det-1'], rsr_thresh)}

viirsm_data = {'M12': [203.000] + get_wave_range(rsr.rsr['M12']['det-1'], rsr_thresh),
               'M13': [191.920] + get_wave_range(rsr.rsr['M13']['det-1'], rsr_thresh),
               'M14': [165.503] + get_wave_range(rsr.rsr['M14']['det-1'], rsr_thresh),
               'M15': [163.733] + get_wave_range(rsr.rsr['M15']['det-1'], rsr_thresh),
               'M16': [165.274] + get_wave_range(rsr.rsr['M16']['det-1'], rsr_thresh)}

# Now GOES-17 / ABI
rsr = RelativeSpectralResponse('GOES-17', 'abi')
abi_data = {'C07': [197.534] + get_wave_range(rsr.rsr['ch7']['det-1'], rsr_thresh),
            'C08': [182.390] + get_wave_range(rsr.rsr['ch8']['det-1'], rsr_thresh),
            'C09': [178.648] + get_wave_range(rsr.rsr['ch9']['det-1'], rsr_thresh),
            'C10': [178.982] + get_wave_range(rsr.rsr['ch10']['det-1'], rsr_thresh),
            'C11': [175.723] + get_wave_range(rsr.rsr['ch11']['det-1'], rsr_thresh),
            'C12': [209.068] + get_wave_range(rsr.rsr['ch12']['det-1'], rsr_thresh),
            'C13': [174.917] + get_wave_range(rsr.rsr['ch13']['det-1'], rsr_thresh),
            'C14': [173.923] + get_wave_range(rsr.rsr['ch14']['det-1'], rsr_thresh),
            'C15': [173.160] + get_wave_range(rsr.rsr['ch15']['det-1'], rsr_thresh),
            'C16': [175.779] + get_wave_range(rsr.rsr['ch16']['det-1'], rsr_thresh)}

# FY-4A / AGRI
rsr = RelativeSpectralResponse('FY-4A', 'agri')
agri_data = {'C07': [200.0] + get_wave_range(rsr.rsr['ch7']['det-1'], rsr_thresh*100.),
             'C08': [193.2] + get_wave_range(rsr.rsr['ch8']['det-1'], rsr_thresh*100.),
             'C09': [181.1] + get_wave_range(rsr.rsr['ch9']['det-1'], rsr_thresh*100.),
             'C10': [182.0] + get_wave_range(rsr.rsr['ch10']['det-1'], rsr_thresh*100.),
             'C11': [183.1] + get_wave_range(rsr.rsr['ch11']['det-1'], rsr_thresh*100.),
             'C12': [178.0] + get_wave_range(rsr.rsr['ch12']['det-1'], rsr_thresh*100.),
             'C13': [178.4] + get_wave_range(rsr.rsr['ch13']['det-1'], rsr_thresh*100.)}

# Himawari-8 / AHI
rsr = RelativeSpectralResponse('Himawari-8', 'ahi')
ahi_data = {'B07': [195.136] + get_wave_range(rsr.rsr['ch7']['det-1'], rsr_thresh),
            'B08': [176.996] + get_wave_range(rsr.rsr['ch8']['det-1'], rsr_thresh),
            'B09': [174.775] + get_wave_range(rsr.rsr['ch9']['det-1'], rsr_thresh),
            'B10': [174.951] + get_wave_range(rsr.rsr['ch10']['det-1'], rsr_thresh),
            'B11': [172.630] + get_wave_range(rsr.rsr['ch11']['det-1'], rsr_thresh),
            'B12': [201.853] + get_wave_range(rsr.rsr['ch12']['det-1'], rsr_thresh),
            'B13': [172.503] + get_wave_range(rsr.rsr['ch13']['det-1'], rsr_thresh),
            'B14': [171.572] + get_wave_range(rsr.rsr['ch14']['det-1'], rsr_thresh),
            'B15': [171.881] + get_wave_range(rsr.rsr['ch15']['det-1'], rsr_thresh),
            'B16': [175.106] + get_wave_range(rsr.rsr['ch16']['det-1'], rsr_thresh)}

# Otherwise, load data from files:
# This is slow, so commented out unless needed (different lat/lon, f.ex)
# topdir = '/gf2/eodg/SRP001_PROUD_TURBREP/Supercold/'
# sat_data_dict = ls.retrieve_satellite_chans(topdir, storm_lat, storm_lon)

# VIIRS viewing geometry.
virs_vza = 68.055565
virs_vaa = 99.93242

ice_sch = 1  # Baum
ice_parm = 1

# tropopause bin in the NWP profile, manually determined
trop_bin = 55

# Number of threads for RTTOV, set to 7 for local desktop (8 thread machine)
n_threads = 7

satangs_virs = np.array([virs_vza, virs_vaa, sza, saa], dtype=np.float64).transpose()

# Pre-made meteorological profiles from ERA5
infile_ml = '../data/Meteo/ML.nc'
infile_sfc = '../data/Meteo/SFC.nc'

# And the pre-determined optimal cloud properties
prof_arr = xr.open_dataset('../data/Meteo/good_profs_iband.nc')
cirr_arr = prof_arr['CIRR'].values
cfra_arr = prof_arr['CFRAC'].values

# How many profiles do we need? 81 for a 9x9 array
nprofs = len(cfra_arr)
print(nprofs)

# This is a generic IASI-like sensor at the VIIRS geometry
in_prof_iasi, nlevs = pp.make_generic_profile(infile_ml,
                                              infile_sfc,
                                              storm_lat,
                                              storm_lon,
                                              satangs_virs,
                                              ice_sch,
                                              ice_parm,
                                              1,
                                              storm_time)

# Modify the profile based on a typical atmospheric lapse rate. Only points
# above the tropopause are modified. This accounts for overshooting tops.
for i in range(trop_bin, trop_bin - 3, -1):
    in_prof_iasi.T[:, i] = in_prof_iasi.T[:, i + 1] - 7.34

# Now scale everything
sizer = int(np.sqrt(nprofs))
cpix = int((sizer - 1) / 2.)
scir = np.reshape(cirr_arr, (sizer, sizer, 137))
scfr = np.reshape(cfra_arr, (sizer, sizer, 137))
in_prof_iasi.Cirr[0] = scir[cpix, cpix, :]
in_prof_iasi.Cfrac[0] = scfr[cpix, cpix, :]

# Channel lists for various satellites
chan_list_iasi = np.arange(1, 8461)

# Create RTTOV objects
iasi_i_rttov, nchan_iasi = pp.setup_rttov(chan_list_iasi, iasi_rtc, iasi_cld, in_prof_iasi, n_threads)
irAtlas, brdfAtlas = pp.init_atlas(ir_atlas_dir, brdf_atlas_dir, storm_time.month)

# Initialise surface conditions
# RTTOV needs this, but as we're looking at a thick cloud it's essentially unused.
pp.setup_surf_r_e(iasi_i_rttov, 1, nchan_iasi, irAtlas, brdfAtlas)

# Run RTTOV itself.
iasi_i_rttov.runDirect()
print("Done IASI")

# Compute wavelengths from the per-sensor wavenumbers
iasi_wvl = 10000 / iasi_i_rttov.WaveNumbers
# Remove singular length dimensions and convert to Kelvin
iasi_i_bt = np.squeeze(iasi_i_rttov.BtRefl)
iasi_out_c = iasi_i_bt - 273.15

# Now, finally, we create the output plot

# Set up some plot variables
# Axis limits, x is micron, y is temperature in C
x_min = 7.8
x_max = 13.8
y_min = -115
y_max = -62

# DPI of output
dpier = 100

# Colours to plot various instruments
agri_col = '#ff7f0e'
abi_col = '#2ca02c'
ahi_col = '#d62728'
viri_col = '#9467bd'
virm_col = '#8c564b'

# Plotting parameters
capsize = 10
marker_sim = 'x'
marker_act = 'o'
msize_sim = 60
msize_act = 80
act_dict = {'zorder': 4, 'marker': marker_act, 's': msize_sim}
sim_dict = {'zorder': 2, 'marker': marker_sim, 's': msize_sim}
err_dict = {'zorder': 2, 'linestyle': "None", 'capsize': capsize, 'linewidth': 0.8}

# Set up the figure
fig, ax = plt.subplots(figsize=(16, 9))

# Plot simulated data
plt.plot(iasi_wvl, iasi_out_c, zorder=1, label='IASI', c='grey', linewidth=1)

# Plot actual data with errorbars
agri_data = np.array(list(agri_data.values()))
agri_data[:, 1] = np.abs(agri_data[:, 1] - agri_data[:, 2])
agri_data[:, 3] = np.abs(agri_data[:, 3] - agri_data[:, 2])
plt.errorbar(agri_data[:, 2],
             agri_data[:, 0] - 273.15,
             xerr=np.vstack((agri_data[:, 1], agri_data[:, 3])),
             c=agri_col,
             **err_dict)

abi_data = np.array(list(abi_data.values()))
abi_data[:, 1] = np.abs(abi_data[:, 1] - abi_data[:, 2])
abi_data[:, 3] = np.abs(abi_data[:, 3] - abi_data[:, 2])
plt.errorbar(abi_data[:, 2],
             abi_data[:, 0] - 273.15,
             xerr=np.vstack((abi_data[:, 1], abi_data[:, 3])),
             c=abi_col,
             **err_dict)

ahi_data = np.array(list(ahi_data.values()))
ahi_data[:, 1] = np.abs(ahi_data[:, 1] - ahi_data[:, 2])
ahi_data[:, 3] = np.abs(ahi_data[:, 3] - ahi_data[:, 2])
plt.errorbar(ahi_data[:, 2],
             ahi_data[:, 0] - 273.15,
             xerr=np.vstack((ahi_data[:, 1], ahi_data[:, 3])),
             c=ahi_col,
             **err_dict)

viirsi_data = np.array(list(viirsi_data.values()))
viirsi_data[:, 1] = np.abs(viirsi_data[:, 1] - viirsi_data[:, 2])
viirsi_data[:, 3] = np.abs(viirsi_data[:, 3] - viirsi_data[:, 2])
plt.errorbar(viirsi_data[:, 2],
             viirsi_data[:, 0] - 273.15,
             xerr=np.vstack((viirsi_data[:, 1], viirsi_data[:, 3])),
             c=viri_col,
             **err_dict)

viirsm_data = np.array(list(viirsm_data.values()))
viirsm_data[:, 1] = np.abs(viirsm_data[:, 1] - viirsm_data[:, 2])
viirsm_data[:, 3] = np.abs(viirsm_data[:, 3] - viirsm_data[:, 2])
plt.errorbar(viirsm_data[:, 2],
             viirsm_data[:, 0] - 273.15,
             xerr=np.vstack((viirsm_data[:, 1], viirsm_data[:, 3])),
             c=virm_col,
             **err_dict)

# Plot pixel values as points with labels.
plt.scatter(agri_data[:, 2], agri_data[:, 0]-273.15, c=agri_col, **act_dict, label='AGRI')
plt.scatter(abi_data[:, 2], abi_data[:, 0]-273.15, c=abi_col, **act_dict, label='ABI')
plt.scatter(ahi_data[:, 2], ahi_data[:, 0]-273.15, c=ahi_col, **act_dict, label='AHI')
plt.scatter(viirsi_data[:, 2], viirsi_data[:, 0]-273.15, c=viri_col, **act_dict, label='VIIRS_I')
plt.scatter(viirsm_data[:, 2], viirsm_data[:, 0]-273.15, c=virm_col, **act_dict, label='VIIRS_M')

# Finish setting up the plot
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel('Wavelength ($\mu$m)', fontsize=16)
ax.set_ylabel('Temperature ($\degree$C)', fontsize=16)
ax.set_facecolor('white')

plt.grid(b=True, which='major', color='black', linewidth=0.5, linestyle='-')

plt.legend(loc=1, framealpha=1., fontsize=14)
plt.tight_layout()

plt.savefig('../Figures/Figure_S1_Spectral.png', dpi=dpier, facecolor='white')
plt.savefig('../Figures/Figure_S1_Spectral.svg', dpi=dpier, facecolor='white')
