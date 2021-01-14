"""Load satellite data for the various instruments to plot on Fig_S1."""

from satpy import Scene
from glob import glob
import numpy as np


def retrieve_satellite_chans(top_dir, storm_lat, storm_lon):
    """Load the satellite data for various instruments.

    For each sensor a dict is created with the band name as a key.
    Each entry in the dict is a list of:
     - Minimum BT for that band
     - Minimum wavelength for the band
     - Central wavelength for the band
     - Maximum wavelength for the band

     The wavelength values are retrieved from Satpy.
     """

    # Glob patterns for loading sat data
    abi_globber = top_dir + 'ABI/*s20183631330343*'
    agri_globber = top_dir + 'AGRI/*20181229130000*'
    ahi_globber = top_dir + 'AHI_F/*20181229_1330*'
    viirsi_globber = top_dir + 'VIIRSI/*d20181229_t1336*'
    viirsm_globber = top_dir + 'VIIRSM2/*d20181229_t1336*'

    # Channels to load for each instrument, in the 3.5 -> 14 micron range (IR).
    abi_chans = ['C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'C13', 'C14', 'C15', 'C16']
    agri_chans = ['C07', 'C08', 'C09', 'C10', 'C11', 'C12', 'C13']
    ahi_chans = ['B07', 'B08', 'B09', 'B10', 'B11', 'B12', 'B13', 'B14', 'B15', 'B16']
    viirsi_chans = ['I04', 'I05']
    viirsm_chans = ['M12', 'M13', 'M14', 'M15', 'M16']

    viirs_i_data = load_sat_polar('viirs_sdr', viirsi_globber, viirsi_chans, storm_lat, storm_lon, 1)
    viirs_m_data = load_sat_polar('viirs_sdr', viirsm_globber, viirsm_chans, storm_lat, storm_lon, 1)
    abi_data = load_sat_geo('abi_l1b', abi_globber, abi_chans, storm_lat, storm_lon, 1)
    agri_data = load_sat_geo('agri_l1', agri_globber, agri_chans, storm_lat, storm_lon, 1)
    ahi_data = load_sat_geo('ahi_hsd', ahi_globber, ahi_chans, storm_lat, storm_lon, 1)

    data_dict = {'VIIRSI': viirs_i_data,
                 'VIIRSM': viirs_m_data,
                 'ABI': abi_data,
                 'AGRI': agri_data,
                 'AHI': ahi_data}
    return data_dict


def load_sat_polar(in_reader, globber, chans, latpos, lonpos, deller):
    """Load data for a given sensor, find BT at particular lat/lon."""
    storm_bbox = [lonpos-deller, latpos-deller,
                  lonpos+deller, latpos+deller]

    files = glob(globber)
    if len(files) < 1:
        return None
    scn = Scene(files, reader=in_reader)
    scn.load(chans)

    lons = None
    lats = None
    output_dict = {}
    for chan in chans:
        try:
            if lons is None or lats is None:
                lons, lats = scn[chan].attrs['area'].get_lonlats()
                try:
                    lats = lats.compute()
                    lons = lons.compute()
                except:
                    pass
            btdata = np.array(scn[chan])
            pts = (lons < storm_bbox[0]).nonzero()
            btdata[pts] = np.nan
            pts = (lats < storm_bbox[1]).nonzero()
            btdata[pts] = np.nan
            pts = (lons > storm_bbox[2]).nonzero()
            btdata[pts] = np.nan
            pts = (lats > storm_bbox[3]).nonzero()
            btdata[pts] = np.nan
            output_dict[chan] = [np.nanmin(btdata),
                                 scn[chan].wavelength.min,
                                 scn[chan].wavelength.central,
                                 scn[chan].wavelength.max]
        except KeyError:
            output_dict[chan] = [None, None, None, None]
    return output_dict


def load_sat_geo(in_reader, globber, chans, latpos, lonpos, deller):
    """Load data for a given sensor, find BT at particular lat/lon."""
    storm_bbox = [lonpos-deller, latpos-deller,
                  lonpos+deller, latpos+deller]

    files = glob(globber)
    if len(files) < 1:
        return None
    scn = Scene(files, reader=in_reader)
    scn.load(chans)

    scn2 = scn.crop(ll_bbox=storm_bbox)
    output_dict = {}
    for chan in chans:
        try:
            btdata = np.array(scn2[chan])
            output_dict[chan] = [np.nanmin(btdata),
                                 scn[chan].wavelength.min,
                                 scn[chan].wavelength.central,
                                 scn[chan].wavelength.max]
        except KeyError:
            output_dict[chan] = [None, None, None, None]
    return output_dict
