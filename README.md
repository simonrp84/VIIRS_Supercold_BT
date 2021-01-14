# Supercold_BT_Paper

This repository provides the scripts and some of the data used in the 'Record-low cloud temperatures associated with a tropical deep convective event' paper submitted to AGU's GRL journal.

The `./Figures/` directory contains the output figures in vector and raster format, as submitted to GRL.
The `./data/` directory contains some input and processed output data:
 - `LOW_BT.csv` is a list of low BT values measured by VIIRS.
 - `Processed_BT.csv` is a similar list, but including VIIRS data from both NASA and NOAA. There are small differences in the BTs measured between data types.
 - `extre_frame.pkl` contains data on MODIS low BT events derived from data collected by the Aqua satellite.
 - `sonde.csv` contains processed radiosonde data from a location close to the overpass.
 - The `sat_csv` directory houses `csv` files detailing time series of BTs measured by various satellite sensors.
 - The `meteo` directory contains meteorological data derived from ERA5, as well as data used by RTTOV for one of the supplementary figures.
The `./Scripts/` directory contains all scripts used in the analysis presented in this manuscript. Please note that many of these scripts contain hardcoded directory names that will not be appropriate to your machine. You should change these to point to the correct locations. This is primarily an issue for the raw satellite data and, where possible, relative directory structures are used for data sources included in this repo.

# Usage
Most of these scripts are jupyter notebooks that can be used in any jupyter client. You will need a variety of libraries in order to run these scripts:
- Satpy: https://github.com/pytroll/satpy
- Pygac: https://github.com/pytroll/pygac
- Pandas: https://github.com/pandas-dev/pandas
- Matplotlib: https://github.com/matplotlib/matplotlib
- PyRTTOV: https://nwp-saf.eumetsat.int/site/software/rttov/
- Numpy: https://github.com/numpy/numpy
