import numpy as np

def do_writing(outfid, data, dtstr):
    """Save some satellite data stats to a file."""
    ndata = np.copy(data)
    pts = (ndata > 273.15-65.15).nonzero()
    ndata[pts] = np.nan
    outstr = dtstr + ','
    outstr = outstr + str(np.nanmin(ndata)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 1)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 2)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 5)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 10)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 25)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 50)) + ','
    outstr = outstr + str(np.nanpercentile(ndata, 90)) + ','
    outstr = outstr + str(np.nanmax(ndata)) + ','
    outstr = outstr + str(np.nanmean(ndata)) + ','
    outstr = outstr + str(np.nanmedian(ndata)) + ','
    outstr = outstr + str(np.nanstd(ndata)) + '\n'
    outfid.write(outstr)
