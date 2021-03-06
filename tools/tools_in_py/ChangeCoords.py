import numpy as np
from astropy.io import fits
from astropy import wcs
import sys

DATA_PATH = "/home/queenmab/DATA/M31/DM/jfactor/"

filename = DATA_PATH+"annihil_rs28_gamma12D_FOVdiameter20.0deg_alphaint0.36deg_nside256-JFACTOR-Jsmooth-image.fits"

hdul = fits.open(filename)

#Create a new WCS object, number of axes set from start

w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [105.0, 105.0]
w.wcs.cdelt = np.array([0.02863961813842482, 0.02863961813842482])
w.wcs.crval = [10.76, 41.19]
w.wcs.ctype = ["RA---CAR", "DEC--CAR"]

w.wcs.set_pv([(2, 1, 45.0)])

# Now, write out the WCS object as a FITS header
header = w.to_header()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = fits.PrimaryHDU(header=header,data=hdul[0].data)
# Save to FITS file
hdu.writeto(filename,overwrite=True)
