

from __future__ import print_function, division, absolute_import

import os

import numpy as np
from astropy.io import fits

from RMS.Formats.FFStruct import FFStruct



def read(directory, filename, array=False):
    """ Read a FF structure from a FITS file. 
    
    Arguments:
        directory: [str] Path to directory containing file
        filename: [str] Name of FF*.fits file (either with FF and extension or without)

    Keyword arguments:
        _array: [ndarray] True in order to populate structure's array element (default is False)
    
    Return:
        [ff structure]

    """

    # Make sure the file starts with "FF_"
    if filename[:3] == "FF_":
        fid = open(os.path.join(directory, filename), "rb")
    else:
        fid = open(os.path.join(directory, "FF_" + filename + ".fits"), "rb")

    # Init an empty FF structure
    ff = FFStruct()

    # Read in the FITS
    hdulist = fits.open(fid)

    # Read the header
    head = hdulist[0].header

    # Read in the data from the header
    ff.nrows = head['NAXIS2']
    ff.ncols = head['NAXIS1']
    ff.nbits = head['NBITS']
    ff.nframes = head['NFRAMES']
    ff.first = head['FIRST']
    ff.camno = head['CAMNO']
    ff.fps = head['FPS']

    # Read in the image data
    ff.maxpixel = np.flipud(hdulist[0].data[2].astype('uint8'))
    ff.maxframe = np.flipud(hdulist[0].data[3].astype('uint8'))
    ff.avepixel = np.flipud(hdulist[0].data[0].astype('uint8'))
    ff.stdpixel = np.flipud(hdulist[0].data[1].astype('uint8'))

    if array:
        ff.array = np.dstack([ff.maxpixel, ff.maxframe, ff.avepixel, ff.stdpixel])

    # CLose the FITS file
    hdulist.close()

    return ff



def write(ff, directory, filename):
    """ Write a FF structure to a FITS file in specified directory.
    
    Arguments:
        ff: [ff bin struct] FF bin file loaded in the FF structure
        directory: [str] path to the directory where the file will be written
        filename: [str] name of the file which will be written
    
    Return:
        None

    """

    # Make sure the file starts with "FF"
    if filename[:3] == "FF_":
        file_path = os.path.join(directory, filename)

    else:
        file_path = os.path.join(directory, "FF_" + filename + ".fits")

    # Create a new FITS file
    
    # Create the header
    head = fits.Header()
    head['NBITS'] = ff.nbits
    head['DATE-OBS'] = ff.nfd
    head['MJD-OBS'] = ff.mjd
    head['EXPTIME'] = ff.dt[-1]-ff.dt[0]
    head['NFRAMES'] = ff.nframes
    head['FIRST'] = ff.first
    head['CAMNO'] = ff.camno
    head['FPS'] = ff.fps
    head['CRPIX1'] = float(ff.ncols)/2.0
    head['CRPIX2'] = float(ff.nrows)/2.0
    head['CRVAL1'] = 0.0
    head['CRVAL2'] = 0.0
    head['CD1_1'] = 1.0
    head['CD1_2'] = 0.0
    head['CD2_1'] = 0.0
    head['CD2_2'] = 1.0
    head['CTYPE1'] = "RA---TAN"
    head['CTYPE2'] = "DEC--TAN"
    head['CUNIT1'] = "deg"
    head['CUNIT2'] = "deg"
    head['CRRES1'] = 0.0
    head['CRRES2'] = 0.0
    head['EQUINOX'] = 2000.0
    head['RADECSYS'] = "ICRS"
    head['COSPAR'] = 4171
    head['OBSERVER'] = "Cees Bassa"
    for i in xrange(len(ff.dt)):
        head['DT%04d'%i] = ff.dt[i]
    for i in xrange(10):
        head['DUMY%03d'%i] = 0.0

    # Reshape, reformat and flip
    ff.maxpixel, ff.maxframe, ff.avepixel, ff.stdpixel = np.split(ff.array, 4, axis=0)
    zmax=np.flipud(ff.maxpixel[0].astype('float32'))
    znum=np.flipud(ff.maxframe[0].astype('float32'))
    zavg=np.flipud(ff.avepixel[0].astype('float32'))
    zstd=np.flipud(ff.stdpixel[0].astype('float32'))
    z=np.array([zavg,zstd,zmax,znum])
    
    # Create the primary part
    hdu = fits.PrimaryHDU(header=head,data=z)
    
    # Save the FITS
    hdu.writeto(file_path,overwrite=True)

if __name__ == "__main__":

    dir_path = '.'
    file_name = 'FF_test.fits'

    wid = 720
    ht = 576


    ff = FFStruct()

    ff.ncols = wid
    ff.nrows = ht

    # ff.maxpixel = np.zeros((ht, wid), dtype=np.uint8)
    # ff.avepixel = np.zeros((ht, wid), dtype=np.uint8) + 10
    # ff.stdpixel = np.zeros((ht, wid), dtype=np.uint8) + 20
    # ff.maxframe = np.zeros((ht, wid), dtype=np.uint8) + 30

    maxpixel = np.zeros((ht, wid), dtype=np.uint8)
    avepixel = np.zeros((ht, wid), dtype=np.uint8) + 10
    stdpixel = np.zeros((ht, wid), dtype=np.uint8) + 20
    maxframe = np.zeros((ht, wid), dtype=np.uint8) + 30

    ff.array = np.stack([maxpixel, maxframe, avepixel, stdpixel], axis=0)

    # Write the FF to FITS
    write(ff, dir_path, file_name)

    # Read the FITS
    ff = read(dir_path, file_name)

    print(ff)
    print(ff.maxpixel)
    print(ff.maxframe)
    print(ff.avepixel)
    print(ff.stdpixel)
