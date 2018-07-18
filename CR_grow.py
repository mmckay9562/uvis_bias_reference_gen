import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import matplotlib.patches as mpatches
import os
from multiprocessing import Pool
from astropy.table import Table
import parser
import pandas as pd
import skimage.morphology as morph
from skimage.morphology import disk
import argparse


def CR_grow_main(fitsName1):
    """This function enlarges the cosmicray flag by a specified pixel radius.

    Author: Matthew Bourque & Myles McKay
    Date: April 30, 2017

    Function:
        CR_grow_main(fitsName1)

    Parameters:
        fitName1: 
            Calibrated .flt.fits files

    Output:
        Calibrated files with enlarged cosmic rays flags in the data quality arrays (ext 3 and 6) 
    Example:

        python CR_grow.py --path='/bias/flt/file/location/' 

    """
    #Open fits file
    hdu=fits.open(fitsName1, mode='update')

    dq3 = fits.getdata(fitsName1, ext=3)
    dq6 = fits.getdata(fitsName1, ext=6)

    dq3_orig = fits.getdata(fitsName1, ext=3)
    dq6_orig = fits.getdata(fitsName1, ext=6)

    #Sets all other data quality flags to zero
    dq3[np.where(dq3 != 8192)] = 0
    dq6[np.where(dq6 != 8192)] = 0

    #Enlarge cosmic rays
    dq3_grown = morph.dilation(dq3.byteswap().newbyteorder('='), disk(5))
    dq6_grown = morph.dilation(dq6.byteswap().newbyteorder('='), disk(5))

    #Populating enlarged cosmic ray extentions with other flagged pixels
    for i in np.arange(0,14):
        value=[2**i]
        dq3_grown[np.where(dq3_orig & value)] += value
        dq6_grown[np.where(dq6_orig & value)] += value

    #Writes a new file with enlarged cosmic ray flags
    hdu.close()
    hdulist = fits.open(fitsName1)
    hdulist[3].data = dq3_grown
    hdulist[6].data = dq6_grown
    hdulist[3].header['EXTNAME'] = 'DQ'
    hdulist[6].header['EXTNAME'] = 'DQ'
    hdulist[3].header['EXTVER'] = '1'
    hdulist[6].header['EXTVER'] = '2'
    hdulist.writeto('crrg_{}_flt.fits'.format(hdu[0].header['rootname'],overwrite=True))
    hdulist.close()


       

def parse_args():
    """Parses command line arguments.

    Parameters:
        nothing

    Returns:
        args : argparse.Namespace object
            An argparse object containing all of the added arguments.

    Outputs:
        nothing
    """

    #Create help string:
    path_help = 'Path to the folder with files to run tweakreg.'
    # Add arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', '-path', dest = 'path', action = 'store',
                        type = str, required = True, help = path_help)


    # Parse args:
    args = parser.parse_args()

    return args
# -------------------------------------------------------------------
if __name__ == '__main__':
    args = parse_args()
    os.chdir(args.path)
    list_of_files=glob.glob('*flt.fits')
    p=Pool(5)
    p.map(CR_grow_main,list_of_files)
    os.system('mkdir nongrown_flts')
    os.system('mv i*.fits nongrown_flts') 
   