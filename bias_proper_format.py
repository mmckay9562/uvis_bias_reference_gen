import os
from astropy.io import fits
import numpy as np
import argparse

def bias_format(path):

	os.chdir(path)
	os.system('cp /grp/hst/cdbs/iref/u6n1741hi_bia.fits ./')
	hdu_ref=fits.open('u6n1741hi_bia.fits',mode='update')
	hdu_new=fits.open('bias_crr_err_stacked_files.fits')
	
	hdu_ref.info()
	hdu_new.info()
	
	ref_chip1=hdu_ref[4].data
	err_chip1=hdu_ref[5].data
	chip1_dq=hdu_ref[6].data
	
	ref_chip2=hdu_ref[1].data
	err_chip2=hdu_ref[2].data
	chip2_dq=hdu_ref[3].data
	
	
	med_chip1=hdu_new[0].data
	new_err_chip1=hdu_new[1].data
	med_chip2=hdu_new[2].data
	new_err_chip2=hdu_new[3].data
	
	
	#--------------------------------------------------	
	#Replacing current data with new masterbias data
	#--------------------------------------------------	
	where_are_nans1 = np.isnan(med_chip1)
	where_are_nans2 = np.isnan(med_chip2)
	med_chip1[where_are_nans1] = 0
	med_chip2[where_are_nans2] = 0
	
	where_are_nans_err1 =np.isnan(new_err_chip1)
	where_are_nans_err2 =np.isnan(new_err_chip2)
	new_err_chip1[where_are_nans_err1] = 0
	new_err_chip2[where_are_nans_err2] = 0
	
	ref_chip1[19:2070,25:2073]=med_chip1[0:2051,0:2048]
	ref_chip1[19:2070,2133:4181]=med_chip1[0:2051,2048:]
	ref_chip2[0:2051,25:2073]=med_chip2[0:2051,0:2048]
	ref_chip2[0:2051,2133:4181]=med_chip2[0:2051,2048:]
	
	chip1_dq=np.zeros([2070,4206], dtype=np.int16)
	chip2_dq=np.zeros([2070,4206], dtype=np.int16)
	
	err_chip1[19:2070,25:2073]=new_err_chip1[0:2051,0:2048]
	err_chip1[19:2070,2133:4181]=new_err_chip1[0:2051,2048:]
	err_chip2[0:2051,25:2073]=new_err_chip2[0:2051,0:2048]
	err_chip2[0:2051,2133:4181]=new_err_chip2[0:2051,2048:]
	
	print(ref_chip1[0:2,24:26])	
	print(ref_chip1[0:2,2132:2134])	
	print(ref_chip2[18:20,24:26])
	print(ref_chip2[18:20,2132:2134])
	
	hdu_ref.close()
	hdu_new.close()
	
	hdu_ref=fits.open('u6n1741hi_bia.fits')
	hdu_ref.info()
	hdu_ref.close()
	
	os.system('mv u6n1741hi_bia.fits final_bias_ref.fits')


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
    bias_format(args.path)

