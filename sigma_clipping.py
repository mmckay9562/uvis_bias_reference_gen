from astropy.io import fits
import numpy as np
import glob
from astropy.table import Table
import argparse
import os
from wfc3tools import calwf3
import matplotlib.pyplot as plt
from multiprocessing import Pool
import shutil
from random import randint
import csv
import skimage.morphology as morph
from skimage.morphology import disk 
import pandas as pd
import matplotlib.dates as mdates

def sigma_clipping(path):
    """This function runs a 3 sigma clip of the average values of flt.fits files in a specified directory.

    Author: Myles McKay
    Date: April 30, 2017

    Function:
        sigma_clipping(path)

    Parameters:
        path: str
            The path to the directory containing raw.fits files.

    Output:
        New directory 1st_interation_sigma_clip, 2nd_interation_sigma_clip, 3rd_interation_sigma_clip
        containing clipped files and a comprihensive plots of the observation date vs. Average values of the bias files 

    Example:

        python sigma_clipping.py --path='/bias/flt/file/location/' 

    """

    #Sets the path to the specified directory
    os.chdir(path)
    #Makes directory for each iteration 
    os.system('mkdir 1st_interation_sigma_clip 2nd_interation_sigma_clip 3rd_interation_sigma_clip')

    #Reading in files and setting the local arrays
    list_of_files=sorted(glob.glob('*flt.fits'))
    file_number=[]
    file_date=[]
    file_mean1=[]
    file_mean2=[]
    #Masking the ext 1 & 4 with ext 3 & 6 then calcualting the staistics 
    for im in list_of_files:
        h=fits.open(im)
        name=h[0].header['rootname']
        date=h[0].header['date-obs']
        filter_name=h[0].header['filter']
        sci_chip1=h[4].data
        sci_chip2=h[1].data
        dq_chip1=h[6].data
        dq_chip2=h[3].data
        sci_chip1[dq_chip1 !=0]=np.nan
        sci_chip2[dq_chip2 !=0]=np.nan
        file_date=np.append(file_date,date)
        file_mean1=np.append(file_mean1,np.nanmean(sci_chip1))
        file_mean2=np.append(file_mean2,np.nanmean(sci_chip2))
    
    file_date = [pd.to_datetime(d,format='%Y-%m-%d') for d in file_date]
    Mean_c1=np.mean(file_mean1)
    STD_c1 =np.std(file_mean1)
    
    Mean_c2=np.mean(file_mean2)
    STD_c2 =np.std(file_mean2)
    
    upper_sigma_c1=Mean_c1 + 3.0*STD_c1
    lower_sigma_c1=Mean_c1 - 3.0*STD_c1
    
    upper_sigma_c2=Mean_c2 + 3.0*STD_c2
    lower_sigma_c2=Mean_c2 - 3.0*STD_c2
    
    #Plotting the Observations Date vs. Average value for UVIS chip1
    plt.scatter(file_date,file_mean1)
    plt.xlabel('Date-Obs')
    plt.ylabel('Mean Values')
    plt.title('Bias Chip1 Statistics\n Mean Value: {}  Std.Dev Values: {}'.format("%.3f" %Mean_c1, "%.3f" %STD_c1))
    plt.xticks(rotation=30)
    plt.axhline(y=upper_sigma_c1, xmin=-100,xmax=100,linewidth=2, color='red')
    plt.axhline(y=lower_sigma_c1, xmin=-100,xmax=100,linewidth=2, color='red')
    plt.axhline(y=Mean_c1, xmin=-100,xmax=100,linewidth=1, color='blue')
    plt.savefig('Statistics chip1 data plot.png')
#    plt.show()
    plt.clf()
    #Plotting the Observations Date vs. Average value for UVIS chip1
    plt.scatter(file_date,file_mean2)
    plt.xlabel('Date-Obs')
    plt.ylabel('Mean Values')
    plt.title('Bias Chip2 Statistics\n Mean Value: {}  Std.Dev Values: {}'.format("%.3f" %Mean_c2, "%.3f" %STD_c2))
    plt.xticks(rotation=30)
    plt.axhline(y=upper_sigma_c2, xmin=-100,xmax=100,linewidth=2, color='red')
    plt.axhline(y=lower_sigma_c2, xmin=-100,xmax=100,linewidth=2, color='red')
    plt.axhline(y=Mean_c2, xmin=-100,xmax=100,linewidth=1, color='blue')
    plt.savefig('Statistics chip2 data plot.png')
    plt.xlabel('Date-Obs')
    plt.ylabel('Mean Values')
 #   plt.show()
    plt.clf()   

    #Remove the sigma cliiped images from directory to new directory
    list_of_files= glob.glob('*flt.fits')
    for im in list_of_files:
        h = fits.open(im)
        sci_chip1=h[4].data
        sci_chip2=h[1].data
        dq_chip1=h[6].data
        dq_chip2=h[3].data
        sci_chip1[dq_chip1 !=0]=np.nan
        sci_chip2[dq_chip2 !=0]=np.nan
#        print(h[0].header['Rootname'],'   ','chip1','   ',np.nanmax(sci_chip1),'  ',np.nanmin(sci_chip1),'  ', np.nanmean(sci_chip1),'     ',np.nanstd(sci_chip1),'   ',np.nanmedian(sci_chip1))
        if np.nanmean(sci_chip1) >= upper_sigma_c1 or np.nanmean(sci_chip1) <= lower_sigma_c1:
            print(h[0].header['Rootname'],'   ','chip1','   ',np.nanmax(sci_chip1),'  ',np.nanmin(sci_chip1),'  ', np.nanmean(sci_chip1),'     ',np.nanstd(sci_chip1),'   ',np.nanmedian(sci_chip1))
            os.system('mv {}_flt.fits {}'.format(h[0].header['Rootname'],Folder))
        if np.nanmean(sci_chip2) >= upper_sigma_c2 or np.nanmean(sci_chip2) <= lower_sigma_c2: 
            print(h[0].header['Rootname'],'   ','chip2','   ',np.nanmax(sci_chip2),'  ',np.nanmin(sci_chip2),'  ', np.nanmean(sci_chip2),'     ',np.nanstd(sci_chip2),'   ',np.nanmedian(sci_chip2))
            os.system('mv {}_flt.fits {}'.format(h[0].header['Rootname'],Folder))
        h.close()
    








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

    Folder = './1st_interation_sigma_clip'
    sigma_clipping(args.path)
    os.system('mv *.png ./1st_interation_sigma_clip')

    Folder = './2nd_interation_sigma_clip'
    sigma_clipping(args.path)
    os.system('mv ./*.png ./2nd_interation_sigma_clip')

    Folder = './3rd_interation_sigma_clip'
    sigma_clipping(args.path)
    os.system('mv ./*.png ./3rd_interation_sigma_clip')










