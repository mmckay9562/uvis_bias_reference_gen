import numpy as np
import matplotlib
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import argparse
import shutil
import os
from random import randint
import csv
from astropy.table import Table
from wfc3tools import calwf3


def create_asn_file(path):
    """This function create an association table of all raw.fits files in the specified directory.

    Author: Heather Kurtz & Myles McKay
    Date: April 30, 2017

    Parameters:
        path: str
            The path to the directory containing raw.fits files.

    Output:
        calibrated flt.fits files
        New directory containing .tra, raw.fits, .csv  and asn.fits files 

    Example:

        python create_asn_file.py --path='/raw/file/location/' 

    """

#Renames all files with in the given directory to associate
    os.chdir(path) 
    base_path = path
    name_matches=[]
    check=[]
    number=[]
    
    setOfNumbers=set()
    while len(setOfNumbers) < 600:
        setOfNumbers.add(str(randint(111,999)))
    
    my_list=list(setOfNumbers)
    
    list_of_files=glob.glob('*raw.fits')
    
    for i in range(len(list_of_files)):
        files=list_of_files[i]
        curentName= os.path.join(base_path, files)
        name=files[9:]
        num=my_list[i]
        new_name='imam11'+num+name
        check.append(new_name)
        name_matches.append((files,new_name))
        fileName= os.path.join(base_path, new_name)
        os.rename(curentName, fileName)
    
       
    with open('name_matches.csv','w') as out:
        csv_out=csv.writer(out)
        csv_out.writerow(['origanal','new'])
        for row in name_matches:
            csv_out.writerow(row)
#-------------------------------------------------------------------------------------------------------------------------------------------------

#Making the association
    files_wo_asn=[]
    files_con_asn=[]
    asns=[]
    obs=[]
    
    
    os.chdir(path)
    #Uses the vists to set new associations
    visits=[]
    list_of_files=glob.glob('*raw.fits')
    for files in list_of_files:
        name=files
        visit=name[:6]
        visits.append(visit)
        hdu=fits.open(files, mode='update')
        hdu[0].header['PCTECORR'] = 'OMIT'
        hdu.close()

    unique_visits= set(visits)
    obs_date=[]
    print (unique_visits)
    for items in unique_visits:
        asn_files=[]
        for j in range(len(list_of_files)):
            files=list_of_files[j]
            comp=files[:6]
            write_file=files[:9]
            if items==comp:
                asn_files.append(write_file)
                output=files[:6]+'010_asn.fits'
                files=files
            else:
                continue
                
        if len(asn_files)>1:
            row = [i+1 for i in range(len(asn_files))]
            memtype = ['EXP-CRJ' for f in asn_files]
            memprsnt = [True for f in asn_files]
            xoffset = [0. for f in asn_files]
            yoffset = [0. for f in asn_files]
            xdelta = [0. for f in asn_files]
            ydelta = [0. for f in asn_files]
            rotation = [0. for f in asn_files]
            scale = [1. for f in asn_files]
            # Change the last row to have a memtype of PROD-CRJ instead of EXP-CRJ
            #memtype[-1] = 'PROD-CRJ'
        
            # Build the table
            columns = [row, asn_files, memtype, memprsnt, xoffset, yoffset, xdelta, ydelta, rotation, scale]
            names = ['row', 'MEMNAME', 'MEMTYPE', 'MEMPRSNT', 'XOFFSET', 'YOFFSET', 'XDELTA', 'YDELTA', 'ROTATION', 'SCALE']
            asn_table = Table(columns, names=names)
            
            # Save the table to a FITS file
            asn_table.write(output)
            
            # Update header to have required keywords
            hdulist = fits.open(output, mode='update')
            hdulist[0].header['INSTRUME'] = 'WFC3'
            hdulist[0].header['DETECTOR'] = 'UVIS'
            hdulist[0].header['PCTECORR'] = 'OMIT'
            hdulist.close()
            
        else:
            continue
    
    calwf3(output)


#----------------------------------------------------------------------------------------------------------------------------------------------------

#Changing the filename back to original filename    
    os.chdir(path)
    
    name_new=[]
    old_name=[]
    name_matches=[]
    
    fileSB=open('name_matches.csv', 'r')
        
    for line in fileSB:
        q= [d for d in line.split(',')]
        old_name.append(q[0])
        name_new.append(q[1]) #reads second columns
    
    
    list_of_raw=glob.glob('*raw.fits')
    
    for i in range(len(list_of_raw)):
        #for num in setOfNumbers:
        File=list_of_raw[i]
        files=File[:9]
        for j in range(len(name_new)):
            name_comp=name_new[j]
            comp_name=name_comp[:9]
            if files==comp_name:
                curentName= os.path.join(base_path, File)
                name=File[9:]
                old=old_name[j]
                OLD=old[:9]
                file_name=OLD+name
                #check.append(file_name)
                name_matches.append((File,file_name))
                fileName= os.path.join(base_path, file_name)
                os.rename(curentName, fileName)
            else:
                continue
    
    list_of_flt=glob.glob('*flt.fits')
    
    for i in range(len(list_of_flt)):
        #for num in setOfNumbers:
        File=list_of_flt[i]
        files=File[:9]
        for j in range(len(name_new)):
            name_comp=name_new[j]
            comp_name=name_comp[:9]
            if files==comp_name:
                curentName= os.path.join(base_path, File)
                name=File[9:]
                old=old_name[j]
                OLD=old[:9]
                file_name=OLD+name
                #check.append(file_name)
                name_matches.append((File,file_name))
                fileName= os.path.join(base_path, file_name)
                os.rename(curentName, fileName)
            else:
                continue

    list_of_flc=glob.glob('*flc.fits')
    
    for i in range(len(list_of_flc)):
        #for num in setOfNumbers:
        File=list_of_flc[i]
        files=File[:9]
        for j in range(len(name_new)):
            name_comp=name_new[j]
            comp_name=name_comp[:9]
            if files==comp_name:
                curentName= os.path.join(base_path, File)
                name=File[9:]
                old=old_name[j]
                OLD=old[:9]
                file_name=OLD+name
                #check.append(file_name)
                name_matches.append((File,file_name))
                fileName= os.path.join(base_path, file_name)
                os.rename(curentName, fileName)
            else:
                continue

    list_of_tra=glob.glob('*.tra')
    
    for i in range(len(list_of_tra)):
        #for num in setOfNumbers:
        File=list_of_tra[i]
        files=File[:9]
        for j in range(len(name_new)):
            name_comp=name_new[j]
            comp_name=name_comp[:9]
            if files==comp_name:
                curentName= os.path.join(base_path, File)
                name=File[9:]
                old=old_name[j]
                OLD=old[:9]
                file_name=OLD+name
                #check.append(file_name)
                name_matches.append((File,file_name))
                fileName= os.path.join(base_path, file_name)
                os.rename(curentName, fileName)
            else:
                continue
    fileSB.close()
        
    with open('matches.csv','w') as out:
        csv_out=csv.writer(out)
        csv_out.writerow(['origanal','new'])
        for row in name_matches:
            csv_out.writerow(row)
    
    os.system('mkdir raw_crj_tra_csv_files')
    os.system('mv *raw.fits raw_crj_tra_csv_files')
    os.system('mv *.tra raw_crj_tra_csv_files')
    os.system('mv *crj.fits raw_crj_tra_csv_files')
    os.system('mv *.csv raw_crj_tra_csv_files')
    os.system('mv *asn.fits raw_crj_tra_csv_files')


    
    
    
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

    create_asn_file(args.path)


