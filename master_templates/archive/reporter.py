# import all essential libraries
import argparse
import glob
import numpy as np
import romspline as romSpline
import h5py
import os
from pycbc import pnutils
import lal
import hashlib
import pycbc.types as types
from pycbc.waveform import utils as wfutils
import subprocess


def format_h5_file_name(datafile, rad):
    pieces = datafile.split('/')
    method = pieces[-3].split(':')[0]
    simnum = pieces[-3].split(':')[1]
    res    = pieces[-2]
    radius = rad
    return 'CoRe_%s_%s_%s_%s.h5' %(method, simnum, res, radius)

def drill_down(dest, datafile, rad, eos, mass1, mass2):
    pieces = datafile.split('/')
    method = pieces[-3].split(':')[0]
    simnum = str(pieces[-3].split(':')[1])
    res    = str(pieces[-2].replace('R', ''))
    radius = str(rad.replace('r', ''))
    mass1  = str("%.2f" % mass1)
    mass2  = str("%.2f" % mass2)
    drillpath = (dest + '/eos_' + eos + '/mgrav_' + mass1 + '_' + mass2 + 
                 '/' + method + '_' + simnum + '/res_' + res + '/radius_' + rad)
    if not os.path.isdir(drillpath):
        os.makedirs(drillpath)
    return drillpath

# takes data from metadata.txt and formats it to be used in this conversion script
def format_metadata(name):
    data = ''
    with open(name, 'r') as file:
        # read a list of lines into data
        line = file.readline()
        while line:
            linea = ''
            lineb = ''
            if '= ' not in line:
                if line[0] != '#':
                    line = "# %s" %(line)
            else:
                line = line.replace('\'', '')
                line = line.replace('\"', '')
                line = line.replace('\n','')
                linea, lineb = line.split('= ')
                if ',' in lineb:
                    lineb = lineb.replace(' ', '')
                    lineb_new = lineb.split(',')
                    i = 0
                    while i < len(lineb_new):
                        try:
                            float(lineb_new[i])
                        except ValueError:
                            lineb_new[i] = "\'%s\'" %(lineb_new[i])
                        i += 1
                    lineb = ', '.join(lineb_new)
                elif (lineb.lower() == 'nan'):
                    lineb = "\'%s\'" %(lineb)
                else:
                    try:
                        float(lineb)
                    except ValueError:
                        lineb = "\'%s\'" %(lineb)
                line = "%s= %s\n" %(linea,lineb)
            data += line
            line = file.readline()
    return data

def write_metadata(metadata):
    data = format_metadata(metadata)
    with open(metadata, 'w') as file:
        file.writelines(data)
    file.close()

parser = argparse.ArgumentParser(
    description='This file creates porperly formatted HDF5 numerical relativity simulations for use with lalsuite/PyCBC.', 
    usage='python %(prog)s --src=[source path] --dest=[destination path]', 
    formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument(
    '--src', 
    help='Top level directory of directory tree to recursively look for files within.\n(default = [current working directory])', 
    metavar='[source directory]', 
    default=os.getcwd()
)
parser.add_argument(
    '--dest', 
    help='Top level directory of location to save newly formatted files.\nIf no directory with the specified name exists, one will be created.\n(default = [current working directory])\n***It is highly recommended that you specify a unique directory***', 
    metavar='[destination directory]', 
    default=os.getcwd()
)
parser.add_argument(
    '--type', 
    help='Type of file extension of source file.\nOptions are either text or h5\n(use "text" if your files are in .dat/.txt format)', 
    metavar='[text/h5]', 
    required=True
)
args = parser.parse_args()

# datafiles = glob.glob('/home/dwhite/GWPAC/h5_conversions/original_files/CoRe/Public/BAM:0001/R01/data.h5')
datafiles = []
path = args.src

for (dirpath, dirnames, filenames) in os.walk(path):
    for file in filenames:
#         if (('.txt' in file) or ('.dat' in file)) and ('metadata.txt' not in file):
#             datafiles.append(os.path.join(dirpath, file))
        if args.type == 'h5':
            if '.h5' in file:
                datafiles.append(os.path.join(dirpath, file))
        elif args.type == 'text':
            if ('.txt' in file) and ('metadata' not in file):
                datafiles.append(os.path.join(dirpath, file))
            if ('.dat' in file):
                datafiles.append(os.path.join(dirpath, file))

# go through list and convert each file as necessary into its own .h5
for datafile in datafiles:
    filepath = os.path.dirname(datafile)
    metadata = filepath + '/metadata.txt'
    
    # format metadata file for our needs and overwrite the old version
    write_metadata(metadata)      
    
    # read metadata.txt in as code, so as to be able to call in its stored values
    # as if they were variables
    exec(open(metadata, 'r'))
    
    # create some default variables based on metadata.txt
    # easier to do this here than to change the variable name each other time it's used
    grav_mass1      = max(id_mass_starA, id_mass_starB)
    grav_mass2      = min(id_mass_starA, id_mass_starB)
    total_grav_mass = id_mass
    spins1          = id_spin_starA
    spins2          = id_spin_starB
    eos             = id_eos
    eccentricity    = id_eccentricity
    spin1x          = id_spin_starA[0]
    spin1y          = id_spin_starA[1]
    spin1z          = id_spin_starA[2]
    spin2x          = id_spin_starB[0]
    spin2y          = id_spin_starB[1]
    spin2z          = id_spin_starB[2]
    
    h = h5py.File(datafile, 'r')
    h_group = ''
    keys = []
    radii = []
    lmax = 0
    if 'rh_22' in h:
        h_group = h['rh_22']
        keys = h_group.keys()
    
    for key in keys:
        pieces = key.split('_')
        radius = pieces[-1].replace('.txt','')
        if (int(str(pieces[1])[-1]) > lmax):
            lmax = int(str(pieces[1])[-1])
        if radius not in radii:
            radii.append(radius)
            
    for rad in radii:
        # naming convention for the newly converted file (may change per use case)
        h5file = format_h5_file_name(datafile, rad)
        drill = drill_down(args.dest, datafile, rad, eos, grav_mass1, grav_mass2)
            
        if ((datafile.split('/'))[-2] == 'R01'):
            if ((rad == radii[-1]) and ('Inf' not in radii[-1])):
                os.chdir(drill)
                comline = "waverep --gwfile " + h5file + " energy overview spectra strain"
                subprocess.call(comline, shell=True)
            elif ((rad == radii[-2]) and ('Inf' in radii[-1])):
                os.chdir(drill)
                comline = "waverep --gwfile " + h5file + " energy overview spectra strain"
                subprocess.call(comline, shell=True)

print 'All done!'