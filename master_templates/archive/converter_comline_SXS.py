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

# unit conversion functions

# converts time to units of Mtotal, based on what format it was already saved in...
# ...as defined by variable "timeFormat"
def convert_time(timeFormat, times, total_grav_mass):
    if (timeFormat == "s"):
        for time in times:
            time /= lal.MTSUN.SI
            timeFormat = "Msun"
    if (timeFormat == "Msun"):
        for time in times:
            time /= total_grav_mass
            timeFormat = "Mtotal"
    if (timeFormat == "Mtotal"):
        return times
    else:
        print "Time is in an unrecognized format. Please edit the metadata file and try again."
        sys.exit()

# converts strain to units of rhOverM, based on what format it was already saved in...
# ...as defined by variable "strainFormat"
def convert_strain(strainFormat, strainVal1, strainVal2, total_grav_mass):
    if (strainFormat == "rh"):
        for i in range(0, len(strainVal1)):
            strainVal1[i] /= total_grav_mass
        # ONLY if in PlusCross
        if (dataFormat == "PlusCross"):
            for i in range(0, len(strainVal2)):
                strainVal2[i] /= total_grav_mass
        return strainVal1, strainVal2
    elif (strainFormat == "rhOverM"):
        return strainVal1, strainVal2
    else:
        print "Strain is in an unrecognized format. Please edit the metadata file and try again."
        sys.exit()

# converts masses to units of Msun, based on what format they were already saved in...
# ...as defined by variable "massFormat"        
def convert_mass(massFormat, grav_mass1, grav_mass2, total_grav_mass):
    if (massFormat == "kg"):
        print "Old mass:"
        print mass1
        grav_mass1 /= lal.MSUN_SI
        grav_mass2 /= lal.MSUN_SI
        massFormat = "Msun"
        print "New mass:"
        print mass1
    elif (massFormat == "Mtotal"):
        grav_mass1 *= total_grav_mass
        grav_mass2 *= total_grav_mass
        massFormat = "Msun"
    if (massFormat == "Msun"):
        return grav_mass1, grav_mass2
    else:
        print "Mass1/mass2 is in an unrecognized format. Please edit the metadata file and try again."
        sys.exit()
        
# # convert distance to Mpc
# def convert_distance(distFormat, distance):
#     if (distFormat == "Msun"):
#         distance *= lal.MRSUN_SI
#         distFormat = "m"
#     elif (distFormat == "km"):
#         distance *= 1000
#         distFormat = "m"
#     if (distFormat == "m"):
#         distance /= lal.PC_SI
#         distFormat = "pc"
#     if (distFormat == "pc"):
#         distance /= 10e6
#         distFormat = "Mpc"
#     if (distFormat == "Mpc"):
#         return distance
#     else:
#         print "Distance is in an unrecognized format. Please edit the metadata file and try again."
#         sys.exit()

# takes data from metadata.txt and formats it to be used as if it were code
# useful for formatting before trying to use exec(read())
# works with all metadata files
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
                if '-' in linea:
                    linea.replace('-','_')
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

#creates new metadata file once it assures it is formatted correctly
def write_metadata(metadata):
    data = format_metadata(metadata)
    with open(metadata, 'w') as file:
        file.writelines(data)
    file.close()

# simply returns the name of the datafile, removing the file path and just leaving the file name
def format_h5_file_name(datafile):
    pieces = datafile.split('/')
#     method = pieces[-3].split(':')[0]
#     simnum = pieces[-3].split(':')[1]
#     res    = pieces[-2]
#     radius = rad
    return pieces[-1]

# returns the directory defined by "dest". this will be where the new file is written and saved
def drill_down(dest):
#     pieces = datafile.split('/')
#     method = pieces[-3].split(':')[0]
#     simnum = str(pieces[-3].split(':')[1])
#     res    = str(pieces[-2].replace('R', ''))
#     radius = str(rad.replace('r', ''))
#     mass1  = str("%.2f" % mass1)
#     mass2  = str("%.2f" % mass2)
#     drillpath = (dest + '/eos_' + eos + '/mgrav_' + mass1 + '_' + mass2 + 
#                  '/' + method + '_' + simnum + '/res_' + res + '/radius_' + rad)
#     if not os.path.isdir(drillpath):
#         os.makedirs(drillpath)
    return dest
    
# necessary variables for conversion to put things in correct units
# sometimes these are saved in metadata. they are here only as a backup...
# ...if not already in metadata
massFormat   = "Mtotal"
timeFormat   = "Mtotal"
strainFormat = "rhOverM"
dataFormat   = "PlusCross"
delta_t      = 1./(4096*16)

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
    eccentricity    = relaxed_eccentricity
    spin1x          = 0
    spin1y          = 0
    spin1z          = 0
    spin2x          = 0
    spin2y          = 0
    spin2z          = 0
    lmax            = 2
                                    
    h = h5py.File(datafile, 'r')
    h_group = h['Extrapolated_N4.dir']
    keys = h_group.keys()
            
    # naming convention for the newly converted file (may change per use case)
    h5file = format_h5_file_name(datafile)
    drill = drill_down(args.dest)

    # finite or infinite? used for 'NR-techniques' attribute
    waveRad = 'Extrapolated-Waveform'

    with h5py.File(drill + '/' + h5file,'w') as fd:
        print 'creating file %s from %s...' %(h5file, datafile)

        # store metadata.txt in our new .h5
        aux = fd.create_group('auxiliary-info')
        mdata = open(metadata, 'r')
        aux.create_dataset('metadata.txt', data=mdata.read())
        mdata.close()

#             # create a readme that explains the file naming convention
#             # e.g. 'CoRe_BAM_SLy_150_150_nospin_R01_r00550.h5'
#             nc = open('/home/dwhite/GWPAC/h5_conversions/original_files/CoRe/h5_naming_convention.txt', 'r')
#             aux.create_dataset('h5_naming_convention.txt', data=nc.read())
#             nc.close()

        # check the spin values to give the correct simulation type
        # also used for "NS-spins-meaningful" .h5 attribute
        simtype = 'non-spinning'
        spinval = True
        if ((spin1x != 0) or (spin1y != 0) 
                        or (spin2x != 0) or (spin2y != 0)):
            simtype = 'precessing'
        elif ((spin1z != 0) or (spin2z != 0)):
            simtype = 'aligned-spin'
        else:
            spinval = False

        # calculate the value for mean_anomaly (assuming this hasn't been calculated)
        # 0 is default for low eccentricity values (<= .001)
        mean_anom = 0
        if (eccentricity > .001):
            # default value for when it hasn't been properly calculated
            mean_anom = -1

        # create all relevant attributes for converted .h5 file

        # default Format 1 attributes
        mchirp, eta = pnutils.mass1_mass2_to_mchirp_eta(initial_mass1, initial_mass2)
        hashtag = hashlib.md5()
        fd.attrs.create('Format', 1)
        fd.attrs.create('type', 'BBH')
        fd.attrs.create('name',simulation_name)
        fd.attrs.create('alternative-names', alternative_names)
        fd.attrs.create('NR-group', 'SXS')
        fd.attrs.create('NR-code', 'SpEC')
        fd.attrs.create('modification-date', '')
        fd.attrs.create('point-of-contact-email', '')
        fd.attrs.create('simulation-type', '')
        fd.attrs.create('INSPIRE-bibtex-keys', '')
        fd.attrs.create('license', 'public')
        fd.attrs.create('Lmax', lmax)
        fd.attrs.create('files-in-error-series', '')
        fd.attrs.create('comparable-simulation', '')
        fd.attrs.create('production-run', 1)
        fd.attrs.create('object1', 'BH')
        fd.attrs.create('object2', 'BH')
        fd.attrs.create('mass1', initial_mass1)
        fd.attrs.create('mass2', initial_mass2)
        fd.attrs.create('eta', eta)
        fd.attrs.create('f_lower_at_1MSUN', initial_orbital_frequency / (lal.TWOPI * lal.MTSUN_SI))
        fd.attrs.create('spin1x', spin1x)
        fd.attrs.create('spin1y', spin1y)
        fd.attrs.create('spin1z', spin1z)
        fd.attrs.create('spin2x', spin2x)
        fd.attrs.create('spin2y', spin2y)
        fd.attrs.create('spin2z', spin2z)
        # HARDCODING for non-spinning / aligned-spin
        # this, too, could one day be in metadata.txt, if we found it worthy
        fd.attrs.create('LNhatx', 0.0)
        fd.attrs.create('LNhaty', 0.0)
        fd.attrs.create('LNhatz', 1.0)
        fd.attrs.create('nhatx', 1.0)
        fd.attrs.create('nhaty', 0.0)
        fd.attrs.create('nhatz', 0.0)
        fd.attrs.create('Omega', initial_orbital_frequency)
        fd.attrs.create('eccentricity', eccentricity)
        fd.attrs.create('mean-anomaly', '')

        # hashtag stuff
        # no idea if this is essential, so I'm just leaving it here
        fd.attrs.create('hashtag', hashtag.digest())
        hashtag.update(fd.attrs['type'])

        tempkey = h_group[keys[1]][:,0]
        dmax = 0
        maxloc = 0
        strainMag = [0]*(len(tempkey)+10)

        # find the largest amplitude and set its time stamp to t=0; adjust all others
        #    this should be a resonably accurate way of finding the moment of "merger" 
        #    and making it t=0, in keeping with PyCBC conventions

        for key in keys:
            if 'History' not in key:
                timeval = h_group[key][:,0]
                strain1 = h_group[key][:,1]
                strain2 = h_group[key][:,2]
                for i in range(0, len(timeval)):
                    strainMag[i] += np.sqrt((strain1[i])**2 + (strain2[i])**2)
                    if strainMag[i] > dmax:
                        maxloc = i
                        dmax = strainMag[i]

        timeAdjust = tempkey[maxloc]

        for key in keys:
            if 'l2_m2' in key:
                chopchop = key.split('_')
                ell = chopchop[1]
                em = chopchop[2].replace('.dat','')

                # pull time and strain data from appropriate file(s) and run appropriate conversions
                times      = h_group[key][:,0]
                strainVal1 = h_group[key][:,1]
                strainVal2 = h_group[key][:,2]

                for i in range(0, len(times)):
                    times[i] = (times[i] - timeAdjust)

#                 strainVal1, strainVal2 = convert_strain(strainFormat, strainVal1, strainVal2, total_grav_mass)
#                 times = convert_time(timeFormat, times, total_grav_mass)

                times       = np.array(types.TimeSeries(times, delta_t=delta_t))
                strainVal1  = types.TimeSeries(strainVal1, delta_t=delta_t)
                strainVal2  = types.TimeSeries(strainVal2, delta_t=delta_t)

                strainAmp   = []
                strainPhase = []

                # run romSpline to convert into reduced order spline, then assign final .h5 values
                # and write all data to .h5 file
                # handled independently for Magnitude/Argument vs. Pluss/Cross data, based on 
                # unique needs for each format
                if (dataFormat == "MagArg"):
                    strainAmp = np.array(strainVal1)
                    strainPhase = np.array(strainVal2)
                elif (dataFormat == "PlusCross"):
                    strainAmp = wfutils.amplitude_from_polarizations(strainVal1, strainVal2).data
                    strainPhase = wfutils.phase_from_polarizations(strainVal1, strainVal2).data 
                else:
                    print "dataFormat is incorrect or is not specified. Edit the metadata file and try again."

                print '    key = %s' %(key)
                print '        fitting spline...'

                try:
                    # when a mode has nothing but zeros for strain, we don't want to add it to the .h5
                    sAmpH = romSpline.ReducedOrderSpline(times, strainAmp, rel=True, verbose=False)
                    sPhaseH = romSpline.ReducedOrderSpline(times, strainPhase, rel=True, verbose=False)

                    grAmp = fd.create_group('amp_' + ell + '_' + em)
                    grPhase = fd.create_group('phase_' + ell + '_' + em)

                    sAmpH.write(grAmp)
                    sPhaseH.write(grPhase)

                    print '        spline created.'
                except AssertionError:
                    print '        SPLINE SKIPPED (no strain data).'

        print 'file created.\n'

print 'All done!'