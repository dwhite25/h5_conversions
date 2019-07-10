# ----------------------------------------------------------------
# if you need to run things in command line rather than in Jupyter
# (like, for instance, if you need to submit a job),
# this is the file to use. Specific to CoRe formatting, of course.
# ----------------------------------------------------------------

# import all essential libraries
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
    elif (massFormat == "mTotal"):
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

def write_metadata(metadata):
    data = format_metadata(metadata)
    with open(metadata, 'w') as file:
        file.writelines(data)
    file.close()

def format_h5_file_name(datafile, rad):
    pieces = datafile.split('/')
    method = pieces[-3].split('_')[0]
    simnum = pieces[-3].split('_')[1]
    res    = pieces[-2]
    radius = rad
        
    return 'CoRe_%s%s_%s_%s.h5' %(method, simnum, res, radius)

# necessary variables for conversion to put things in correct units
# sometimes these are saved in metadata. they are here only as a backup...
# ...if not already in metadata
massFormat   = "Msun"
timeFormat   = "Mtotal"
strainFormat = "rhOverM"
dataFormat   = "PlusCross"
delta_t      = 1./(4096*16)

datafiles = []
# path = os.getcwd()
path ="/home/dwhite/GWPAC/h5_conversions/original_files/CoRe/Public/"

for (dirpath, dirnames, filenames) in os.walk(path):
    for file in filenames:
#         if (('.txt' in file) or ('.dat' in file)) and ('metadata.txt' not in file):
#             datafiles.append(os.path.join(dirpath, file))
        if '.h5' in file:
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
        
        # finite or infinite? used for 'NR-techniques' attribute
        waveRad = 'Finite-Radius-Waveform'
        if 'Inf' in rad:
            waveRad = 'Extrapolated-Waveform'

        with h5py.File('/home/dwhite/GWPAC/h5_conversions/converted_files/CoRe/' + h5file,'w') as fd:
            print 'creating file %s from %s...' %(h5file, datafile)
            
            # run mass conversion function for correct PyCBC conventions
            grav_mass1, grav_mass2 = convert_mass(massFormat, grav_mass1, grav_mass2, total_grav_mass)

            # store metadata.txt in our new .h5
            aux = fd.create_group('auxiliary-info')
            mdata = open(metadata, 'r')
            aux.create_dataset('metadata.txt', data=mdata.read())
            mdata.close()

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
            if (id_eccentricity > .001):
                # default value for when it hasn't been properly calculated
                mean_anom = -1

            # create all relevant attributes for converted .h5 file
            
            # default Format 1 attributes
            mchirp, eta = pnutils.mass1_mass2_to_mchirp_eta(id_mass_starA, id_mass_starB)
            hashtag = hashlib.md5()
            fd.attrs.create('Format', 1)
            fd.attrs.create('type', 'BNS')
            fd.attrs.create('name','CoRe:' + database_key)
            fd.attrs.create('alternative-names', simulation_name)
            fd.attrs.create('NR-group', 'CoRe')
            fd.attrs.create('NR-code', id_code + ', ' + evolution_code)
            fd.attrs.create('modification-date', '')
            fd.attrs.create('point-of-contact-email', 'computational.relativity@gmail.com')
            fd.attrs.create('simulation-type', simtype)
            fd.attrs.create('INSPIRE-bibtex-keys', reference_bibkeys)
            fd.attrs.create('license', 'public')
            fd.attrs.create('Lmax', 2)  ## CoRe only has 2-2 mode
            fd.attrs.create('files-in-error-series', '')
            fd.attrs.create('comparable-simulation', '')
            fd.attrs.create('production-run', 1)
            fd.attrs.create('object1', 'NS')
            fd.attrs.create('object2', 'NS')
            fd.attrs.create('mass1', grav_mass1/total_grav_mass)
            fd.attrs.create('mass2', grav_mass2/total_grav_mass)
            fd.attrs.create('eta', eta)
            fd.attrs.create('f_lower_at_1MSUN', id_gw_frequency_Hz)
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
            fd.attrs.create('Omega', id_gw_frequency_Momega22/2.0)
            fd.attrs.create('eccentricity', eccentricity)
            fd.attrs.create('mean_anomaly', mean_anom)
            fd.attrs.create('NR-techniques', ('Quasi-Equilibrium-ID, ' 
                            + metric_scheme + ', Psi4-integrated, ' + waveRad))
            
            # hashtag stuff
            # no idea if this is essential, so I'm just leaving it here
            fd.attrs.create('hashtag', hashtag.digest())
            hashtag.update(fd.attrs['type'])
            
            # attributes unique to NS sims
            fd.attrs.create('file-format-version', 2)
            fd.attrs.create('mass1-msol', grav_mass1)
            fd.attrs.create('mass2-msol', grav_mass2)
            fd.attrs.create('baryonic-mass1-msol', id_rest_mass_starA)
            fd.attrs.create('baryonic-mass2-msol', id_rest_mass_starB)
            fd.attrs.create('NS-spins-meaningful', spinval)
            fd.attrs.create('EOS-name', id_eos)
            fd.attrs.create('EOS-references', 'http://computational-relativity.org/EOS.html')
            fd.attrs.create('EOS-remarks', '')
            fd.attrs.create('have-ns-tidal-lambda', True)
            fd.attrs.create('tidal-lambda1-l2', id_Lambdaell_starA[0])
            fd.attrs.create('tidal-lambda1-l3', id_Lambdaell_starA[1])
            fd.attrs.create('tidal-lambda1-l4', id_Lambdaell_starA[2])
            fd.attrs.create('tidal-lambda2-l2', id_Lambdaell_starB[0])
            fd.attrs.create('tidal-lambda2-l3', id_Lambdaell_starB[1])
            fd.attrs.create('tidal-lambda2-l4', id_Lambdaell_starB[2])
            fd.attrs.create('chirp-mass', ((grav_mass1 * grav_mass2)**(.6))/((grav_mass1 + grav_mass2)**(.2)))

            # all other attributes from original metadata
            # (stuff that didn't otherwise have somewhere to go)
            fd.attrs.create('id_code',id_code)
            fd.attrs.create('id_type',id_type)
            fd.attrs.create('id_mass',id_mass)
            fd.attrs.create('id_rest_mass',id_rest_mass)
            fd.attrs.create('id_mass_ratio',id_mass_ratio)
            fd.attrs.create('id_ADM_mass',id_ADM_mass)
            fd.attrs.create('id_ADM_angularmomentum',id_ADM_angularmomentum)
            fd.attrs.create('id_gw_frequency_Momega22',id_gw_frequency_Momega22)
            fd.attrs.create('id_kappa2T',id_kappa2T)
            fd.attrs.create('id_LoveNum_l2_starA', id_LoveNum_kell_starA[0])
            fd.attrs.create('id_LoveNum_l3_starA', id_LoveNum_kell_starA[1])
            fd.attrs.create('id_LoveNum_l4_starA', id_LoveNum_kell_starA[2])
            fd.attrs.create('id_LoveNum_l2_starB', id_LoveNum_kell_starA[0])
            fd.attrs.create('id_LoveNum_l3_starB', id_LoveNum_kell_starA[1])
            fd.attrs.create('id_LoveNum_l4_starB', id_LoveNum_kell_starA[2])
            fd.attrs.create('evolution_code', evolution_code)
            fd.attrs.create('grid_refinement_levels', grid_refinement_levels)
            fd.attrs.create('grid_refinement_levels_moving', grid_refinement_levels_moving)
            fd.attrs.create('grid_refinement_levels_npoints', grid_refinement_levels_npoints)
            fd.attrs.create('grid_refinement_levels_moving_npoints', grid_refinement_levels_moving_npoints)
            fd.attrs.create('grid_spacing_min', grid_spacing_min)
            fd.attrs.create('grid_symmetries', grid_symmetries)
            fd.attrs.create('grid_shells', grid_refinement_levels)
            fd.attrs.create('grid_shells_radial_npoints', grid_shells_radial_npoints)
            fd.attrs.create('grid_shells_angular_npoints', grid_shells_angular_npoints)
            fd.attrs.create('grid_conservative_amr', grid_conservative_amr)
            fd.attrs.create('metric_scheme', metric_scheme)
#             fd.attrs.create('metric_boundary_conditions', metric_boundary_conditions)
            fd.attrs.create('hydro_flux', hydro_flux)
            fd.attrs.create('hydro_reconstruction', hydro_reconstruction)
            fd.attrs.create('hydro_atmosphere_level', hydro_atmosphere_level)
            fd.attrs.create('hydro_atmosphere_factor', hydro_atmosphere_factor)
            fd.attrs.create('number_of_orbits', number_of_orbits)
            fd.attrs.create('evolution_mol_scheme', evolution_mol_scheme)
            fd.attrs.create('eos_evolution_Gamma_thermal', eos_evolution_Gamma_thermal)

            tempkey = h_group[keys[0]][:,0]
            dmax = 0
            maxloc = 0
            strainMag = [0]*(len(tempkey)+10)

            # find the largest amplitude and set its time stamp to t=0; adjust all others
            #    this should be a resonably accurate way of finding the moment of "merger" 
            #    and making it t=0, in keeping with PyCBC conventions

            for key in keys:
                if rad in key:
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
                if rad in key:
                    pieces = key.split('_')
                    # EX: Rh_l2_m0_r00400.txt

                    # pull time and strain data from appropriate file(s) and run appropriate conversions
                    times      = h_group[key][:,0]
                    strainVal1 = h_group[key][:,1]
                    strainVal2 = h_group[key][:,2]

                    for i in range(0, len(times)):
                        times[i] = (times[i] - timeAdjust)

                    strainVal1, strainVal2 = convert_strain(strainFormat, strainVal1, strainVal2, total_grav_mass)
                    times = convert_time(timeFormat, times, total_grav_mass)

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

                        grAmp = fd.create_group('amp_' + pieces[1] + '_' + pieces[2])
                        grPhase = fd.create_group('phase_' + pieces[1] + '_' + pieces[2])

                        sAmpH.write(grAmp)
                        sPhaseH.write(grPhase)

                        print '        spline created.'
                    except AssertionError:
                        print '        SPLINE SKIPPED (no strain data).'

            print 'file created.\n'

print 'All done!'