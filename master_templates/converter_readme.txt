# GW_h5_converter
Converts existing gravitational wave simulation data into the HDF5 format to be used with PyCBC.

To use, place converter.py and a properly formatted version of metadata.txt in the same folder as any gravitational wave data files you wish to convert, then run converter.py. It will convert all files within the folder that have either .txt or .dat as its expension, but all files must use the same metadata.

Gravitational wave data files must have a minimum of three columns of data: column 1 must be timestamps, and columns two and three must be either the plus and cross values for the magnitude of strain, or else column two must be magnitude and column three must be argument (angle).

The process implemented for converting gravitational wave simulations into HDF5 files is heavily dependent upon the proper formatting of a "metadata.txt" file, which stores important information about the gravitational wave to be converted. Metadata.txt is imported and turned into code when coverter.py is run, so its careful formatting is crucial. A sample metadata.txt file is attached within this repository; further details regarding formatting will be discussed in this README file. Using the two in conjunction for reference should lead to a properly formatted file. 

## Style Guide

Since the metadata file will be read in as code, the file should contain data in the standard variable/value format. 
```
  Ex: 
    variableText   = "text"
    variableInt    = 11
    variableFloat  = 1.0
```    
Anything that is included in the metadata file but is not intended to be read into converter.py as code (in other words, anything useful to have but that is not inlcluded by default should be commented out with a pound (#), per Python convention
```
  Ex:
    # this line is commented out
    this line is not
```
## Metadata Critical Fields
The following fields must exist in the metadata file, and must have a meaningful value from the list of possible values for each:

#### dataFormat = "PlusCross"
This field explains what format strain data was saved in. **Possible values:** 
```
PlusCross - data columns two and three in original data file represent 
            two strain values, along opposite directions
MagArg    - data column two in original file represents strain magnitude
            (amplitude); column three represents argument (angle)
```
#### timeFormat = "Msun"
Indicates what ocnvention of units was used to represent time. **Possible values:**
```
Msun    - Solar units. 1 Msun = the mass of one Sol
Mtotal  - Units as a fraction of mass of binary system, in units of Msun
          Ex: If mass of binary system is 2.5 Sols, then 1 Msun in time would be .4 Mtotal
s       - Seconds.
```
#### strainFormat = "rhOverM"
This field idenitfies the units that were used to express the magnitude of strain values. If dataFormat is PlusCross and strainFormat is anything other than rhOverM, columns 1 and 2 of data file will be converted to rhOverM. If dataFormat is MagArg, only column 2 will be converted. **Possible Values:**
```
rhOverM - Strain value, multiplied by observer distance (I believe?) and divided by system mass
rh      - Strain value, multiplied by observer distance (I think)
```
#### massFormat = "Msun"
Indicates which units were used to save numerical data for the mass of each stellar body. **Possible values:**
```
Msun     - solar units. 1 Msun is the mass of one Sol
Mtotal   - Units as a fraction of mass of binary system, in units of Msun
           Ex: If mass of binary system is 2.5 Sols, then 1 Msun in mass would be .4 Mtotal
```
#### grav_mass1 = 1.528
Mass of object 1, in whatever units are defined by massFormat

#### grav_mass2 = 1.222
Mass of object 2, in whatever units are defined by massFormat

#### total_grav_mass = 2.75
Total mass of the system, i.e. the sum of grav_mass1 and grav_mass2. This field is only necessary if any data was saved in the Mtotal format (if in any other format, the total mass of the system will be calculated automatically by adding grav_mass1 and 2) 

#### l = 2
#### m = 2

## Metadata Essential Fields
These fields must be included in the metadata.txt file, but do not directly play a role in the proper formatting of the file. If necessary, values of 0 can be falsely inserted into these fields in the absense of real data
```
# spinning values per stellar body, broken into x/y/z vectors
spin1x = 0.0 
spin1y = 0.0 
spin1z = 0.0 
spin2x = 0.0 
spin2y = 0.0 
spin2z = 0.0

# alternative mass formats
ADM_mass1 = 0
ADM_mass2 = 0
total_ADM_mass = 0
Baryon_mass1 = 0
Baryon_mass2 = 0
total_Baryon_mass = 0

# tidal paramer per stellar body
lambda_1 = 0
lambda_2 = 0
```
