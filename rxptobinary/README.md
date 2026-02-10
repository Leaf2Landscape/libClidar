# README #


This is a tool to read Riegl's .rxp format and convert to either an ASCII .pts file or a custom binary format, for use in [voxeTLS](https://bitbucket.org/StevenHancock/voxelate). The code makes use of Riegl's RiVlib library.


## Dependencies ##

This software requires Riegl's [RiVlib](http://www.riegl.com/index.php?id=224) library. Choose the version for your gcc, either gcc44 (version 4) or gcc55 (version 5 or higher).

## Compilation ##

Set an environment variable to point to **RiVlib** and another, **CXXFLAGS**, to ensure that C++ 11 support is included. For example, in bash

    export RiVLib_DIR="$HOME/src/rivlib-2_5_11-x86_64-linux-gcc5"
    export CXXFLAGS="-std=c++11"


Then run

    cmake .
    make

## Usage ##

The command is **readRXP**. It reads a Riegl .rxp file. To align with other scans it requires a 4 by 4 matrix to put it to the plot coordinate system, known as an SOP matrix in RiScan. To but the plot in a world coordinate system it requires a second matrix, known as a POP matrix in RiScan. These matrices are both optional and if left out, the scan coordinate system will be used.

##### Input output filenames and format
    -input name;       input rxp filename
    -trans name;       plot transformation matrix filename
    -globMat name;     global transformation matrix if needed. Unity otherwise
    -output name;      output filename
    -ascii name;       ascii (.pts) output and filename
    -noBin;            don't output binary

##### Filtering options
    -maxZen zen;       max zenith angle to output, degrees
    -maxDev dev;       maximum return shape deviation to accept
    -dec n;            rate to decimate pointcloud by
    -bounds minX minY minZ maxX maxY maxZ;    bounds for ASCII output

There are a number of filtering options. 

**-maxDev**: Riegl scanners record a deviation parameter, which is a measure of the spreading of the return if it has hit multiple targets. This can lead to ghost points, spreading between two nearby surfaces. These ghost points can be removed by setting **maxDev** to a value. Generally a value between 10 and 100 is appropriate.

**-dec**: To reduce the point cloud density, the data can be decimated by only outputting one point per dec points.

**-bounds**: If there is a particular area of interest for a pts file, bounds can be set to only output points within those bounds.

**-maxZen**: The Riegl TLS instruments cannot see objects within 50 cm - 1.5 m (depending on setting used). Any beams which hit objects within that limit will give no return and be measured as a gap. This can cause bias in gap fraction based analyses. To avoid the false gaps created when using a tilt mount, a maximum zenith angle can be used. If the scanner is at 90 degrees on a tilt mount, a maximum zenith of 90 degrees is recommended.


## Output

The data is output as either a .pts file (point cloud ASCII file) or as a custom binary file for use in [voxeTLS](https://bitbucket.org/StevenHancock/voxelate), either in a flat binary of HDF5 format.

The pts file contains the columns:

     1 x:       x coordinate
     2 y:       y coordinate
     3 z:       z coordinate
     4 refl:    reflectance, as estimates bu the Riegl's own calibration
     5 zen:     zenith angle of laser shot
     6 az:      azimuth angle of laser shot
     7 range:   range to hit
     8 hitN:    hit number
     9 nHits:   total number of hits for this laser shot

The binary files contain the same information.



### Contribution guidelines ###

Please talk to svenhancock@gmail.com to suggest edits.


### License ###

Gnu Public License

### Who do I talk to? ###

svenhancock@gmail.com
                          

