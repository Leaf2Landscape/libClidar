# Voxelate #

This is a program to make voxel maps from either full-waveform ALS or TLS data. The algorithms is described in an open access paper:

[Hancock, S., Anderson, K., Disney, M. and Gaston, K.J., 2017. Measurement of fine-spatial-resolution 3D vegetation structure with airborne waveform lidar: Calibration and validation with voxelised terrestrial lidar. Remote Sensing of Environment, 188, pp.37-50.](http://www.sciencedirect.com/science/article/pii/S0034425716304205)

This directory contains the programs to generate voxels from either ALS or TLS data. For comparisons between TLS and ALS use the programs in a TLS to ALS voxel [library](https://bitbucket.org/StevenHancock/voxel_lidar).


# Tools #

The library contains a number of tools for voxelising either ALS or TLS data. Note that the ALS data assumes full-waveform, and the TLS data must first be processed to form a polar projected gap-fraction file (a point cloud that also includes the misses). The programs are:

**voxelate**: generates a voxel map from full-waveform ALS data

**voxelTLS**: generates a voxel map from TLS data, one scan at a time

**combinePAI.py**: Combines voxel maps from multiple TLS scans to produce a single map


### Dependencies

All programs depend on these libraries:

* Gnu Scientific Library
* HDF5
* libgeotiff
* [Minpack (Levenberg-Maruqardt)](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html)


They also require two libraries without package managers, which must be cloned:

* [C-tools](https://bitbucket.org/StevenHancock/tools)
* [libClidar](https://bitbucket.org/StevenHancock/libclidar)


Clone these libraries and point to their locations with the environment variables:

    HANCOCKTOOLS_ROOT
    LIBCLIDAR_ROOT


To compile, type:

  **make THIS=voxelTLS**

Make sure that ~/bin/$ARCH exists and is in the path, where $ARCH is the result of `uname -m`.

  **make THIS=voxelTLS install**

The python script, **combinePAI.py** depends upon the following packages.

* numpy
* argparse
* sys
* csv
* h5py


## Voxelate ##
voxelate.c generates a voxel map from full-waveform ALS data. This requires las 1.3 or higher data and knowledge of the pulse width and noise parameters. In the above paper these instrument specific parameters were found by calibrating against TLS data, using functions in this [library](https://bitbucket.org/StevenHancock/voxel_lidar).

## voxelTLS ##

voxelTLS.c generates a voxel map from TLS data. The input must be in a custom binary format. This can be be generated from Riegl's rxp format usiong this [tool](https://bitbucket.org/StevenHancock/rxptobinary/overview), or from [SALCA](https://www.sciencedirect.com/science/article/pii/S0168192314001749) using this [tool](https://bitbucket.org/StevenHancock/salcamake/src/default/).


### usage ###

The code expects an input filename (custom binary polar format, which it is hoped to replace soon with a more useful format), a voxel extent, voxel resolution and an output filename. The output will either be in an ASCII file (**warning**, large file size) or a compressed HDF5 format. An example command would be:

    *voxelTLS -input input.bin -output output.hdf -hdf -bounds -10 -10 0 10 10 30 -res 0.75*

The command line options are given below:

    -input name;     input TLS filename. rxp or ptx
    -inList list;    input file list for multiple files
    -output name;    output filename
    -hdf;    write output as a compressed HDF5 file
    -bounds minX minY minZ maxX maxY maxZ:    voxel bounds
    -res res;        voxel resolution in metrs
    -vRes x y z;     voxel resolution in each dimension
    -coord x y z;    do a single beam at this origin
    -grad x y z;     beam vector
    -ang zen az;     beam angle in degrees
    -pSigma sigma;   pulse width
    -fSigma sigma;   footprint width
    -gauss;          Gaussian footprint
    -voxGap;         just do the voxel gap fractions
    -readBounds;     read bounds from data
    -maxRefl x;      DN for reflectance=1
    -beamDiv ang;    beam divergence, in degrees
    -correctR;       correct reflectance for range in gap fraction
    -minGap gap;     minimum trustable gaop fraction when scaling
    -appRefl rho;    scalar between TLS and ALS reflectance
    -printFilled;    only write out filled voxel coordinates

    # DTM
    -dtm name;       use a DTM to remove ground hits
    -demTol tol;     tolerance for counting points as ground from DTM


The results can be output in either ASCII or HDF5 format. If ASCII, the file contains the following variables in each column.

    1 x;            x coordinate of the centre of the voxel
    2 y;            y coordinate of the centre of the voxel
    3 z;            z coordinate of the centre of the voxel
    4 gap0;         gap fraction within this voxel from scan 0. This is inMiss/(inHit+inMiss)
    5 inHit0;       number of beams hitting objects within the voxel
    6 inMiss0;      number of beams that passed through the voxel without hitting targers
    7 gapTo0;       gap fraction from scan 0 to this voxel. This is hits/(miss+hits)
    8 hits0;        number of beams that reached this voxel
    9 miss0;        number of beams that were blocked before this voxel
    10 PAIbel0;     PAI from Beland et al (2011)
    11 volSamp0;    volume of voxel sampled by unobstructed beams
    12 volTrans0;   volume of voxel that would have been sampled by beams if the voxel were empty
    13 PAIrad0;     PAI from radius method in Hancock et al (2017)
    14 meanRefl0;   mean reflectance returned by all beams reaching this voxel
    15 meanZen0;    mean zenith angle of all beams reaching this voxel

Columns 4-15 will repeat if multiple scans have been processed at the same time.

HDF file output contains the same information in arrays, although the coordinates (x,y,z) must be inferred from the known resolution.

Negative numbers indicate missing data. A gapTo of *-2* means that that voxel was undeground relative to the TLS, and so indicates a voxel with no usable TLS samples. A gapTo of *0* means that no TLS beams managed to reach that voxel. If the gapTo is 0, all other values for that voxel will be set to the missing data flag of -1.


**References**

[Béland, M., Widlowski, J.L., Fournier, R.A., Côté, J.F. and Verstraete, M.M., 2011. Estimating leaf area distribution in savanna trees from terrestrial LiDAR measurements. Agricultural and Forest Meteorology, 151(9), pp.1252-1266.](https://www.sciencedirect.com/science/article/pii/S0168192311001481)

[Hancock, S., Anderson, K., Disney, M. and Gaston, K.J., 2017. Measurement of fine-spatial-resolution 3D vegetation structure with airborne waveform lidar: Calibration and validation with voxelised terrestrial lidar. Remote Sensing of Environment, 188, pp.37-50.](http://www.sciencedirect.com/science/article/pii/S0034425716304205)


## combinePAI.py ##
combinePAI.py reads a list of output files from **voxelTLS** and combines them in order to reduce the impact of attenuation and variable TLS sampling density with distance from scan centres. It calculates the mean gap fraction and PAIvalues for each voxel, weighted by an estimate of the uncertainty of each scan's estimate. This measure of uncertainty used for the weighting is set as either one minus the gap fraction from the scan location to that voxel (all values except for PAI from Beland), or one minus the fraction of volume sampled by the TLS beams that reached the voxel (PAI from Beland).

### Usage ###

Command line options are:

    -h, --help          show this help message and exit
    --input INNAMEN     Input of combined scans No default
    --inList LISTNAMEN  Input list of single scans No default
    --output OUTNAMEN   Outpt filename Default = meanVox.vox
    --minGap MINGAP     Minimum trustable gap fraction to a voxel Default = 0
    --hdf               Write output in HDF5 Default si ASCII


To run, type something like:

    python3 $bin/combinePAI.py --inList list.txt --output combinedVoxels.txt

Where *bin* is the directory containing combinePAI.py, *list.txt* is an ASCII file containing a list of absolute paths to the filenames output by voxelTLS and *combinedVoxels.txt* is the output filename.


The output format is similar to that of voxelTLS, except that there is now only a single gap fraction and PAI per voxel per method, rather than a value per voxel per scan per method. The gapTo metric (measure of uncertainty) is replaced by the mean gap fraction to that voxel across all scans (mGapTo) and the maximum gap fraction to that voxel from any one scans (maxGapTo). The latter (maxGapTo) and the total number of beams intersecting (totBeam) are perhaps the most useful measure of uncertainty for the combined voxel map.

The output variables are:

    1 x:         x coordinate of that voxel
    2 y:         y coordinate of that voxel
    3 z:         z coordinate of that voxel
    4 mGapTo:    mean gap fraction to that voxel from all scans
    5 maxGapTo:  maximum gap fraction to that voxel for single scan
    6 mGapIn:    mean gap fraction within that voxel across all scans
    7 mPAIb:     mean PAI according to Beland et al (2011)
    8 mPAIrad:   mean PAI using the radius method
    9 mPAIw:     not yet working
    10 totBeam:  total number of TLS beams reaching that voxel
    11 bPAIb:    PAI from the single scan with the best view of that voxel, as used by Beland et al (2011)


### Who do I talk to? ###

svenhancock@gmail.com
