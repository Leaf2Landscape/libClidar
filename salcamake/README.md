# salcaMake

This repository contains C code to process data from the SALCA terrestrial laser scanner (TLS). The SALCA instrument is described in [Danson et al (2014)](https://www.sciencedirect.com/science/article/pii/S0168192314001749) and this code library was first used in [Hancock et al (2017)](https://ieeexplore.ieee.org/abstract/document/7879210).

SALCA is a full-waveform, dual-wavelength TLS and the raw data comes in a custom binary format (.bin).


The main tools are

**hedgehog**: A C program to process SALCA waveform data to extract either a point cloud or a polar projection file.

**hedgehog.csh**: A wrapper around **hedgehog** to ease batch processing.

The other programs are either for calibrating SALCA returns to reflectance or for calculating the angle offsets of SALCA's squinty lasers.


## Installing

### Dependencies

All programs depend on a library a without package manager, which must be cloned:

* [C-tools](https://bitbucket.org/StevenHancock/tools)

And point to its locations with the environment variables:

    HANCOCKTOOLS_ROOT


### Compiling

Compile by typing:

  **make THIS=hedgehog**
  **make THIS=hedgehog install**

Then make sure that **hedgehog** and **hedgehog.csh** are in your PATH.





## Usage 

Command line options are

    -inRoot root;          input filename root
    -outRoot root;         output filename root
    -writeBin;             write output as binary instead of ASCII
    -nAz n;                number of azimuth steps (number of binary files)
    -azStep ang;           azimuth step in degrees
    -maxR range;           maximum recorded range, metres
    -calibrate;            calibrate to reflectance, need calibration file
    -calFile name;         calibrate to reflectance, need calibration file
    -laser p0 p1;          outgoing laser power for each wavelength
    -filt n;               filter optical depth (0, 0.6, 1 or 1.6)
    -nFilt t1 t2;          specify transmissions rather than use defaults
    -coarsen n;            coarsen by a factor. CAUTIION, intensities will not scale properly yet
    -maxZen zen;           maximum zenith angle, degrees. 190 by default
    -azSquint ang0 ang1;   azimuth squint angle in degrees for each band
    -zenSquint ang0 ang1;  zenith squint angle, degrees for each band
    -mSquint ang           mirror squint anbgle
    -azStart angle;        azimuth start, degrees
    -header file;          read header
    -back DN;              background DN value
 
    -sWidth width;         smoothing width for quadratic fitting
    -quad;                 quadratic fitting
 
    ### Geolocation
    -trans name;           collocation matrix filename
    -globName name;        geolocation matrix filename


A usage example is:

    hedgehog.csh -inRoot /exports/csce/datastore/geos/groups/3d_env/data/TLS/SALCA/cockle_park/lai.2.06_2832/lai.2.06_2832 -header /exports/csce/datastore/geos/groups/3d_env/data/TLS/SALCA/headers/lai.2.06_2832.hdr -trans /exports/csce/datastore/geos/groups/3d_env/data/msc_results/2019/feng_zemin/matrix/2832.1.mat.txt

Here *lai.2.06\_2832.hdr* is a header file describing the SALCA binary files (resolution, number of shots, number of waveform bins etc.), and 2832.1.mat.txt is a geolocation matrix to but the data in to the appropriate coordinate system.


### Output

By default the code outputs in .pts format. If the *-writeBin* flag is selected at run-time, the data will be output in a more compact custom binary format, including the beams which did not return a hit. This can be used for calculating gap fraction with this [tool](https://bitbucket.org/StevenHancock/voxelate).


# Any issues

Contact svenhancock@gmail.com

