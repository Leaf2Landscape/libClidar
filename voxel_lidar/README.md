# README #

This is a set of programs to produce maps of plant area index and voxelised gap fraction from airborne laser scanning (ALS) data. The programs in this repository calibrate full-wavewform ALS processing against the higher resolution terrestrial laser scanne (TLS) data. Programs to create voxelised plant area index (PAI) maps from ALS or TLS data alone are available in the [voxelate](https://bitbucket.org/StevenHancock/voxelate) repository. The method is described in:

[Hancock, S., Anderson, K., Disney, M. and Gaston, K.J., 2017. Measurement of fine-spatial-resolution 3D vegetation structure with airborne waveform lidar: Calibration and validation with voxelised terrestrial lidar. Remote Sensing of Environment, 188, pp.37-50.](http://www.sciencedirect.com/science/article/pii/S0034425716304205)


This depends upon the [libClidar](https://bitbucket.org/StevenHancock/libclidar) library, which also underpins the [GEDI simulator](https://bitbucket.org/StevenHancock/gedisimulator) set of tools. The latter is described in [Hancock et al (2019)](https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2018EA000506).



## Usage ##

The main tool is **optimiseTLS**. This calibrates a full-waveform airborne lidar (ALS) processing method using terrestrial lidar (TLS) for a subset of ALS beams. It required full-waveform ALS data input in .las format.

**rxpVox** assesses the accuracy over the full ALS and TLS datasets.


The TLS processing requires TLS files to be converted into a custom binary using the programs in [rxptobinary](https://bitbucket.org/StevenHancock/rxptobinary).

Currently this only works for Riegl's ".rxp" format. Other formats could be added, but the beams that do not record interactions must be recorded (or be possible to infer) in order to calculate gap fractions.

The command line options are given below, and can be seen with the -help command.

    # Input/output options
    -input name;     input filename
    -inList name;    list of input files to solve simultaneously
    -output name;    output filename
    -pFile name;     pulse filename
    -statsOut name;  stats for multi-fit output
    -optDenoise;     optimise ALS denoising
    -denWidth s;     smoothing width when testing denoising accuracy
    -atten;          optimise with attenuation correction
    -GBIC;           apply GBIC
    -readTLSpar;     read TLS par from file rather than optimise
    -denoise;        perform denoising without optimisation

    # Denoising options
    -varNoise;       use variable noise levels
    -richLucy;       Richardson-Lucy deconvolution rather than Gold's

    # Initial parameter guess
    -thresh x;
    -meanN x;
    -tailThresh x;
    -sWidth x;
    -psWidth x;
    -minWidth x;
    -medLen x;
    -pScale x;
    -deChang x;
    -varScale x;     variable noise threshold scale
    -gauss;          fit a Gaussian
    -gWidth sig;     smooth to choose Gaussian features
    -fixPS;          fix pulse scaling
    -noWeight;       do not weight fit to avoid subterranean
    -fixSmoo;        fix smoothing width
    -optPsmoo;       optimise pre-smoothing
    -firstLast;      optimise to find first and last return
    -noTrack;        do not use noise tracking
    -flWeight w;     scale factor for first-last weighting
    -lockGap;       set min gap fraction to 1
    -findHard;      find hard targets




## Compiling

All programs depend on these libraries:

* Gnu Scientific Library
* Geotiff
* [Minpack (Levenberg-Maruqardt)](https://www.physics.wisc.edu/~craigm/idl/cmpfit.html)

Once they're installed, it also requires three libraries without package managers:

* [C-tools](https://bitbucket.org/StevenHancock/tools)
* [libClidar](https://bitbucket.org/StevenHancock/libclidar)
* [voxelate](https://bitbucket.org/StevenHancock/voxelate)


Clone these from bitbucket and point to their locations with the environment variables:

    HANCOCKTOOLS_ROOT
    LIBCLIDAR_ROOT
    CMPFIT_ROOT
    VOXELATE_ROOT

To compile, type:

  **make THIS=optimiseTLS**

Make sure that ~/bin/$ARCH exists and is in the path, where $ARCH is the result of `uname -m`.

  **make THIS=optimiseTLS install**

Replace "**optimiseTLS**" with each of the commands above to compile and install.






### Who do I talk to? ###

svenhancock@gmail.com
