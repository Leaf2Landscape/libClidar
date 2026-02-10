This converts Riegl's .rxp format to the binary format needed by:

https://bitbucket.org/StevenHancock/voxelate

https://bitbucket.org/StevenHancock/voxel_lidar



Set an environment variable:

setenv RiVLib_DIR "/home/shancock/src/downloaded/rivlib-2_5_11-x86_64-linux-gcc5"
setenv CXXFLAGS "-std=c++11"

Note that gcc version 5 or later needs the above file. Version 4 needs rivlib-2_5_7-x86_64-linux-gcc44


# Run 

cmake .
make

