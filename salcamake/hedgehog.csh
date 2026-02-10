#!/bin/csh -f


###############################
# A wrapper for salcaMake.c   #
# S. Hancock 15th August 2013 #
###############################



# defaults
set nAz=3050
set azStep=0.06
set maxR=30
set calibrate=" "
set calFile=" "
set filt=" "
set nFilt=" "
set coarsen=" "
set maxZen=" "
set azSquint=" "
set zenSquint=" "
set mSquint=" "
set azStart=" "
set outRoot="salcaTest"
set back=" "
set quad=" "
set writeBin=" "
set inRoot="teast"
set trans=" "
set globName=" "
set laser=" "


set sWidth=" "
set lmWidth=" "

@ readHead=0
set bin="$HOME/src/salcaMake"




# read command line
while ($#argv>0)
  switch("$argv[1]")

  case -inRoot
    set inRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -outRoot
    set outRoot="$argv[2]"
  shift argv;shift argv
  breaksw

  case -nAz
    set nAz=$argv[2]
  shift argv;shift argv
  breaksw

  case -azStep
    set azStep="$argv[2]"
  shift argv;shift argv
  breaksw

  case -maxR
    set maxR="$argv[2]"
  shift argv;shift argv
  breaksw

  case -calibrate
    set calibrate="-calibrate"
  shift argv
  breaksw

  case -calFile
    set calFile="-calFile $argv[2]"
  shift argv;shift argv
  breaksw

  case -filt
    set filt="-filt $argv[2]"
  shift argv;shift argv
  breaksw

  case -nFilt
    set nFilt="-nFilt $argv[2] $argv[3]"
  shift argv;shift argv;shift argv
  breaksw
    
  case -dualOut
    set dualOut="-dualOut"
  shift argv
  breaksw

  case -coarsen
    set coarsen="-coarsen $argv[2]"
  shift argv;shift argv
  breaksw

  case -maxZen
    set maxZen="-maxZen $argv[2]"
  shift argv;shift argv
  breaksw

  case -azSquint
    set azSquint="-azSquint $argv[2] $argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -zenSquint
    set zenSquint="-zenSquint $argv[2] $argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -mSquint
    set mSquint="-mSquint $argv[2]"
  shift argv;shift argv
  breaksw

  case -azStart
    set azStart="-azStart $argv[2]"
  shift argv;shift argv
  breaksw

  case -header
    @ readHead=1
    set header="$argv[2]"
  shift argv;shift argv
  breaksw

  case -background
    set back="-background $argv[2]"
  shift argv;shift argv
  breaksw

  case -sWidth
    set sWidth="-sWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -quad
    set quad="-quad"
  shift argv
  breaksw

  case -writeBin
    set writeBin="-writeBin"
  shift argv
  breaksw

  case -trans
    set trans="-trans $argv[2]"
  shift argv;shift argv
  breaksw

  case -laser
    set laser="-laser $argv[2] $argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -globName
    set globName="-trans $argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inRoot root;          input filename root"
    echo "-outRoot root;         output filename root"
    echo "-writeBin;             write output as binary instead of ASCII"
    echo "-nAz n;                number of azimuth steps (number of binary files)"
    echo "-azStep ang;           azimuth step in degrees"
    echo "-maxR range;           maximum recorded range, metres"
    echo "-calibrate;            calibrate to reflectance, need calibration file"
    echo "-calFile name;         calibrate to reflectance, need calibration file"
    echo "-laser p0 p1;          outgoing laser power for each wavelength"
    echo "-filt n;               filter optical depth (0, 0.6, 1 or 1.6)"
    echo "-nFilt t1 t2;          specify transmissions rather than use defaults"
    echo "-coarsen n;            coarsen by a factor. CAUTIION, intensities will not scale properly yet"
    echo "-maxZen zen;           maximum zenith angle, degrees. 190 by default"
    echo "-azSquint ang0 ang1;   azimuth squint angle in degrees for each band"
    echo "-zenSquint ang0 ang1;  zenith squint angle, degrees for each band"
    echo "-mSquint ang           mirror squint anbgle"
    echo "-azStart angle;        azimuth start, degrees"
    echo "-header file;          read header"
    echo "-back DN;              background DN value"
    echo " "
    echo "-sWidth width;         smoothing width for quadratic fitting"
    echo "-quad;                 quadratic fitting"
    echo " "
    echo "Geolocation"
    echo "-trans name;           collocation matrix filename"
    echo "-globName name;        geolocation matrix filename"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

# if specified, read header
if( $readHead )then
  set bits=`gawk -f $bin/readHead.awk` < $header
  set maxR=`echo $bits[1]|sed -e s%m%" "%`
  set azStep=$bits[2]
  set nAz=$bits[3]
  set azStart="-azStart $bits[4]"
  set maxZen="-maxZen $bits[5]"
endif

hedgehog -inRoot $inRoot -outRoot $outRoot -nAz $nAz -azStep $azStep -maxR $maxR $calibrate $calFile $filt $nFilt $coarsen $maxZen $azSquint $zenSquint $mSquint $azStart $back $sWidth $quad $writeBin $trans $globName $laser


echo "Ping"

