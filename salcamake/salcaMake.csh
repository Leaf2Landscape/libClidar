#!/bin/csh -f


###############################
# A wrapper for salcaMake.c   #
# S. Hancock 15th August 2013 #
###############################



# defaults
set nAz=3050
set azStep=0.06
set maxR=30
set subBin=" "
set intTrack=" "
set calibrate=" "
set calFile=" "
set filt=" "
set nFilt=" "
set dualOut=" "
set coarsen=" "
set zenOffset=" "
set maxZen=" "
set azSquint=" "
set zenSquint=" "
set azStart=" "
set laserP=" "
set outRoot="salcaTest"
set back=" " # background level
set noQuad=" "
set noLM=" "

@ readHead=0
set bin="$HOME/src/salcaMake"


set comm="salcaMake"
set sWidth=" "
set lmWidth=" "


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

  case -subBin
    set subBin="-subBin"
  shift argv
  breaksw

  case -intTrack
    set intTrack="-intTrack"
  shift argv
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

  case -zenOffset
    set zenOffset="-zenOffset"
  shift argv
  breaksw

  case -maxZen
    set maxZen="-maxZen $argv[2]"
  shift argv;shift argv
  breaksw

  case -azSquint
    set azSquint="-azSquint $argv[2]"
  shift argv;shift argv
  breaksw

  case -zenSquint
    set zenSquint="-zenSquint $argv[2]"
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

  case -startAgain
    set comm="startAgain"
  shift argv
  breaksw

  case -hedgehog
    set comm="hedgehog"
  shift argv
  breaksw

  case -sWidth
    set sWidth="-sWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -lmWidth
    set lmWidth="-lmWidth $argv[2]"
  shift argv;shift argv
  breaksw

  case -noLM
    set noLM="-noLM"
  shift argv
  breaksw

  case -noQuad
    set noQuad="-noQuad"
  shift argv
  breaksw

  case -laser
    set laserP="-laser $argv[2] $argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-inRoot root;      input filename root"
    echo "-outRoot root;     output filename root"
    echo "-nAz n;            number of azimuth steps (number of binary files)"
    echo "-azStep ang;       azimuth step in degrees"
    echo "-maxR range;       maximum recorded range, metres"
    echo "-subBin;           sub-bin range resolution;"
    echo "-intTrack;         track to background noise"
    echo "-calibrate;        calibrate to reflectance, need calibration file"
    echo "-calFile name;     calibrate to reflectance, need calibration file"
    echo "-laser p0 p1;      laser power correction factors. Replaces filters below."
    echo "-filt n;           filter optical depth (0, 0.6, 1 or 1.6)"
    echo "-nFilt t1 t2;      specify transmissions rather than use defaults"
    echo "-dualOut;          output a combined point cloud"
    echo "-coarsen n;        coarsen by a factor. CAUTIION, intensities will not scale properly yet"
    echo "-zenOffset;        offset 1064nm zenith by 1"
    echo "-maxZen zen;       maximum zenith angle, degrees. 190 by default"
    echo "-azSquint angle;   azimuth squint angle in degrees"
    echo "-zenSquint angle;  zenith squint angle, degrees"
    echo "-azStart angle;    azimuth start, degrees"
    echo "-header file;      read header"
    echo "-back DN;          background DN value"
    echo " "
    echo "-hedgehog;         final version of code"
    echo "-startAgain;       do quadratic and Gaussian fitting"
    echo "-sWidth width;     smoothing width for quadratic fitting"
    echo "-lmWidth width;    smoothing width for LM fitting"
    echo "-noLM;             no LM fitting"
    echo "-noQuad;             no quadratic fitting"
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

nice +19 $comm -inRoot $inRoot -outRoot $outRoot -nAz $nAz -azStep $azStep -maxR $maxR $subBin $intTrack $calibrate $calFile $filt $nFilt $dualOut $coarsen $zenOffset $maxZen $azSquint $zenSquint $azStart $back $sWidth $lmWidth $noQuad $noLM $laserP

echo "Ping"

