#!/bin/csh -f


#######################
# Creates a histogram #
# from a subset       #
#######################


set bin="/home/server/users/bakgrp3/nsh103/src/salcaMake"
set input="teast_0_2161.band.0.coords"

set minX=-1000
set maxX=1000
set minY=-1000
set maxY=1000
set minZ=-1000
set maxZ=1000
set minR=0
set maxR=200
set minAz=-1
set maxAz=300000
set minZen=-1
set maxZen=320000
@ nBins=20
@ col=4

while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -x
    set minX=$argv[2]
    set maxX=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -y
    set minY=$argv[2]
    set maxY=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -z
    set minZ=$argv[2]
    set maxZ=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -r
    set minR=$argv[2]
    set maxR=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -azInd
    set minAz=$argv[2]
    set maxAz=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -zenInd
    set minZen=$argv[2]
    set maxZen=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -nBins
    @ nBins=$argv[2]
  shift argv;shift argv
  breaksw

  case -integral
    @ col=12
  shift argv
  breaksw

  case -col
    @ col=$argv[2]
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-input name;      input filename"
    echo "-x min max;       x bounds"
    echo "-y min max;       y bounds"
    echo "-z min max;       z bounds"
    echo "-r min max;       range bounds"
    echo "-azInd min max;   azimuth bounds"
    echo "-zenInd min max;  zenith bounds"
    echo "-nBins n;         number of histogram bins"
    echo "-integral;        use integral rather than max"
    echo "-col n;           column to make histogram"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

set output="$input:r.hist"

gawk -f $bin/subHist.awk -v col=$col -v minX=$minX -v maxX=$maxX -v minY=$minY -v maxY=$maxY -v minZ=$minZ -v maxZ=$maxZ -v nBins=$nBins -v minR=$minR -v maxR=$maxR -v minAz=$minAz -v maxAz=$maxAz -v minZen=$minZen -v maxZen=$maxZen < $input > $output

echo "Written to $output"

