#!/bin/csh -f

set bin="/home/server/users/bakgrp3/nsh103/src/salcaMake"
set input="teast_0_2161.band.0.coords"

set minX=-1000
set maxX=1000
set minY=-1000
set maxY=1000
set minZ=-1000
set maxZ=1000
set maxSat=2     # don't filter for saturation by default
set minRefl=-300 # output all returns

set minR=0
set maxR=1000
@ minAzI=-1
@ maxAzI=300000
@ minZenI=-1
@ maxZenI=300000

@ joy=1


while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -satTest
    set maxSat=0.5
  shift argv
  breaksw

  case -minRefl
    set minRefl="$argv[2]"
  shift argv;shift argv
  breaksw

  case -r
    set minR=$argv[2]
    set maxR=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -azInd
    @ minAzI=$argv[2]
    @ maxAzI=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -zenInd
    @ minZenI=$argv[2]
    @ maxZenI=$argv[3]
  shift argv;shift argv;shift argv
  breaksw

  case -joyless
    @ joy=0
  shift argv
  breaksw

  case -help
    echo " "
    echo "-input name;     input filename"
    echo "-r min max;      range bounds"
    echo "-azInd min max;  az index bounds"
    echo "-zenInd min max; zen index bounds"
    echo "-satTest;        filter out saturated"
    echo "-minRefl r;      minimum signal to output"
    echo "-joy;            no progress messages"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

set temp="seperateWorkspace.$$.dat"


gawk -f $bin/extractBounds.awk -v minX=$minX -v maxX=$maxX -v minY=$minY -v maxY=$maxY -v minZ=$minZ -v maxZ=$maxZ -v minRefl=$minRefl -v maxSat=$maxSat -v minZenI=$minZenI -v maxZenI=$maxZenI -v minAzI=$minAzI -v maxAzI=$maxAzI -v minR=-10 -v maxR=1000 < $input > $temp

gawk -f $bin/avoidMult.awk -v minZenI=$minZenI -v maxZenI=$maxZenI -v minAzI=$minAzI -v maxAzI=$maxAzI -v minR=$minR -v maxR=$maxR -v root="$input:r" < $temp 

if( $joy )echo Written to $input:r.hard and $input:r.soft""

if( -e $temp )rm $temp

