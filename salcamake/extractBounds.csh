#!/bin/csh -f

set bin="$HOME/src/salcaMake"
set input="teast_0_2161.band.0.coords"

set minX=-1000
set maxX=1000
set minY=-1000
set maxY=1000
set minZ=-1000
set maxZ=1000
set maxSat=200000     # don't filter for saturation by default
set minRefl=-300 # output all returns

set minR=-100000
set maxR=100000
@ minAzI=-1000000
@ maxAzI=300000
@ minZenI=-100000
@ maxZenI=300000

@ joy=1


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
    echo "-x min max;      x bounds"
    echo "-y min max;      y bounds"
    echo "-z min max;      z bounds"
    echo "-r min max;      range bounds"
    echo "-azInd min max;  az index bounds"
    echo "-zenInd min max; zen index bounds"
    echo "-satTest;        filter out saturated"
    echo "-minRefl r;      minimum signal to output"
    echo "-joyless;        don't print out status"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

set output="$input:r.filt"

gawk -f $bin/extractBounds.awk -v minX=$minX -v maxX=$maxX -v minY=$minY -v maxY=$maxY -v minZ=$minZ -v maxZ=$maxZ -v minRefl=$minRefl -v maxSat=$maxSat -v minR=$minR -v maxR=$maxR -v minZenI=$minZenI -v maxZenI=$maxZenI -v minAzI=$minAzI -v maxAzI=$maxAzI < $input > $output

if( $joy )echo "Written to $output"

