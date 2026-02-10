#!/bin/csh -f

set bin="$HOME/src/voxelate"
set minX=-10000000
set maxX=100000000
set minY=-10000000
set maxY=100000000
set minZ=-10000000
set maxZ=100000000

@ setName=1

# read command line
while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
    @ setName=0
  shift argv;shift argv
  breaksw

  case -x
    set minX="$argv[2]"
    set maxX="$argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -y
    set minY="$argv[2]"
    set maxY="$argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -z
    set minZ="$argv[2]"
    set maxZ="$argv[3]"
  shift argv;shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-input name;  input filename"
    echo "-output name; output filename, if needed"
    echo "-x min max;   x bounds"
    echo "-y min max;   y bounds"
    echo "-z min max;   z bounds"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;

  endsw
end

if( $setName )then
  set output="$input:r.strip"
endif

gawk -f $bin/trimCloud.awk -v minX=$minX -v maxX=$maxX -v minY=$minY -v maxY=$maxY -v minZ=$minZ -v maxZ=$maxZ < $input > $output
echo "Written to $output"

