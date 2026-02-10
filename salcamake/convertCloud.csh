#!/bin/csh -f


#########################
# converts a salcaMake  #
# cloud to matlab style #
#########################


set bin="/home/server/users/bakgrp3/nsh103/src/salcaMake"
set input="del_0407_10_06_1441.band.0.coords"
set output="del_0407_10_06_1441.band.0.matlab"

# read the command line
while ($#argv>0)
  switch("$argv[1]")

  case -input
    set input="$argv[2]"
  shift argv;shift argv
  breaksw

  case -output
    set output="$argv[2]"
  shift argv;shift argv
  breaksw

  case -help
    echo " "
    echo "-output name;   output filename"
    echo "-input name;    input filename"
    echo " "
    exit

  default:
    echo "Unknown argument $argv[1]"
    echo "Who knows"
    echo "Type echi.rat -help"
    exit;
  
  endsw
end


gawk -f $bin/convertCloud.awk < $input > $output

echo "Written to $output"

