#!/bin/csh -f
  
set list="CMakeCache.txt Tutorial CMakeFiles Makefile cmake_install.cmake"

foreach file( $list )
  if( -e $file )then
    echo "Cleaning $file"
    rm -r $file
  endif
end

