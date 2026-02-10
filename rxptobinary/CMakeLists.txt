# $Id: CMakeLists.txt 1059 2014-07-09 08:26:09Z RS $
# This is the CMake project file for configuration of the RiVLib examples.
#
#           Copyright (c) 2010 Riegl Laser Measurement Systems.

# The authors hereby  grant permission  to use,  and copy this software  for the
# purpose  of enhancing the useability  of  Riegl Laser Measurement Systems GmbH
# instruments you own. You may NOT distribute or modify the software for the use
# in commercial applications without the written consent of RLMS.
#
# In any case copyright notices and  this notice must be  retained in all copies
# verbatim in any distributions.

# IN  NO EVENT SHALL  THE AUTHORS OR  DISTRIBUTORS BE  LIABLE TO  ANY PARTY  FOR
# DIRECT, INDIRECT, SPECIAL,  INCIDENTAL, OR  CONSEQUENTIAL DAMAGES  ARISING OUT
# OF THE USE OF THIS SOFTWARE, ITS  DOCUMENTATION,  OR ANY DERIVATIVES  THEREOF,
# EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM  ANY WARRANTIES, INCLUDING,
# BUT NOT LIMITED TO, THE IMPLIED WARRANTIES  OF MERCHANTABILITY,  FITNESS FOR A
# PARTICULAR PURPOSE,  AND  NON-INFRINGEMENT. THIS SOFTWARE  IS  PROVIDED  ON AN
# "AS IS" BASIS, AND THE AUTHORS  AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE
# MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

# This file is the input description for CMake, a tool which is able to
# generate build environments for various compilers.
# You can obtain a copy of CMake at http://www.cmake.org
# It is recommended to use the graphical tool from CMake to do the
# configuration. First you have to specify the source directory i.e. where
# the RiVlib examples are. Next you specify a directory where you want your
# build to take place. Typically this is the directory where the project files
# of your IDE will be created.
# Example:
# Source dir: ...../examples
# Build dir:  ...../examples-build
# The next step is to press the configure button. CMake will prompt you
# with a message and a variable in red, requesting you to specify the
# root directory of your variant of rivlib you will use. Use the ... navigate
# button to specify the subdirectory within your downloaded copy of RiVLib that
# fits your compiler.
# Example:
# RiVLib_DIR = .../rivlib-1_30-x86-windows-vc120
# Repeat pressing the configuration button and resolve any more red lines,
# unitil there are no unsatisfied constraints remaining.
# Press generate, and the project files will be created.
# Now you can launch your favorite IDE or invoke the genrated Makefile.
#
# NOTE:
# RiVLib needs a compiler with C++11 support. Most recent compilers
# already do have such support. If using GCC you might need to add the
# flag -std=gnu++0x to the CMAKE_CXX_FLAGS variable.
#

if(UNIX)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11 ")
endif()

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lhdf5")
SET(CXXFLAGS "${CXXFLAGS} -std=c++11 -lhdf5")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(HDF5 hdf5)


cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

include( InstallRequiredSystemLibraries )

project(readRXP)

find_package(RiVLib
    COMPONENTS scanlib scanifc ctrllib ctrlifc
)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/bin
        CACHE PATH
        "Install path prefix, prepended onto install directories."
        FORCE
    )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

include_directories(
    ${RiVLib_INCLUDE_DIRS} $ENV{LIBCLIDAR_ROOT} $ENV{HANCOCKTOOLS_ROOT}
)


add_executable( readRXP
    readRXP.cpp
)

target_include_directories(readRXP PUBLIC
                           $ENV{HANCOCKTOOLS_ROOT}
                           $ENV{LIBCLIDAR_ROOT}
)

target_link_libraries( readRXP
    ${RiVLib_SCANLIB_LIBRARY}
    ${HDF5}
)
install(
    TARGETS readRXP
    RUNTIME DESTINATION .
)

install( 
    PROGRAMS ${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS} 
    DESTINATION . 
    COMPONENT System 
)
