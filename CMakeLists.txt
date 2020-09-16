#  Copyright &copy; 2020, Dr. Stefan Sicklinger, Munich
#
#  All rights reserved.
#
#  EVAA: Efficient Vehicle dynAmics simulAtor
#
#  This file is part of EVAA.
#
#  EVAA is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  EVAA is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with EVAA. If not, see http://www.gnu.org/licenses/.
#------------------------------------------------------------------------------------#
cmake_minimum_required (VERSION 3.14)
#------------------------------------------------------------------------------------#
project (EVAA VERSION 1.0 LANGUAGES CXX)
set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo")
#------------------------------------------------------------------------------------#)
#------------------------------------------------------------------------------------#)
# Check OS and host system
IF((${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")) 
   SET(OS_LINUX_X86_64 TRUE)
   MESSAGE("EVAA is going to be configured for Linux x86_64!")
ELSEIF((${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "AMD64")) 
   SET(OS_WIN_X86_64 TRUE)
   MESSAGE("EVAA is going to be configured for Windows x64!")
ELSEIF((${CMAKE_SYSTEM_NAME} MATCHES "Darwin" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")) 
   SET(OS_MAC_X86_64 TRUE)
   MESSAGE("EVAA is going to be configured for Mac x86_64!")
ELSE() 
   MESSAGE(FATAL_ERROR "EVAA is supported on Linux x86_64 and Windows x64 and Mac x86_64 only")   
ENDIF()
IF(CMAKE_OUTPUT_VERBOSE)
   MESSAGE("Current source dir is " ${EVAA_SOURCE_DIR})
ENDIF()
#------------------------------------------------------------------------------------#
# 1. Options
# 2. Macros
# 3. Compiler
# 5. Configure
# 4. Subdirectories
#------------------------------------------------------------------------------------#
######################################################################################
#1. Option
######################################################################################
option(BUILD_UNITTESTS           "Build unit test and micro benchmarks" ON)
option(CMAKE_OUTPUT_VERBOSE      "Add Cmake debug output" ON)
######################################################################################
#2. Macros
######################################################################################
# Add EVAA cmake modules
SET(CMAKE_MODULE_PATH ${EVAA_SOURCE_DIR}/cmake/modules)
include(MacroOutOfSourceBuild)
include(MacroAppendGlobalVariable)
include(GetGitRevisionDescription)
# Download automatically, you can also just copy the conan.cmake file
if(NOT EXISTS "${CMAKE_MODULE_PATH}/conan.cmake")
   message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
   file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
                  "${CMAKE_MODULE_PATH}/conan.cmake")
endif()
include(conan)
######################################################################################
#3. Compiler
######################################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS NO)
######################################################################################
#4. Configure
######################################################################################
#------------------------------------------------------------------------------------#
# Avoid in source builds
macro_ensure_out_of_source_build("In source build is not allowed for ${PROJECT_NAME}.\
Please run CMake from your build directory. You may want to remove CMakeCache.txt.")
#------------------------------------------------------------------------------------#
# Build unit test and micro benchmarks
IF (${BUILD_UNITTESTS}) 
   MESSAGE("Using Catch2 for unit tests and micro benchmarks")
   conan_cmake_run(REQUIRES catch2/2.13.0 BASIC_SETUP CMAKE_TARGETS)
   include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
ENDIF()
#------------------------------------------------------------------------------------#
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_TAG --tags)
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/tools/AuxiliaryParameters.cpp.in" 
"${CMAKE_CURRENT_SOURCE_DIR}/src/tools/AuxiliaryParameters.cpp" @ONLY)
######################################################################################
#5. Subdirectories THE ORDER IS IMPORTANT, BECAUSE THEY ARE COMPILED IN THAT ORDER
######################################################################################
#------------------------------------------------------------------------------------#
# Where to put EXECUTABLE
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin64 )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin64 )
#------------------------------------------------------------------------------------#
MESSAGE("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
# INCLUDE dirs
add_subdirectory(src)
IF(${BUILD_UNITTESTS})
enable_testing()
add_subdirectory(unitTests)
ENDIF()
