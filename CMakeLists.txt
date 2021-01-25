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
#  along with EVAA. If not, see http://www.gnu.org/licenses/
#------------------------------------------------------------------------------------#
CMAKE_MINIMUM_REQUIRED(VERSION 3.18)
#------------------------------------------------------------------------------------#
PROJECT(EVAA VERSION 1.0 LANGUAGES CXX)
SET(CMAKE_CONFIGURATION_TYPES "Release;Debug")
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
# 1. Options
# 2. Macros and Checks
# 3. Compiler
# 5. Configure
# 4. Subdirectories
#------------------------------------------------------------------------------------#
######################################################################################
#1. Option
######################################################################################
OPTION(BUILD_UNITTESTS           "Build unit test and micro benchmarks" ON)
OPTION(USE_HDF5                  "Use HDF5"                             ON)
OPTION(CMAKE_OUTPUT_VERBOSE      "Add Cmake debug output"               ON)
######################################################################################
#2. Macros and Checks
######################################################################################
#------------------------------------------------------------------------------------#
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
   MESSAGE("Current source dir is " ${PROJECT_SOURCE_DIR})
ENDIF()
#------------------------------------------------------------------------------------#
# Add EVAA cmake modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/helperScripts/cmake/modules)
INCLUDE(MacroOutOfSourceBuild)
INCLUDE(MacroAppendGlobalVariable)
INCLUDE(GetGitRevisionDescription)
# Download automatically, you can also just copy the conan.cmake file
IF(NOT EXISTS "${CMAKE_MODULE_PATH}/conan.cmake")
   message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
   file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
                  "${CMAKE_MODULE_PATH}/conan.cmake")
ENDIF()
INCLUDE(conan)
#------------------------------------------------------------------------------------#
######################################################################################
#3. Compiler
######################################################################################
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
######################################################################################
#4. Configure
######################################################################################
#------------------------------------------------------------------------------------#
# Avoid in source builds
macro_ensure_out_of_source_build("In source build is not allowed for ${PROJECT_NAME}.\
Please run CMake from your build directory. You may want to remove CMakeCache.txt.")
#------------------------------------------------------------------------------------#
# Where to put EXECUTABLE
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
#------------------------------------------------------------------------------------#
# Build list of conan packages
IF(${BUILD_UNITTESTS})
   LIST(APPEND CONAN_PACKAGES "gtest/1.10.0")
   LIST(APPEND CONAN_PACKAGES "benchmark/1.5.2")
ENDIF()
IF(${USE_HDF5})
   LIST(APPEND CONAN_PACKAGES "hdf5/1.12.0")
ENDIF()
LIST(APPEND CONAN_PACKAGES "xerces-c/3.2.2")
LIST(APPEND CONAN_PACKAGES "boost/1.74.0")
#------------------------------------------------------------------------------------#
IF(OS_MAC_X86_64)
   # see https://github.com/conan-io/conan-center-index/issues/4097
   conan_cmake_run(REQUIRES ${CONAN_PACKAGES} BASIC_SETUP CMAKE_TARGETS BUILD missing SETTINGS compiler.cppstd=14)
ELSE()
   conan_cmake_run(REQUIRES ${CONAN_PACKAGES} BASIC_SETUP CMAKE_TARGETS BUILD missing)
ENDIF()
#------------------------------------------------------------------------------------#
# Configure options
IF(${BUILD_UNITTESTS}) 
   MESSAGE("Using GTest and GBench for unit tests and micro benchmarks")
   INCLUDE(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
ENDIF()
IF(${USE_HDF5})
   MESSAGE("Using HDF5")
   ADD_DEFINITIONS(-DUSE_HDF5)
   IF(OS_WIN_X86_64)
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} /ignore:4099")
   ENDIF()
ENDIF()
#------------------------------------------------------------------------------------#
# Install CODESYNTHESIS XSD
IF(OS_WIN_X86_64)
   SET(CODESYNTHESIS_PATH "https://codesynthesis.com/download/xsd/4.0/windows/i686/xsd-4.0.0-i686-windows.zip")
   SET(CODESYNTHESIS_FILE "xsd-4.0.0-i686-windows")
   SET(CODESYNTHESIS_FILE_END ".zip")
   SET(CODESYNTHESIS_FILE_HASH "d129469e109784c663387ca8bee5ac627434cfca")
ELSEIF(OS_LINUX_X86_64)
   SET(CODESYNTHESIS_PATH "https://codesynthesis.com/download/xsd/4.0/linux-gnu/x86_64/xsd-4.0.0-x86_64-linux-gnu.tar.bz2")
   SET(CODESYNTHESIS_FILE "xsd-4.0.0-x86_64-linux-gnu")
   SET(CODESYNTHESIS_FILE_END ".tar.bz2")
   SET(CODESYNTHESIS_FILE_HASH "5eeb2eeca0d893949e3677bb374e7b96f19770d6")
ELSEIF(OS_MAC_X86_64)
   SET(CODESYNTHESIS_PATH "https://sicklinger.com/xsd/xsd-4.0.0-i686-macosx.tar.bz2")
   SET(CODESYNTHESIS_FILE "xsd-4.0.0-i686-macosx")
   SET(CODESYNTHESIS_FILE_END ".tar.bz2")
   SET(CODESYNTHESIS_FILE_HASH "f3f850e4c335f375ddda07c75902c439b0d6b0cc")  
   ADD_DEFINITIONS(-DXSD_CXX11)
   ADD_DEFINITIONS(-DXSD_CXX11_TEMPLATE_ALIAS)
   ADD_DEFINITIONS(-DXSD_CXX11_NULLPTR)
ELSE() 
   MESSAGE(FATAL_ERROR "EVAA is supported on Linux x86_64 and Windows x64 and Mac x86_64 only")   
ENDIF()

IF(NOT EXISTS ${PROJECT_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END})
   FILE(DOWNLOAD ${CODESYNTHESIS_PATH} ${PROJECT_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END}
      TIMEOUT 60
      EXPECTED_HASH SHA1=${CODESYNTHESIS_FILE_HASH}
      TLS_VERIFY ON)
ENDIF()
FILE(ARCHIVE_EXTRACT INPUT ${PROJECT_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END} DESTINATION ${PROJECT_SOURCE_DIR}/xsd)
FIND_PROGRAM(XSD_EXE xsd PATHS ${PROJECT_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}/bin REQUIRED NO_DEFAULT_PATH)
IF(CMAKE_OUTPUT_VERBOSE)
   MESSAGE("XSD_EXE ${XSD_EXE}")
ENDIF()
# Run CODESYNTHESIS XSD
#https://www.codesynthesis.com/projects/xsd/documentation/xsd.xhtml
EXECUTE_PROCESS(  
   COMMAND ${XSD_EXE} cxx-tree --std c++11 --hxx-suffix .h --cxx-suffix .cpp "${PROJECT_SOURCE_DIR}/xsd/InputSchemaEVAA.xsd"
   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/src/database)
IF(NOT EXISTS ${PROJECT_SOURCE_DIR}/src/database/InputSchemaEVAA.h )
   MESSAGE(FATAL_ERROR  "XSD has not generated ${PROJECT_SOURCE_DIR}/src/database/InputSchemaEVAA.h")
ENDIF()
SET(CODESYNTHESIS_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}/libxsd")
#------------------------------------------------------------------------------------#
# Configure files
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
git_describe(GIT_TAG --tags)
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/src/tools/AuxiliaryParameters.cpp.in" 
"${CMAKE_CURRENT_SOURCE_DIR}/src/tools/AuxiliaryParameters.cpp" @ONLY)
######################################################################################
#5. Subdirectories THE ORDER IS IMPORTANT, BECAUSE THEY ARE COMPILED IN THAT ORDER
######################################################################################
#------------------------------------------------------------------------------------#
MESSAGE("CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
# INCLUDE dirs
ADD_SUBDIRECTORY(src)
IF(${BUILD_UNITTESTS})
   ENABLE_TESTING()
   ADD_SUBDIRECTORY(unitTests)
   ADD_SUBDIRECTORY(benchmarkTests)
ENDIF()
