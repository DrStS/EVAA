#  Copyright &copy; 2019, Dr. Stefan Sicklinger, Munich
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
SET(CMAKE_VERBOSE_MAKEFILE off)
cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0091 NEW)#MSVC runtime library flags are selected by an abstraction
#------------------------------------------------------------------------------------#
project (EVAA VERSION 0.1 LANGUAGES CXX)
set(CMAKE_CONFIGURATION_TYPES "Release;RelWithDebInfo")
#------------------------------------------------------------------------------------#)
# Check OS and host system
IF((${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "x86_64")) 
   SET(OS_LINUX_X86_64 TRUE)
   MESSAGE("EVAA is going to be configured for Linux x86_64!")
ELSEIF((${CMAKE_SYSTEM_NAME} MATCHES "Windows" AND ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "AMD64")) 
   SET(OS_WIN_X86_64 TRUE)
   MESSAGE("EVAA is going to be configured for Windows x64!")  
ELSE() 
   MESSAGE(FATAL_ERROR "EVAA is supported on Linux x86_64 and Windows x64 only")   
ENDIF()
# Diag output
MESSAGE("Current source dir is" ${EVAA_SOURCE_DIR})
#------------------------------------------------------------------------------------#
# 1. Options
# 2. Macros
# 3. Compiler
# 4. Subdirectories
# 5. Configure
# 6. MPI
#------------------------------------------------------------------------------------#
######################################################################################
# THESE GLOBAL VARIABLES ARE SET DURING THE CONFIGURATION
######################################################################################
#EVAA_thirdparty_SOURCES
#EVAA_thirdparty_INCLUDES
#EVAA_SOURCES
#EVAA_INCLUDES
######################################################################################
#1. Option
######################################################################################
option(USE_INTEL_MKL           "Use Intel MKL"  ON)
option(USE_EIGEN               "Use Eigen"  ON)
option(USE_BLAZE               "Use Blaze"  ON)
option(USE_XSD                 "Use XSD with XERCES"  OFF)
option(USE_HDF5                "Use HDF5"  OFF)
option(USE_MATLAB              "Build Matlab API"  OFF)
option(EVAA_COMMANDLINE        "Build EVAA COMMANDLINE"  ON)
option(BUILD_UNITTEST          "Build EVAA UNITEST"  OFF)
option(BUILD_BENCHMARK         "Build EVAA MICRO BENCHMARKS"  ON)
######################################################################################
#2. Macros
######################################################################################
# Add EVAA cmake modules
SET(CMAKE_MODULE_PATH ${EVAA_SOURCE_DIR}/cmake/modules)
include(MacroOutOfSourceBuild)
include(MacroAppendGlobalVariable)
include(ExternalProject)
include(GetGitRevisionDescription)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
######################################################################################
#3. Compiler
######################################################################################
IF (${EVAA_COMMANDLINE})
   MESSAGE("Build EVAA commandline only!")
   add_definitions(-DEVAA_COMMANDLINE_ON)
   SET(USE_OCCT  false)
   SET(USE_VTK   false)
   SET(USE_QT5   false)
ELSE ()
   MESSAGE("Build EVAA with GUI QT Quick!")
   SET(USE_OCCT  false)
   SET(USE_VTK   true)
   SET(USE_QT5   true)
ENDIF()
#------------------------------------------------------------------------------------#
# Avoid in source builds
macro_ensure_out_of_source_build("In source build is not allowed for ${PROJECT_NAME}. Please run CMake from your build directory. You may want to remove CMakeCache.txt.")
#------------------------------------------------------------------------------------#
# Check compilers, their versions
IF ( ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel") )
   SET(USE_INTEL_COMPILER ON)
   MESSAGE("C++ compiler is ${CMAKE_CXX_COMPILER_ID} and has version ${CMAKE_CXX_COMPILER_VERSION} " ) 
ELSEIF ( ( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU") )
   SET(USE_GCC_COMPILER ON)
   MESSAGE("C++ compiler is ${CMAKE_CXX_COMPILER_ID} and has version ${CMAKE_CXX_COMPILER_VERSION} " ) 
ELSEIF ( ( ${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC ) )
   SET(USE_MICROSOFT_COMPILER ON)
   MESSAGE("C++ compiler is ${CMAKE_CXX_COMPILER_ID} and has version ${CMAKE_CXX_COMPILER_VERSION} " ) 
ELSE()
  MESSAGE(FATAL_ERROR "Your CXX compiler is ${CMAKE_CXX_COMPILER_ID} which is not supported by the build system yet") 
ENDIF() 
#------------------------------------------------------------------------------------#
# Try to find Opencascade
IF (${USE_OCCT})
   find_package(OCCT REQUIRED)
   add_definitions(-DOCCT_ON)
   MESSAGE("OCCT status: ${OCCT_FOUND}")
   IF (${USE_MICROSOFT_COMPILER})
      add_definitions(-DWNT)
   ENDIF()
   IF (${USE_GCC_COMPILER})
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -lglut -lGL -lGLU -lX11")
   ENDIF()
   include_directories(${OCCT_INCLUDE_DIR})
   SET(EVAA_3rd_LIBS ${EVAA_3rd_LIBS}  ${OCCT_LIBRARIES})
ENDIF()
#------------------------------------------------------------------------------------#
# Try to find VTK
IF (${USE_VTK})
 add_definitions(-DUSE_VTK)
 SET(VTK_BASE_DIR "C:/software/libs/VTK-8.2.0")
 IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
   SET(VTK_DIR "${VTK_BASE_DIR}/build")
 ENDIF()
 IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
     SET(VTK_DIR "/opt/software/libs/VTK-8.1.1/gcc/lib/cmake/vtk-8.1")
 ENDIF()
  MESSAGE("VTK_DIR: ${VTK_DIR}")
  find_package(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})
  LIST(REMOVE_ITEM VTK_INCLUDE_DIRS ${VTK_BASE_DIR}/ThirdParty/hdf5)
  LIST(REMOVE_ITEM VTK_INCLUDE_DIRS ${VTK_BASE_DIR}/build/ThirdParty/hdf5)
  LIST(REMOVE_ITEM VTK_INCLUDE_DIRS ${VTK_BASE_DIR}/build/ThirdParty/hdf5/vtkhdf5)
  LIST(REMOVE_ITEM VTK_INCLUDE_DIRS ${VTK_BASE_DIR}/ThirdParty/hdf5/vtkhdf5/hl/src)
  LIST(REMOVE_ITEM VTK_INCLUDE_DIRS ${VTK_BASE_DIR}/ThirdParty/hdf5/vtkhdf5/src)
  MESSAGE("VTK status: ${VTK_FOUND} ")#${VTK_INCLUDE_DIRS}
  include_directories(${VTK_INCLUDE_DIR})
  SET(EVAA_3rd_LIBS ${EVAA_3rd_LIBS}  ${VTK_LIBRARIES})
ENDIF()
#------------------------------------------------------------------------------------#
# Try to find Intel MKL
SET(CMAKE_EXE_LINKER_FLAGS   "")
IF (${USE_INTEL_MKL})
  SET(MKL_THREADING_LAYER "sequential")
  IF(MKL_THREADING_LAYER MATCHES "sequential")
  # add_definitions(-DMKL_DIRECT_CALL_SEQ_JIT)
  ENDIF()
  IF (MKL_THREADING_LAYER MATCHES "thread")
  # add_definitions(-DMKL_DIRECT_CALL_JIT)
  ENDIF()
  find_package(MKL REQUIRED)
   #add_definitions(-DMKL_ILP64)
   add_definitions(-DUSE_INTEL_MKL)
MESSAGE("MKL status: ${MKL_FOUND}")
   # Get OpenMP activated and use Intel link library
   IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
      IF(USE_INTEL_COMPILER)
         SET(CMAKE_C_FLAGS   "/MD /Qopenmp")
         SET(CMAKE_CXX_FLAGS "/MD /Qopenmp")
      ELSE()
         SET(CMAKE_C_FLAGS   "/MT /openmp")
         SET(CMAKE_CXX_FLAGS "/MT /openmp")
         SET(CMAKE_EXE_LINKER_FLAGS "/nodefaultlib:vcomp ${MKL_OMP_LIBRARY}")
      ENDIF()
   ENDIF()
   IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
      SET(CMAKE_C_FLAGS    "-fopenmp")
      SET(CMAKE_CXX_FLAGS  "-fopenmp")
      SET(CMAKE_EXE_LINKER_FLAGS "-L${MKL_OMP_LIBRARY}")
   ENDIF()
   include_directories(${MKL_INCLUDE_DIR})
   SET(EVAA_3rd_LIBS ${EVAA_3rd_LIBS} ${MKL_LIBRARIES})
ENDIF()
#------------------------------------------------------------------------------------#
IF (${USE_EIGEN}) 
   add_definitions(-DUSE_EIGEN)
   include_directories("C:/software/libs/EIGEN")
ENDIF()
#------------------------------------------------------------------------------------#
IF (${USE_BLAZE}) 
   add_definitions(-DUSE_BLAZE)
   include_directories("C:/software/libs/BLAZE")
ENDIF()
#------------------------------------------------------------------------------------#
IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
   SET(CMAKE_PREFIX_PATH "C:/software/libs/QT/Qt5123/5.12.3/msvc2017_64/lib/cmake")
ENDIF()
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
   SET(CMAKE_PREFIX_PATH "/opt/software/libs/Qt/5.12.3/gcc_64/lib/cmake")
ENDIF()
# Try to find QT5 Quick
IF (${USE_QT5}) 
   add_definitions(-DUSE_QT5)
   find_package(Qt5 REQUIRED Core Gui Widgets Qml Quick QuickControls2)
   MESSAGE("Qt5 status: ${Qt5_FOUND}")
   SET(CMAKE_AUTOMOC ON)
   SET(CMAKE_AUTOUIC ON)
   SET(CMAKE_AUTORCC ON)
   SET(CMAKE_INCLUDE_CURRENT_DIR ON)
ENDIF()
#------------------------------------------------------------------------------------#
# Try to find XSD with XERCES
IF (${USE_XSD})
   find_package(XSDXERCES REQUIRED)
   MESSAGE("XSDXERCES status: ${XSDXERCES_FOUND}")
   IF (${USE_MICROSOFT_COMPILERS_WINDOWS})
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      ENDIF()
   IF (${USE_GCC_COMPILERS_LINUX})
      SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
   ENDIF()
   include_directories(${XSD_INCLUDE_DIR})
   include_directories(${XERCES_INCLUDE_DIR})
   SET(EVAA_3rd_LIBS ${EVAA_3rd_LIBS} ${XERCES_LIBRARIES})
ENDIF()
#------------------------------------------------------------------------------------#
# Try to find HDF5
IF (${USE_HDF5})
   add_definitions(-DUSE_HDF5)
   find_package(HDF5 REQUIRED)
   MESSAGE("HDF5 status: ${HDF5_FOUND}")
   include_directories(${HDF5_INCLUDE_DIR})
ENDIF()
#------------------------------------------------------------------------------------#
# Try to find Matlab
IF (${USE_MATLAB})
   add_definitions(-DUSE_MATLAB)
   find_package(Matlab REQUIRED)
   MESSAGE("Matlab status: ${MATLAB_FOUND}")
ENDIF()
#------------------------------------------------------------------------------------#
# Global compiler flags
IF (${USE_INTEL_COMPILER})
   SET(CMAKE_C_STANDARD 99) 
   SET(CMAKE_C_STANDARD_REQUIRED ON)
   SET(CMAKE_C_EXTENSIONS OFF)
   SET(CMAKE_CXX_STANDARD 14) 
   SET(CMAKE_CXX_STANDARD_REQUIRED ON)
   SET(CMAKE_CXX_EXTENSIONS OFF)
   SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -O3 -fp-model fast -qopt-report=5 -xHost -D NOALIAS")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fp-model fast -qopt-report=5 -xHost -D NOALIAS")
ENDIF()
IF (${USE_MICROSOFT_COMPILER})
   SET(CMAKE_CXX_FLAGS_DEBUG     "${CMAKE_CXX_FLAGS_DEBUG} -DEBUG")
   SET(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE} ")
ENDIF()
######################################################################################
#4. Configure
######################################################################################
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
# INCLUDE dirs
add_subdirectory(src)
IF(${BUILD_BENCHMARK})
add_subdirectory(benchmarkTests)
ENDIF()
IF(${BUILD_UNITTEST})
add_subdirectory(unitTests)
ENDIF()
IF(${USE_MATLAB})
add_subdirectory(interfaces)
ENDIF()
