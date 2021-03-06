# A cmake module to find XSD and XERCES
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#

set(XSDROOT_DIR $ENV{XSDROOT})
set(XERCESROOT_DIR $ENV{XERCESROOT})
IF (NOT XSDROOT_DIR)
   IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
       IF (EXISTS "/opt/software/libs/XSD/libxsd")
          set(XSDROOT_DIR "/opt/software/libs/XSD/libxsd")
       ENDIF(EXISTS "/opt/software/libs/XSD/libxsd")
   ENDIF()
   IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
       IF (EXISTS "C:/software/libs/XSD/libxsd")
          set(XSDROOT_DIR "C:/software/libs/XSD/libxsd")
       ENDIF()
   ENDIF()
ENDIF ()
message("XSDROOT_DIR is: ${XSDROOT_DIR}")
IF (NOT XERCESROOT_DIR)
   IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
       IF (EXISTS "/opt/software/libs/XERCES/xerces-c")
          set(XERCESROOT_DIR "/opt/software/libs/XERCES/xerces-c")
       ENDIF(EXISTS "/opt/software/libs/XERCES/xerces-c")
   ENDIF()
   IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
       IF (EXISTS "C:/software/libs/XERCES/xerces-c")
          set(XERCESROOT_DIR "C:/software/libs/XERCES/xerces-c")
       ENDIF()
   ENDIF()
ENDIF ()
message("XERCESROOT_DIR is: ${XERCESROOT_DIR}")
#------------------------------------------------------------------------------------#
# Stage 1: find the include directory
#------------------------------------------------------------------------------------#
IF (NOT XSD_INCLUDE_DIR)
set(EXPECT_XSD_INCLUDE_DIR "${XSDROOT_DIR}")	
    if (IS_DIRECTORY ${EXPECT_XSD_INCLUDE_DIR})
      set(XSD_INCLUDE_DIR ${EXPECT_XSD_INCLUDE_DIR})
    endif ()
	IF(XSD_INCLUDE_DIR MATCHES "XSD_INCLUDE_DIR-NOTFOUND")
         UNSET(XSD_INCLUDE_DIR CACHE)
	ENDIF()
ENDIF ()
message("XSD_INCLUDE_DIR is: ${XSD_INCLUDE_DIR}")
IF (NOT XERCES_INCLUDE_DIR)
  set(EXPECT_XERCES_INCLUDE_DIR "${XERCESROOT_DIR}/include")	
    if (IS_DIRECTORY ${EXPECT_XERCES_INCLUDE_DIR})
      set(XERCES_INCLUDE_DIR ${EXPECT_XERCES_INCLUDE_DIR})
    endif ()
	IF(XERCES_INCLUDE_DIR MATCHES "XERCES_INCLUDE_DIR-NOTFOUND")
         UNSET(XERCES_INCLUDE_DIR CACHE)
	ENDIF()
ENDIF ()
message("XERCES_INCLUDE_DIR is: ${XERCES_INCLUDE_DIR}")
#------------------------------------------------------------------------------------#
# Stage 2: find the lib directory
#------------------------------------------------------------------------------------#	
if (NOT XERCES_LIB_DIR)
  if (XERCESROOT_DIR)
        IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
	 SET(EXPECT_XERCES_LIBPATH "${XERCESROOT_DIR}/lib")
        ENDIF()
	IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
         SET(EXPECT_XERCES_LIBPATH "${XERCESROOT_DIR}/lib")
        ENDIF()
    if (IS_DIRECTORY ${EXPECT_XERCES_LIBPATH})
      SET(XERCES_LIB_DIR ${EXPECT_XERCES_LIBPATH})
    endif ()
  endif ()
endif ()
message("XERCES_LIB_DIR is: ${XERCES_LIB_DIR}")
#------------------------------------------------------------------------------------#
# Stage 3: find the libraries
#------------------------------------------------------------------------------------#	
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
set (CMAKE_FIND_LIBRARY_SUFFIXES .so)
set (CMAKE_FIND_LIBRARY_PREFIXES lib)
ENDIF()
IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
set (CMAKE_FIND_LIBRARY_SUFFIXES .lib)
ENDIF()
if (XERCES_LIB_DIR)
    IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
      find_library(XERCES_XERCES-C_LIBRARY     xerces-c     ${XERCES_LIB_DIR})
    ENDIF()
    IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
      find_library(XERCES_XERCES-C_LIBRARY     xerces-c_3     ${XERCES_LIB_DIR})
    ENDIF()
  if (XERCES_XERCES-C_LIBRARY)
    IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
      set (XERCES_LIBRARIES "${XERCES_XERCES-C_LIBRARY}")
    ENDIF()
    IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
      set (XERCES_LIBRARIES ${XERCES_XERCES-C_LIBRARY} )
    ENDIF()
  else ()
    set (XERCES_LIBRARIES "")
  endif ()
endif ()
get_filename_component(XSDEXE_DIR ${XSDROOT_DIR} DIRECTORY)
SET(XSDEXE_DIR "${XSDEXE_DIR}/bin")
FIND_PROGRAM(XSD_EXE xsd PATHS ${XSDEXE_DIR})
# set XSDXERCES_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(XSDXERCES REQUIRED_VARS  XERCES_XERCES-C_LIBRARY)
mark_as_advanced(XSD_INCLUDE_DIR XERCES_INCLUDE_DIR XERCES_LIB_DIR XERCES_LIBRARIES)

