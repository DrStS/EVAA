#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
FILE(GLOB EVAA_src_subdirectories "*")
FOREACH(f ${EVAA_src_subdirectories})
   IF(NOT f MATCHES "\\.txt" AND IS_DIRECTORY ${f} )
      add_subdirectory(${f})
   ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#
get_property(EVAA_INCLUDES GLOBAL PROPERTY EVAA_INCLUDES) 
get_property(EVAA_SOURCES GLOBAL PROPERTY EVAA_SOURCES) 
#------------------------------------------------------------------------------------#
# SET OCCT
IF (${USE_OCCT})
   include_directories(${OCCT_INCLUDE_DIR})
   SET(EVAA_LIBS ${EVAA_LIBS}  ${OCCT_LIBRARIES})
ENDIF()
# Set HDF5
IF (${USE_HDF5})
   include_directories(${HDF5_INCLUDE_DIR})
   # SET(EVAA_LIBS ${EVAA_LIBS} optimized ${HDF5_LIBRARIES} debug ${HDF5_LIBRARIES_DEBUG})
ENDIF()
# SET VTK
IF (${USE_VTK})
   include_directories(${VTK_INCLUDE_DIR})
   SET(EVAA_LIBS ${EVAA_LIBS}  ${VTK_LIBRARIES})
ENDIF()
# SET Intel MKL 
IF (${USE_INTEL_MKL})
   include_directories(${MKL_INCLUDE_DIR})
   SET(EVAA_LIBS ${EVAA_LIBS} ${MKL_LIBRARIES})
ENDIF()
# SET QT5 
IF (${USE_QT5})
ENDIF()
# SET SIMULIA ODB 
IF (${USE_SIMULIA_ODB})
   include_directories(${SIMULIA_INCLUDE_DIR})
   SET(EVAA_LIBS ${EVAA_LIBS} ${SIMULIA_LIBRARIES})
ENDIF()
IF (${USE_XSD})
   include_directories(${XSD_INCLUDE_DIR})
   include_directories(${XERCES_INCLUDE_DIR})
   SET(EVAA_LIBS ${EVAA_LIBS} ${XERCES_LIBRARIES})
ENDIF()
#------------------------------------------------------------------------------------#
include_directories(${EVAA_INCLUDES})
LIST(REMOVE_ITEM EVAA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
IF (${EVAA_COMMANDLINE})
   LIST(REMOVE_ITEM EVAA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/EVAAMainWindow.h)
   LIST(REMOVE_ITEM EVAA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/EVAAMainWindow.cpp)
ENDIF()
#------------------------------------------------------------------------------------#
IF (CMAKE_SYSTEM_NAME MATCHES "Windows") 
# WIN32
   add_executable(EVAA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${EVAA_SOURCES} resources.qrc winResources.rc)
ENDIF()
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
   add_executable(EVAA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${EVAA_SOURCES} resources.qrc)
ENDIF()
target_link_libraries(EVAA ${EVAA_LIBS} debug ${HDF5_LIBRARIES_DEBUG} optimized ${HDF5_LIBRARIES})
IF (${USE_QT5})
   qt5_use_modules(EVAA Core Gui Widgets Qml Quick QuickControls2)
ENDIF()
#------------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------------#
MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
MESSAGE("EVAA_LIBS: ${EVAA_LIBS}")
MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
MESSAGE("EVAA_SOURCES: ${EVAA_SOURCES}")
MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
#------------------------------------------------------------------------------------#
# Where to put EXECUTABLE
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin64)
#------------------------------------------------------------------------------------#
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()
