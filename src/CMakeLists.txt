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
include_directories(${EVAA_INCLUDES})
LIST(REMOVE_ITEM EVAA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
add_library(EVAALibStatic STATIC ${EVAA_SOURCES})
#------------------------------------------------------------------------------------#
IF (CMAKE_SYSTEM_NAME MATCHES "Windows") 
   #WIN32
   add_executable(EVAA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp winResources.rc)
ENDIF()
#------------------------------------------------------------------------------------#
target_link_libraries(EVAA EVAALibStatic)
#------------------------------------------------------------------------------------#
IF(CMAKE_OUTPUT_VERBOSE)
MESSAGE("==== Executable ====")
MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
MESSAGE("EVAA_SOURCES: ${EVAA_SOURCES}")
MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
ENDIF()
#------------------------------------------------------------------------------------#
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
#------------------------------------------------------------------------------------#
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()