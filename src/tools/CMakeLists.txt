#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
file(GLOB INCLUDES *.h LIST_DIRECTORIES false)
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_INCLUDES "${INCLUDES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
#------------------------------------------------------------------------------------#

#------------------------------------------------------------------------------------#
FILE(GLOB EVAA_src_sub_subdirectories "*")
FOREACH(sf ${EVAA_src_sub_subdirectories})
  IF(NOT sf MATCHES "\\.txt" AND IS_DIRECTORY ${sf} )
    add_subdirectory(${sf})
  ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#