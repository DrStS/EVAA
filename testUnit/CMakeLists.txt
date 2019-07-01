#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_UnitTest_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_UnitTest_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
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
get_property(EVAA_UnitTest_INCLUDES GLOBAL PROPERTY EVAA_UnitTest_INCLUDES)
get_property(EVAA_UnitTest_SOURCES GLOBAL PROPERTY EVAA_UnitTest_SOURCES)
#------------------------------------------------------------------------------------#
# Set GTest
SET(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
add_subdirectory("C:/software/libs/GTEST" ${CMAKE_BINARY_DIR}/gtest)
enable_testing()
include_directories(${GTEST_INCLUDE_DIRS})
include_directories(${EVAA_INCLUDES})
include_directories(${EVAA_UnitTest_INCLUDES})
#------------------------------------------------------------------------------------#
add_executable(testEVAA ${EVAA_UnitTest_SOURCES})
#------------------------------------------------------------------------------------#
target_link_libraries(testEVAA EVAAlib gtest gtest_main ${EVAA_3rd_LIBS} debug ${HDF5_LIBRARIES_DEBUG} optimized ${HDF5_LIBRARIES})
add_test(AllTests testEVAA)
#------------------------------------------------------------------------------------#
MESSAGE("Unit Test")
MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
MESSAGE("EVAA_UnitTest_INCLUDES: ${EVAA_UnitTest_INCLUDES}")
MESSAGE("EVAA_UnitTest_SOURCES: ${EVAA_UnitTest_SOURCES}")
MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
#------------------------------------------------------------------------------------#
# Where to put EXECUTABLE
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin64)
#------------------------------------------------------------------------------------#