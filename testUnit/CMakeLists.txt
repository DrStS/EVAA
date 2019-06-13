#------------------------------------------------------------------------------------#
add_subdirectory("C:/software/libs/GTEST" ${CMAKE_BINARY_DIR}/gtest)
enable_testing()
include_directories(${GTEST_INCLUDE_DIRS})
#------------------------------------------------------------------------------------#
add_executable(testEVAA testEVAA.cpp)
#------------------------------------------------------------------------------------#
set (CMAKE_CXX_FLAGS_RELEASE "/MT")
set( CMAKE_CXX_FLAGS_DEBUG "/MTd")
target_link_libraries(testEVAA gtest gtest_main)
add_test(AllTests testEVAA)
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