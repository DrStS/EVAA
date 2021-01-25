#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_BENCHMARK_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_BENCHMARK_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
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
get_property(EVAA_BENCHMARK_INCLUDES GLOBAL PROPERTY EVAA_BENCHMARK_INCLUDES)
get_property(EVAA_BENCHMARK_SOURCES GLOBAL PROPERTY EVAA_BENCHMARK_SOURCES)
#------------------------------------------------------------------------------------#
INCLUDE_DIRECTORIES(${EVAA_INCLUDES})
INCLUDE_DIRECTORIES(${EVAA_BENCHMARK_INCLUDES})
#------------------------------------------------------------------------------------#
IF(OS_WIN_X86_64) 
   ADD_EXECUTABLE(EVAABench ${EVAA_BENCHMARK_SOURCES} winResourcesBenchmarkTest.rc)
ELSE()
   ADD_EXECUTABLE(EVAABench ${EVAA_BENCHMARK_SOURCES})
ENDIF()
#------------------------------------------------------------------------------------#
TARGET_LINK_LIBRARIES(EVAABench EVAALibStatic CONAN_PKG::benchmark)
#------------------------------------------------------------------------------------#
IF(CMAKE_OUTPUT_VERBOSE)
   MESSAGE("==== Benchmark Test Executable ====")
   MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
   MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
   MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
   MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
   MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
   MESSAGE("EVAA_BENCHMARK_INCLUDES: ${EVAA_BENCHMARK_INCLUDES}")
   MESSAGE("EVAA_BENCHMARK_SOURCES: ${EVAA_BENCHMARK_SOURCES}")
   MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
ENDIF()
#------------------------------------------------------------------------------------#
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
   MESSAGE(STATUS "dir='${dir}'")
endforeach()