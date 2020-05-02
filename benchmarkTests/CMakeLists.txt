#------------------------------------------------------------------------------------#
file(GLOB SOURCES *.cpp)
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_BenchTest_SOURCES "${SOURCES}")
#------------------------------------------------------------------------------------#
MACRO_APPEND_GLOBAL_VARIABLE(EVAA_BenchTest_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR};${CMAKE_CURRENT_BINARY_DIR}")
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
get_property(EVAA_BenchTest_INCLUDES GLOBAL PROPERTY EVAA_BenchTest_INCLUDES)
get_property(EVAA_BenchTest_SOURCES GLOBAL PROPERTY EVAA_BenchTest_SOURCES)
#------------------------------------------------------------------------------------#
# Set GBench 
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Suppressing benchmark's tests" FORCE)
add_subdirectory("C:/software/libs/GBENCH" ${CMAKE_BINARY_DIR}/gbench)
include_directories(${EVAA_INCLUDES})
include_directories(${EVAA_BenchTest_INCLUDES})
#MultiThreaded Compile with -MT 
#MultiThreadedDLL Compile with -MD
IF (CMAKE_SYSTEM_NAME MATCHES "Windows") 
	set_property(TARGET benchmark PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
ENDIF()
#------------------------------------------------------------------------------------#
add_executable(EVAABench  ${EVAA_BenchTest_SOURCES})
#------------------------------------------------------------------------------------#
IF (${USE_HDF5})
	target_link_libraries(EVAABench EVAALibStatic benchmark ${EVAA_3rd_LIBS} debug ${HDF5_LIBRARIES_DEBUG} optimized ${HDF5_LIBRARIES})
ELSE()
	target_link_libraries(EVAABench EVAALibStatic benchmark ${EVAA_3rd_LIBS})
ENDIF()
add_dependencies(EVAABench EVAALibStatic)
IF (CMAKE_SYSTEM_NAME MATCHES "Windows") 
	set_property(TARGET EVAABench PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
ENDIF()
#------------------------------------------------------------------------------------#
MESSAGE("==== Benchmark Test ====")
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
