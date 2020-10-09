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
#Catch 2
include(CTest)
include(${CONAN_BUILD_DIRS_CATCH2}/Catch.cmake)
include_directories(${EVAA_INCLUDES})
include_directories(${EVAA_UnitTest_INCLUDES})
#------------------------------------------------------------------------------------#
IF(OS_WIN_X86_64) 
   ADD_EXECUTABLE(EVAATest ${EVAA_UnitTest_SOURCES} winResourcesUnitTest.rc)
ELSE()
   ADD_EXECUTABLE(EVAATest ${EVAA_UnitTest_SOURCES})
ENDIF()
#------------------------------------------------------------------------------------#
target_link_libraries(EVAATest EVAALibStatic CONAN_PKG::catch2)
catch_discover_tests(EVAATest)
ADD_CUSTOM_COMMAND(TARGET EVAATest POST_BUILD
COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${PROJECT_SOURCE_DIR}/inputFiles/inputFileLinear.xml"
        $<TARGET_FILE_DIR:EVAATest>)
ADD_CUSTOM_COMMAND(TARGET EVAATest POST_BUILD		
COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${PROJECT_SOURCE_DIR}/xsd/InputSchemaEVAA.xsd"
        $<TARGET_FILE_DIR:EVAATest>)
#------------------------------------------------------------------------------------#
IF(CMAKE_OUTPUT_VERBOSE)
MESSAGE("==== Unit Test Executable ====")
MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
MESSAGE("EVAA_UnitTest_INCLUDES: ${EVAA_UnitTest_INCLUDES}")
MESSAGE("EVAA_UnitTest_SOURCES: ${EVAA_UnitTest_SOURCES}")
MESSAGE("MKL_LIBRARIES: ${MKL_LIBRARIES}")
ENDIF()
#------------------------------------------------------------------------------------#
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()