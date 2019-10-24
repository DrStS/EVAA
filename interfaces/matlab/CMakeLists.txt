#------------------------------------------------------------------------------------#
get_property(EVAA_INCLUDES GLOBAL PROPERTY EVAA_INCLUDES) 
include_directories(${EVAA_INCLUDES})
#------------------------------------------------------------------------------------#
matlab_add_mex(
    NAME EVAAMatlab
    SHARED
    MatlabInterface.cpp
    LINK_TO EVAALibStatic ${EVAA_3rd_LIBS}
    R2018a
)
IF (CMAKE_SYSTEM_NAME MATCHES "Windows") 
set_property(TARGET EVAAMatlab PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
ENDIF()
#------------------------------------------------------------------------------------#
MESSAGE("Matlab MEX")
MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
MESSAGE("EVAA_SOURCES: ${EVAA_SOURCES}")
#------------------------------------------------------------------------------------#
# Where to put EXECUTABLE
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin64)
#------------------------------------------------------------------------------------#