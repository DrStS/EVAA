#------------------------------------------------------------------------------------#
get_property(EVAA_INCLUDES GLOBAL PROPERTY EVAA_INCLUDES) 
include_directories(${EVAA_INCLUDES})
#------------------------------------------------------------------------------------#
matlab_add_mex(
    NAME EVAAMatlab
    SHARED
    MatlabInterface.cpp
    LINK_TO EVAAlib ${EVAA_3rd_LIBS}
    R2018a
)
#------------------------------------------------------------------------------------#
MESSAGE("Matlab MEX")
MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
MESSAGE("EVAA_SOURCES: ${EVAA_SOURCES}")
#------------------------------------------------------------------------------------#