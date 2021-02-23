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
FILE(GLOB EVAA_src_subdirectories "*")
FOREACH(f ${EVAA_src_subdirectories})
   IF(NOT f MATCHES "\\.txt" AND IS_DIRECTORY ${f} )
      add_subdirectory(${f})
   ENDIF()
ENDFOREACH()
#------------------------------------------------------------------------------------#
get_property(EVAA_INCLUDES GLOBAL PROPERTY EVAA_INCLUDES)
get_property(EVAA_INCLUDE_DIRS GLOBAL PROPERTY EVAA_INCLUDE_DIRS) 
get_property(EVAA_SOURCES GLOBAL PROPERTY EVAA_SOURCES)
get_property(EVAA_INCLUDE_DIRS_AUTOGEN GLOBAL PROPERTY EVAA_INCLUDE_DIRS_AUTOGEN)
get_property(EVAA_SOURCES_AUTOGEN GLOBAL PROPERTY EVAA_SOURCES_AUTOGEN)
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
ADD_LIBRARY(EVAALibStaticGenerated STATIC ${EVAA_SOURCES_AUTOGEN})
TARGET_INCLUDE_DIRECTORIES(EVAALibStaticGenerated PUBLIC ${EVAA_INCLUDE_DIRS_AUTOGEN} ${CODESYNTHESIS_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(EVAALibStaticGenerated CONAN_PKG::xerces-c)
#------------------------------------------------------------------------------------#
ADD_CUSTOM_TARGET(
   clang-tidy
   COMMAND clang-tidy-10 -p=${CMAKE_BINARY_DIR} ${EVAA_SOURCES})
FOREACH(header_file ${EVAA_INCLUDES})
   GET_FILENAME_COMPONENT(header_file_only "${header_file}" NAME)
   LIST(APPEND HEADERS_FOR_CLANG_TIDY "${header_file_only}")
   MESSAGE("${header_file_only}")
ENDFOREACH()
STRING(REPLACE ";" "|" HEADERS_FOR_CLANG_TIDY "${HEADERS_FOR_CLANG_TIDY}")
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/.clang-tidy.in" "${PROJECT_SOURCE_DIR}/.clang-tidy" @ONLY)

ADD_CUSTOM_TARGET(
   clang-format-diffprint
   COMMAND python3 ${PROJECT_SOURCE_DIR}/helperScripts/run-clang-format/run-clang-format.py --style=file  ${EVAA_SOURCES} ${EVAA_INCLUDES})

ADD_CUSTOM_TARGET(
   clang-format-apply
   COMMAND python3 ${PROJECT_SOURCE_DIR}/helperScripts/run-clang-format/run-clang-format.py -i --style=file  ${EVAA_SOURCES} ${EVAA_INCLUDES})

ADD_CUSTOM_TARGET(
   clang-format-check
   COMMAND python3 ${PROJECT_SOURCE_DIR}/helperScripts/run-clang-format/run-clang-format.py -q --style=file  ${EVAA_SOURCES} ${EVAA_INCLUDES})

if(CMAKE_EXPORT_COMPILE_COMMANDS)
    SET(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES} ${CONAN_INCLUDE_DIRS})
endif()
#------------------------------------------------------------------------------------#
LIST(REMOVE_ITEM EVAA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
ADD_LIBRARY(EVAALibStatic STATIC ${EVAA_SOURCES})
TARGET_INCLUDE_DIRECTORIES(EVAALibStatic PUBLIC ${EVAA_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(EVAALibStatic EVAALibStaticGenerated CONAN_PKG::hdf5 CONAN_PKG::boost)
#------------------------------------------------------------------------------------#
IF(OS_WIN_X86_64) 
   ADD_EXECUTABLE(EVAA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp winResources.rc)
ELSE()
   ADD_EXECUTABLE(EVAA ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp)
ENDIF()
#------------------------------------------------------------------------------------#
TARGET_LINK_LIBRARIES(EVAA EVAALibStatic)
#------------------------------------------------------------------------------------#
IF(CMAKE_OUTPUT_VERBOSE)
   MESSAGE("==== Executable ====")
   MESSAGE("C Compiler Flags are: ${CMAKE_C_FLAGS}")
   MESSAGE("C++ Compiler Flags are: ${CMAKE_CXX_FLAGS} ")
   MESSAGE("Linker Flags are: ${CMAKE_EXE_LINKER_FLAGS}")
   MESSAGE("EVAA_3rd_LIBS: ${EVAA_3rd_LIBS}")
   MESSAGE("EVAA_INCLUDES: ${EVAA_INCLUDES}")
   MESSAGE("EVAA_SOURCES: ${EVAA_SOURCES}")
ENDIF()
#------------------------------------------------------------------------------------#
SET(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)
#------------------------------------------------------------------------------------#
get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
   MESSAGE(STATUS "dir='${dir}'")
endforeach()