#------------------------------------------------------------------------------------#
IF(OS_WIN_X86_64)
	SET(CODESYNTHESIS_PATH "https://codesynthesis.com/download/xsd/4.0/windows/i686/xsd-4.0.0-i686-windows.zip")
	SET(CODESYNTHESIS_FILE "xsd-4.0.0-i686-windows")
	SET(CODESYNTHESIS_FILE_END ".zip")
	SET(CODESYNTHESIS_FILE_HASH "d129469e109784c663387ca8bee5ac627434cfca")
ELSEIF(OS_LINUX_X86_64)
	SET(CODESYNTHESIS_PATH "https://codesynthesis.com/download/xsd/4.0/linux-gnu/x86_64/xsd-4.0.0-x86_64-linux-gnu.tar.bz2")
	SET(CODESYNTHESIS_FILE "sd-4.0.0-x86_64-linux-gnu")
	SET(CODESYNTHESIS_FILE_END ".tar.bz2")
	SET(CODESYNTHESIS_FILE_HASH "5eeb2eeca0d893949e3677bb374e7b96f19770d6")
ELSEIF(OS_MAC_X86_64)
	SET(CODESYNTHESIS_PATH "https://codesynthesis.com/download/xsd/4.0/macosx/i686/xsd-4.0.0-i686-macosx.tar.bz2")
	SET(CODESYNTHESIS_FILE "xsd-4.0.0-i686-macosx")
	SET(CODESYNTHESIS_FILE_END ".tar.bz2")
	SET(CODESYNTHESIS_FILE_HASH "bb96454da6acafb93180368220d555e2b9747023")
ELSE() 
   MESSAGE(FATAL_ERROR "EVAA is supported on Linux x86_64 and Windows x64 and Mac x86_64 only")   
ENDIF()
#------------------------------------------------------------------------------------#
IF(CMAKE_OUTPUT_VERBOSE)
	MESSAGE("Downloading... ${EVAA_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END}")
ENDIF()
IF(NOT EXISTS ${EVAA_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END})
FILE(DOWNLOAD ${CODESYNTHESIS_PATH} ${EVAA_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END}
	TIMEOUT 60
	EXPECTED_HASH SHA1=${CODESYNTHESIS_FILE_HASH}
	TLS_VERIFY ON)
ENDIF()
FILE(ARCHIVE_EXTRACT INPUT ${EVAA_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}${CODESYNTHESIS_FILE_END} DESTINATION ${EVAA_SOURCE_DIR}/xsd)
FIND_PROGRAM(XSD_EXE xsd PATHS ${EVAA_SOURCE_DIR}/xsd/${CODESYNTHESIS_FILE}/bin REQUIRED NO_DEFAULT_PATH)
IF(CMAKE_OUTPUT_VERBOSE)
	MESSAGE("XSD_EXE ${XSD_EXE}")
ENDIF()
#https://www.codesynthesis.com/projects/xsd/documentation/xsd.xhtml
EXECUTE_PROCESS(  
	COMMAND ${XSD_EXE} cxx-tree --type-naming ucc  --function-naming lcc --hxx-suffix .h --cxx-suffix .cpp ${EVAA_SOURCE_DIR}/xsd/IP_EVAA_XML.xsd
	WORKING_DIRECTORY ${EVAA_SOURCE_DIR}/src/database
   )
#------------------------------------------------------------------------------------#