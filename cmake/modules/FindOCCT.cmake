# A cmake module to find Opencascade
#------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------#
set(OCCTROOT_DIR $ENV{CASROOT})
IF (NOT OCCTROOT_DIR)
   IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
       IF (EXISTS "/home/stefan/software/tools/opencascade-7.0.0/opencascade-7.0.0")
          set(OCCTROOT_DIR "/home/stefan/software/tools/opencascade-7.0.0/opencascade-7.0.0")
       ENDIF(EXISTS "/home/stefan/software/tools/opencascade-7.0.0/opencascade-7.0.0")
   ENDIF()
   IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
       IF (EXISTS "C:/software/libs/OpenCASCADE7.0.0-vc12-64/opencascade-7.0.0")
          set(OCCTROOT_DIR "C:/software/libs/OpenCASCADE7.0.0-vc12-64/opencascade-7.0.0")
       ENDIF()
	   IF (EXISTS "C:/software/libs/OpenCASCADE7.1.0-vc12-64/opencascade-7.1.0")
          set(OCCTROOT_DIR "C:/software/libs/OpenCASCADE7.1.0-vc12-64/opencascade-7.1.0")
       ENDIF()
   ENDIF()
ENDIF ()
message("OCCTROOT_DIR is: ${OCCTROOT_DIR}")
#------------------------------------------------------------------------------------#
# Stage 1: find the include directory
#------------------------------------------------------------------------------------#
IF (NOT OCCT_INCLUDE_DIR)
    find_path(OCCT_INCLUDE_DIR
    TopoDS_Shape.hxx
    HINTS ${OCCTROOT_DIR}
    PATH_SUFFIXES inc
    )  
	IF(OCCT_INCLUDE_DIR MATCHES "MKL_INCLUDE_DIR-NOTFOUND")
     unset(OCCT_INCLUDE_DIR CACHE)
	ENDIF()
ENDIF ()
message("OCCT_INCLUDE_DIR is: ${OCCT_INCLUDE_DIR}")
#------------------------------------------------------------------------------------#
# Stage 2: find the lib directory
#------------------------------------------------------------------------------------#	
if (NOT OCCT_LIB_DIR)
  if (OCCTROOT_DIR)
    IF (OS_LINUX_X86_64)
	   set(EXPECT_OCCT_LIBPATH "${OCCTROOT_DIR}/linux64/gcc/lib/lib/")
    endif ()	
    IF (OS_WIN_X86_64)
	   set(EXPECT_OCCT_LIBPATH "${OCCTROOT_DIR}/win64/vc14/lib")
    endif ()	
    if (IS_DIRECTORY ${EXPECT_OCCT_LIBPATH})
      set(OCCT_LIB_DIR ${EXPECT_OCCT_LIBPATH})
    endif ()
  endif ()
endif ()
message("OCCT_LIB_DIR is: ${OCCT_LIB_DIR}")
#------------------------------------------------------------------------------------#
# Stage 3: find the libraries
#------------------------------------------------------------------------------------#	
IF (CMAKE_SYSTEM_NAME MATCHES "Linux")
#set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
ENDIF()
IF (CMAKE_SYSTEM_NAME MATCHES "Windows")
#set (CMAKE_FIND_LIBRARY_SUFFIXES .dll)
# lib needs to be used even for dynamic linking on windows 
# http://stackoverflow.com/questions/3250467/what-is-inside-lib-file-of-static-library-statically-linked-dynamic-library-an
ENDIF()
if (OCCT_LIB_DIR)
#find_library(OCCT_LIBRARY_FWOSPlugin        FWOSPlugin                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBin             TKBin                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBinL            TKBinL                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBinTObj         TKBinTObj                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBinXCAF         TKBinXCAF                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBO              TKBO                            ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBool            TKBool                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKBRep            TKBRep                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKCAF             TKCAF                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKCDF             TKCDF                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKD3DHost         TKD3DHost                       ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKDCAF            TKDCAF                          ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKDraw            TKDraw                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKernel           TKernel                         ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKFeat            TKFeat                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKFillet          TKFillet                        ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKG2d             TKG2d                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKG3d             TKG3d                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKGeomAlgo        TKGeomAlgo                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKGeomBase        TKGeomBase                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKHLR             TKHLR                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKIGES            TKIGES                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKIVtk            TKIVtk                          ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKIVtkDraw        TKIVtkDraw                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKLCAF            TKLCAF                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKMath            TKMath                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKMesh            TKMesh                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKMeshVS          TKMeshVS                        ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKOffset          TKOffset                        ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKOpenGl          TKOpenGl                        ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKPrim            TKPrim                          ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKQADraw          TKQADraw                        ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKService         TKService                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKShHealing       TKShHealing                     ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKStd             TKStd                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKStdL            TKStdL                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKSTEP            TKSTEP                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKSTEP209         TKSTEP209                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKSTEPAttr        TKSTEPAttr                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKSTEPBase        TKSTEPBase                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKSTL             TKSTL                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKTObj            TKTObj                          ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKTObjDRAW        TKTObjDRAW                      ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKTopAlgo         TKTopAlgo                       ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKTopTest         TKTopTest                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKV3d             TKV3d                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKVCAF            TKVCAF                          ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKViewerTest      TKViewerTest                    ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKVRML            TKVRML                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXCAF            TKXCAF                          ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKXDEDRAW         TKXDEDRAW                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXDEIGES         TKXDEIGES                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXDESTEP         TKXDESTEP                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXMesh           TKXMesh                         ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXml             TKXml                           ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXmlL            TKXmlL                          ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXmlTObj         TKXmlTObj                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXmlXCAF         TKXmlXCAF                       ${OCCT_LIB_DIR})
find_library(OCCT_LIBRARY_TKXSBase          TKXSBase                        ${OCCT_LIB_DIR})
#find_library(OCCT_LIBRARY_TKXSDRAW          TKXSDRAW                        ${OCCT_LIB_DIR})

IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set (OCCT_LIBRARIES  ${OCCT_LIBRARY_TKBin} ${OCCT_LIBRARY_TKBinL} ${OCCT_LIBRARY_TKBinTObj} ${OCCT_LIBRARY_TKBinXCAF} ${OCCT_LIBRARY_TKBO} ${OCCT_LIBRARY_TKBool} ${OCCT_LIBRARY_TKBRep} ${OCCT_LIBRARY_TKCAF} ${OCCT_LIBRARY_TKCDF} ${OCCT_LIBRARY_TKernel} ${OCCT_LIBRARY_TKFeat} ${OCCT_LIBRARY_TKFillet} ${OCCT_LIBRARY_TKG2d} ${OCCT_LIBRARY_TKG3d} ${OCCT_LIBRARY_TKGeomAlgo} ${OCCT_LIBRARY_TKGeomBase} ${OCCT_LIBRARY_TKHLR} ${OCCT_LIBRARY_TKIGES}  ${OCCT_LIBRARY_TKLCAF} ${OCCT_LIBRARY_TKMath} ${OCCT_LIBRARY_TKMesh} ${OCCT_LIBRARY_TKMeshVS} ${OCCT_LIBRARY_TKOffset} ${OCCT_LIBRARY_TKOpenGl} ${OCCT_LIBRARY_TKPrim} ${OCCT_LIBRARY_TKService} ${OCCT_LIBRARY_TKShHealing} ${OCCT_LIBRARY_TKStd} ${OCCT_LIBRARY_TKStdL} ${OCCT_LIBRARY_TKSTEP} ${OCCT_LIBRARY_TKSTEP209} ${OCCT_LIBRARY_TKSTEPAttr} ${OCCT_LIBRARY_TKSTEPBase} ${OCCT_LIBRARY_TKSTL} ${OCCT_LIBRARY_TKTObj} ${OCCT_LIBRARY_TKTopAlgo} ${OCCT_LIBRARY_TKV3d} ${OCCT_LIBRARY_TKVCAF} ${OCCT_LIBRARY_TKVRML} ${OCCT_LIBRARY_TKXCAF} ${OCCT_LIBRARY_TKXDEIGES} ${OCCT_LIBRARY_TKXDESTEP} ${OCCT_LIBRARY_TKXMesh} ${OCCT_LIBRARY_TKXml} ${OCCT_LIBRARY_TKXmlL} ${OCCT_LIBRARY_TKXmlTObj} ${OCCT_LIBRARY_TKXmlXCAF} ${OCCT_LIBRARY_TKXSBase} )
  ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows")
    set (OCCT_LIBRARIES  ${OCCT_LIBRARY_TKBin} ${OCCT_LIBRARY_TKBinL} ${OCCT_LIBRARY_TKBinTObj} ${OCCT_LIBRARY_TKBinXCAF} ${OCCT_LIBRARY_TKBO} ${OCCT_LIBRARY_TKBool} ${OCCT_LIBRARY_TKBRep} ${OCCT_LIBRARY_TKCAF} ${OCCT_LIBRARY_TKCDF} ${OCCT_LIBRARY_TKD3DHost} ${OCCT_LIBRARY_TKernel} ${OCCT_LIBRARY_TKFeat} ${OCCT_LIBRARY_TKFillet} ${OCCT_LIBRARY_TKG2d} ${OCCT_LIBRARY_TKG3d} ${OCCT_LIBRARY_TKGeomAlgo} ${OCCT_LIBRARY_TKGeomBase} ${OCCT_LIBRARY_TKHLR} ${OCCT_LIBRARY_TKIGES} ${OCCT_LIBRARY_TKIVtk} ${OCCT_LIBRARY_TKLCAF} ${OCCT_LIBRARY_TKMath} ${OCCT_LIBRARY_TKMesh} ${OCCT_LIBRARY_TKMeshVS} ${OCCT_LIBRARY_TKOffset} ${OCCT_LIBRARY_TKOpenGl} ${OCCT_LIBRARY_TKPrim} ${OCCT_LIBRARY_TKService} ${OCCT_LIBRARY_TKShHealing} ${OCCT_LIBRARY_TKStd} ${OCCT_LIBRARY_TKStdL} ${OCCT_LIBRARY_TKSTEP} ${OCCT_LIBRARY_TKSTEP209} ${OCCT_LIBRARY_TKSTEPAttr} ${OCCT_LIBRARY_TKSTEPBase} ${OCCT_LIBRARY_TKSTL} ${OCCT_LIBRARY_TKTObj} ${OCCT_LIBRARY_TKTopAlgo} ${OCCT_LIBRARY_TKV3d} ${OCCT_LIBRARY_TKVCAF} ${OCCT_LIBRARY_TKVRML} ${OCCT_LIBRARY_TKXCAF} ${OCCT_LIBRARY_TKXDEDRAW} ${OCCT_LIBRARY_TKXDEIGES} ${OCCT_LIBRARY_TKXDESTEP} ${OCCT_LIBRARY_TKXMesh} ${OCCT_LIBRARY_TKXml} ${OCCT_LIBRARY_TKXmlL} ${OCCT_LIBRARY_TKXmlTObj} ${OCCT_LIBRARY_TKXmlXCAF} ${OCCT_LIBRARY_TKXSBase})
 ELSE()
    set (OCCT_LIBRARIES "")
  ENDIF()
ENDIF()
# message("OCCT_LIBRARIES are: ${OCCT_LIBRARIES}")
# set OCCT_FOUND
include(FindPackageHandleStandardArgs)
IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
find_package_handle_standard_args(OCCT REQUIRED_VARS  OCCT_LIBRARY_TKBin OCCT_LIBRARY_TKBinL OCCT_LIBRARY_TKBinTObj OCCT_LIBRARY_TKBinXCAF OCCT_LIBRARY_TKBO OCCT_LIBRARY_TKBool OCCT_LIBRARY_TKBRep OCCT_LIBRARY_TKCAF OCCT_LIBRARY_TKCDF OCCT_LIBRARY_TKernel OCCT_LIBRARY_TKFeat OCCT_LIBRARY_TKFillet OCCT_LIBRARY_TKG2d OCCT_LIBRARY_TKG3d OCCT_LIBRARY_TKGeomAlgo OCCT_LIBRARY_TKGeomBase OCCT_LIBRARY_TKHLR OCCT_LIBRARY_TKIGES OCCT_LIBRARY_TKLCAF OCCT_LIBRARY_TKMath OCCT_LIBRARY_TKMesh OCCT_LIBRARY_TKMeshVS OCCT_LIBRARY_TKOffset OCCT_LIBRARY_TKOpenGl OCCT_LIBRARY_TKPrim OCCT_LIBRARY_TKService OCCT_LIBRARY_TKShHealing OCCT_LIBRARY_TKStd OCCT_LIBRARY_TKStdL OCCT_LIBRARY_TKSTEP OCCT_LIBRARY_TKSTEP209 OCCT_LIBRARY_TKSTEPAttr OCCT_LIBRARY_TKSTEPBase OCCT_LIBRARY_TKSTL OCCT_LIBRARY_TKTObj OCCT_LIBRARY_TKTopAlgo OCCT_LIBRARY_TKV3d OCCT_LIBRARY_TKVCAF OCCT_LIBRARY_TKVRML OCCT_LIBRARY_TKXCAF  OCCT_LIBRARY_TKXDEIGES OCCT_LIBRARY_TKXDESTEP OCCT_LIBRARY_TKXMesh OCCT_LIBRARY_TKXml OCCT_LIBRARY_TKXmlL OCCT_LIBRARY_TKXmlTObj OCCT_LIBRARY_TKXmlXCAF OCCT_LIBRARY_TKXSBase)
ELSEIF (CMAKE_SYSTEM_NAME MATCHES "Windows")
find_package_handle_standard_args(OCCT REQUIRED_VARS  OCCT_LIBRARY_TKBin OCCT_LIBRARY_TKBinL OCCT_LIBRARY_TKBinTObj OCCT_LIBRARY_TKBinXCAF OCCT_LIBRARY_TKBO OCCT_LIBRARY_TKBool OCCT_LIBRARY_TKBRep OCCT_LIBRARY_TKCAF OCCT_LIBRARY_TKCDF OCCT_LIBRARY_TKD3DHost  OCCT_LIBRARY_TKernel OCCT_LIBRARY_TKFeat OCCT_LIBRARY_TKFillet OCCT_LIBRARY_TKG2d OCCT_LIBRARY_TKG3d OCCT_LIBRARY_TKGeomAlgo OCCT_LIBRARY_TKGeomBase OCCT_LIBRARY_TKHLR OCCT_LIBRARY_TKIGES OCCT_LIBRARY_TKIVtk  OCCT_LIBRARY_TKLCAF OCCT_LIBRARY_TKMath OCCT_LIBRARY_TKMesh OCCT_LIBRARY_TKMeshVS OCCT_LIBRARY_TKOffset OCCT_LIBRARY_TKOpenGl OCCT_LIBRARY_TKPrim  OCCT_LIBRARY_TKService OCCT_LIBRARY_TKShHealing OCCT_LIBRARY_TKStd OCCT_LIBRARY_TKStdL OCCT_LIBRARY_TKSTEP OCCT_LIBRARY_TKSTEP209 OCCT_LIBRARY_TKSTEPAttr OCCT_LIBRARY_TKSTEPBase OCCT_LIBRARY_TKSTL OCCT_LIBRARY_TKTObj  OCCT_LIBRARY_TKTopAlgo OCCT_LIBRARY_TKV3d OCCT_LIBRARY_TKVCAF OCCT_LIBRARY_TKVRML OCCT_LIBRARY_TKXCAF  OCCT_LIBRARY_TKXDEIGES OCCT_LIBRARY_TKXDESTEP OCCT_LIBRARY_TKXMesh OCCT_LIBRARY_TKXml OCCT_LIBRARY_TKXmlL OCCT_LIBRARY_TKXmlTObj OCCT_LIBRARY_TKXmlXCAF OCCT_LIBRARY_TKXSBase)
ENDIF()