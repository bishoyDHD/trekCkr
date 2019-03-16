# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION 2.4.3 FATAL_ERROR)
# Needed for cmake 2.6
IF(COMMAND cmake_policy)
	cmake_policy(SET CMP0003 NEW)
ENDIF(COMMAND cmake_policy)
##############################
# POI:
# Set project name
PROJECT(CsI)

#################################
# POI:
# Set the library version
SET(CSI_MAJOR_VERSION 0)
SET(CSI_MINOR_VERSION 1)
SET(CSI_PATCH_VERSION 0)
SET(CSI_VERSION "${CSI_MAJOR_VERSION}.${CSI_MINOR_VERSION}.${CSI_PATCH_VERSION}")
LIST(APPEND CSI_LIBRARY_PROPERTIES VERSION "${CSI_VERSION}" SOVERSION "${CSI_MAJOR_VERSION}")

###################################################
# POI:
# find source and header files, set paths, includes...
SET(LIBNAME "CsI")

FILE(GLOB SRCS src/*.cpp)

FILE(GLOB HEADERS src/include/*.h ${CMAKE_SOURCE_DIR}/include/cookerrawtree.h)

SET(INCLUDES ${ROOT_INCLUDE_DIR}  ${GSL_INCLUDE_DIR} "${CMAKE_CURRENT_SOURCE_DIR}/src/include" "${CMAKE_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES(${INCLUDES} )

# create root dictionary
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/src/include/LinkDef.hh)
SET(DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
LIST(APPEND SRCS ${DICTIONARY})

###################################
# POI:
# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS})
TARGET_LINK_LIBRARIES(${LIBNAME}  cookerrawtree Plugin "-L${ROOT_LIBRARY_DIR}" ${ROOT_LIBRARIES})
TARGET_COMPILE_FEATURES(${LIBNAME} PRIVATE cxx_auto_type cxx_range_for)
SET_TARGET_PROPERTIES(${LIBNAME} PROPERTIES ${CSI_LIBRARY_PROPERTIES})



# install the library
INSTALL(TARGETS ${LIBNAME} RUNTIME DESTINATION ${BINDIR} LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})

INSTALL (FILES init/demo.xml DESTINATION "../shared/init")
INSTALL (FILES init/demo.xsd DESTINATION "../shared/init")
