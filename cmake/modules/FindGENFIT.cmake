# find the GENFIT++ constrained fitter wrapper
# simplistic, should mark variables as cached or advanced

message(STATUS "Looking for GENFIT...")


set(GENFIT_SEARCH_PATHS
  "${CMAKE_SOURCE_DIR}/../GenFit"
  "${CMAKE_SOURCE_DIR}/../genfit"
)

MESSAGE(${GENFIT_SEARCH_PATHS})

find_path(GENFIT_BASE_DIR lgpl.txt
  PATHS  ${GENFIT_SEARCH_PATHS}  
  NO_DEFAULT_PATH
  )

set(GENFIT_INCLUDE_DIR ${GENFIT_BASE_DIR}/core/include
                       ${GENFIT_BASE_DIR}/eventDisplay/include/
                       ${GENFIT_BASE_DIR}/fields/include/
                       ${GENFIT_BASE_DIR}/finitePlanes/include/
                       ${GENFIT_BASE_DIR}/fitters/include/
                       ${GENFIT_BASE_DIR}/GBL/include
                       ${GENFIT_BASE_DIR}/GFRave/include
                       ${GENFIT_BASE_DIR}/measurements/include
                       ${GENFIT_BASE_DIR}/trackReps/include
                       ${GENFIT_BASE_DIR}/utilities/include
                       )

message(STATUS "Dirs: ${GENFIT_INCLUDE_DIR}")

if(NOT GENFIT_INCLUDE_DIR)
  Message(STATUS "Looking for GENFIT... - GENFIT.hpp not found")
  if(GENFIT_FIND_REQUIRED)
    message(FATAL_ERROR "GENFIT is required, please make sure cmake finds it")
  endif()
  return()
endif()

include_directories(${GENFIT_INCLUDE_DIR})

FIND_LIBRARY(GENFIT_LIBRARIES NAMES genfit2
  PATHS ${GENFIT_SEARCH_PATHS}
  PATH_SUFFIXES "build/lib"
  NO_DEFAULT_PATH
  )

SET(GENFIT_LIBRARY_DIR "${GENFIT_BASE_DIR}/build/lib")

set(GENFIT_FOUND TRUE)
Message(STATUS "Looking for GENFIT... - Found ${GENFIT_LIBRARIES}")

