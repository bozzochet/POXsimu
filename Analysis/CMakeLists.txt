set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})
#include_directories(${PROJECT_SOURCE_DIR}/include)

# declare the project name
PROJECT(Analysis)

configure_file(DAQ.C DAQ.C COPYONLY)
configure_file(Analysis.C Analysis.C COPYONLY)
configure_file(magneticField.cxx magneticField.cxx COPYONLY)

install(DIRECTORY . DESTINATION Analysis
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       )

###################################################################

FIND_PACKAGE(Eigen3 REQUIRED)

###################################################################

find_package(ROOT CONFIG REQUIRED COMPONENTS EG Geom)

include( ${ROOT_USE_FILE} )
#include("${ROOT_DIR}/modules/RootNewMacros.cmake")

#Install the ROOT 6 PCM and rootmap files.
if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      DESTINATION "${INSTALL_LIB_DIR}")
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

# declare the ROOT library path
LINK_DIRECTORIES(${ROOT_LIBRARY_DIR})

###################################################################

find_path(
	GENFIT_INCLUDE_DIR
	NAMES ProlateSpacepointMeasurement.h
	PATH_SUFFIXES include
	PATHS $ENV{GENFIT_DIR}
	)
if(NOT GENFIT_INCLUDE_DIR)
    message(FATAL_ERROR "GenFit headers not found!")
endif()

include_directories(${GENFIT_INCLUDE_DIR})

find_library(GENFIT_LIBRARY
	NAMES genfit2
	PATH_SUFFIXES lib
	PATHS $ENV{GENFIT_DIR}
	)
if(NOT GENFIT_LIBRARY)
    message(FATAL_ERROR "GenFit was not found!")
endif()

get_filename_component(GENFIT_LIBRARY_DIR ${GENFIT_LIBRARY} DIRECTORY)
link_directories(${GENFIT_LIBRARY_DIR})

###################################################################

INCLUDE_DIRECTORIES(
                BEFORE
                ${ROOT_INCLUDE_DIRS}
                SYSTEM ${EIGEN3_INCLUDE_DIR}
                ${GF_INC_DIRS}
		)

###################################################################

set(GenFitMinimalExample_SRC
	GenFitMinimalExample.cc
)

add_executable(GenFitMinimalExample ${GenFitMinimalExample_SRC})
target_link_libraries(GenFitMinimalExample
	${GENFIT_LIBRARY}
	${ROOT_LIBRARIES}
)

set_target_properties(GenFitMinimalExample PROPERTIES
				     RUNTIME_OUTPUT_DIRECTORY "${BIN_DIR}"
				     )

install(TARGETS GenFitMinimalExample DESTINATION "${INSTALL_BIN_DIR}")
