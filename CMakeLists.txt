project(TestGeometry)

cmake_minimum_required(VERSION 3.15)

SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY TRUE)

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(GGS 2.3.1 REQUIRED)
set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
# Remove /usr/include from list of include folders.
# This folder is added by Geant4 through GGS, and it causes this error:
#   #error "You need a ISO C conforming compiler to use the glibc headers"
# when generating Root dictionaries in a system with a recent version of glibc
# (tested with glibc 2.23).
get_property(include_dirs DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
list(REMOVE_ITEM include_dirs "/usr/include")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

add_subdirectory(src)

# Set install rules
install(DIRECTORY macros DESTINATION .
	PATTERN *~ EXCLUDE
	PATTERN rootlogon.C EXCLUDE
	PATTERN vis.mac EXCLUDE)

install(DIRECTORY Analysis DESTINATION .
	PATTERN *~ EXCLUDE)

install(FILES macros/rootlogon.C DESTINATION .)
install(FILES macros/vis.mac DESTINATION .)
