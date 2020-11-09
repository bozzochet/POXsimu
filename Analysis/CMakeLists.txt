set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(DAQ.C DAQ.C COPYONLY)
configure_file(Analysis.C Analysis.C COPYONLY)
configure_file(magneticField.cxx magneticField.cxx COPYONLY)

install(DIRECTORY . DESTINATION Analysis
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       )
