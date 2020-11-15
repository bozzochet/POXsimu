set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(run_savewrl.mac run_savewrl.mac COPYONLY)
configure_file(run.mac run.mac COPYONLY)
configure_file(geo.mac geo.mac COPYONLY)
configure_file(geo_4x7.mac geo_4x7.mac COPYONLY)
configure_file(drawtgeometry.C drawtgeometry.C COPYONLY)
configure_file(rootlogon.C ../rootlogon.C COPYONLY)
configure_file(vis.mac ../vis.mac COPYONLY)

install(DIRECTORY . DESTINATION macros
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN rootlogon.C EXCLUDE
       PATTERN vis.mac EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       )
