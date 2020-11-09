set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(cmake_and_compile.sh cmake_and_compile.sh COPYONLY)
configure_file(run.sh run.sh COPYONLY)
configure_file(pennyShow.sh pennyShow.sh COPYONLY)
configure_file(show_events.sh show_events.sh COPYONLY)
configure_file(show_geometry.sh show_geometry.sh COPYONLY)

install(DIRECTORY . DESTINATION useful_commands
       USE_SOURCE_PERMISSIONS
       PATTERN *~ EXCLUDE
       PATTERN "CMakeLists.txt" EXCLUDE
       PATTERN "cmake_and_compile.sh" EXCLUDE
       )
