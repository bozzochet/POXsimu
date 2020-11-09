set(CMAKE_CXX_FLAGS "${GGS_CXX_FLAGS} ${CMAKE_CXX_FLAGS}")
set_property(DIRECTORY ${CMAKE_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES ${include_dirs})

configure_file(DetectorConstruction.hh DetectorConstruction.hh COPYONLY)

install(FILES
  DetectorConstruction.hh
  DESTINATION "${INSTALL_INC_DIR}"
  )
