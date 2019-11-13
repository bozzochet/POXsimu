include_directories(${PROJECT_SOURCE_DIR}/include)

####################### Build HerdBaselineGeometry shared library ##################
# Set source

set(TestGeometry_SOURCE
  DetectorConstruction.cc
#  DetectorMessenger.cc
)

# Set target
add_library(TestGeometry SHARED ${TestGeometry_SOURCE} )

target_link_libraries(TestGeometry
                      ${Geant4_LIBRARIES}
                      ${GGS_GEO_LIBRARIES}
                     )

set_target_properties(TestGeometry PROPERTIES VERSION 1)

# Set install rules
install(TARGETS TestGeometry DESTINATION plugins)
