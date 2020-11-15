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

set_target_properties(TestGeometry PROPERTIES
                      VERSION 1
		      LIBRARY_OUTPUT_DIRECTORY "${LIB_DIR}"
		      )

if (APPLE)
#    set_property(TARGET TestGeometry PROPERTY PREFIX "lib")
#    set_property(TARGET TestGeometry PROPERTY OUTPUT_NAME "TestGeometry.so")
    set_property(TARGET TestGeometry PROPERTY SUFFIX ".so")
#    set_property(TARGET TestGeometry PROPERTY SOVERSION "32.1.2.0")
endif()

# Set install rules
install(TARGETS TestGeometry DESTINATION "${INSTALL_LIB_DIR}")
