
set(ShapeIO_SRCS
  SphericalHarmonicCoefficientFileReader.cxx
  SphericalHarmonicCoefficientFileWriter.cxx
  )

add_library(ShapeIO STATIC ${ShapeIO_SRCS})
target_link_libraries(ShapeIO PUBLIC SpatialObject)

if(NOT SPHARM-PDM_BUILD_SLICER_EXTENSION)
  install(TARGETS ShapeIO
    ARCHIVE DESTINATION lib COMPONENT Development
    LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
    RUNTIME DESTINATION bin COMPONENT RuntimeLibraries
    )
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/IO COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
