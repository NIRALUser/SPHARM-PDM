
set(ShapeIO_SRCS
  SphericalHarmonicCoefficientFileReader.cxx
  SphericalHarmonicCoefficientFileWriter.cxx
  )

add_library(ShapeIO STATIC ${ShapeIO_SRCS})
target_link_libraries(ShapeIO PUBLIC SpatialObject)

if(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
  install(TARGETS ShapeIO
    ARCHIVE DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_ARCHIVE_DESTINATION} COMPONENT Development
    LIBRARY DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
    RUNTIME DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION} COMPONENT RuntimeLibraries
    )
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/IO COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
