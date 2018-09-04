
set(SpatialObject_SRCS
  SphericalHarmonicSpatialObject.cxx
  )

add_library(SpatialObject STATIC ${SpatialObject_SRCS})
target_link_libraries(SpatialObject
  PUBLIC
    ${ITK_LIBRARIES}
  PRIVATE
    ShapeAlgorithms
  )

if(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
  install(TARGETS SpatialObject
    ARCHIVE DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_ARCHIVE_DESTINATION} COMPONENT Development
    LIBRARY DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
    RUNTIME DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION} COMPONENT RuntimeLibraries
    )
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/SpatialObject COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
