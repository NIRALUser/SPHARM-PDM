
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
    ARCHIVE DESTINATION lib COMPONENT Development
    LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
    RUNTIME DESTINATION bin COMPONENT RuntimeLibraries
    )
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/SpatialObject COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
