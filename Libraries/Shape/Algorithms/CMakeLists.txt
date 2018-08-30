
set(ShapeAlgorithms_SRCS
  EqualAreaParametricMeshNewtonIterator.cxx
  ParametricMeshToSPHARMSpatialObjectFilter.cxx
  SphericalHarmonicMedialAxisMeshSource.cxx
  SphericalHarmonicMeshSource.cxx
  TSurfaceNet.cxx
  TVoxelVolume.cxx
  )

add_library(ShapeAlgorithms STATIC ${ShapeAlgorithms_SRCS})
target_link_libraries(ShapeAlgorithms
  PUBLIC
    ShapeNumerics
    SpatialObject
    ${ITK_LIBRARIES}
  PRIVATE
    SparseMatrixLib
    ${CLAPACK_LIBRARIES}
  )

if(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
  install(TARGETS ShapeAlgorithms
    ARCHIVE DESTINATION lib COMPONENT Development
    LIBRARY DESTINATION lib COMPONENT RuntimeLibraries
    RUNTIME DESTINATION bin COMPONENT RuntimeLibraries
    )
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/Algorithms COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
