include(spharm_add_executable)

spharm_add_executable(
  NAME SpharmTool
  TARGET_LIBRARIES
    PRIVATE
      ${ITK_LIBRARIES}
      ShapeIO SpatialObject
      ShapeAlgorithms
      SpatialObject
  )
