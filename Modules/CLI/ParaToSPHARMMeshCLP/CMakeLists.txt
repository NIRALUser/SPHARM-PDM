SEMMacroBuildCLI(
  NAME ParaToSPHARMMeshCLP
  EXECUTABLE_ONLY
  ADDITIONAL_SRCS
    itkMeshTovtkPolyData.cxx
    vtkPolyDataToitkMesh.cxx
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
    ShapeIO
    SpatialObject
    ShapeAlgorithms
    SpatialObject
  )
