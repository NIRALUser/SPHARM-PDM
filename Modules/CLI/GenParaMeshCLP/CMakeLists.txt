SEMMacroBuildCLI(
  NAME GenParaMeshCLP
  EXECUTABLE_ONLY
  ADDITIONAL_SRCS
    itkMeshTovtkPolyData.cxx
    vtkPolyDataToitkMesh.cxx
  TARGET_LIBRARIES
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
    ShapeAlgorithms
    SparseMatrixLib
  )
