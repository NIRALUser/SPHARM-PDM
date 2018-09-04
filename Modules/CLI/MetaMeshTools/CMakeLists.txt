include(spharm_add_executable)

# STL -> Meta
spharm_add_executable(
  NAME STL2Meta
  ADDITIONAL_SRCS
    vtkPolyDataToitkMesh.cxx
  TARGET_LIBRARIES PRIVATE
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
  )

# Meta -> STL
spharm_add_executable(
  NAME Meta2STL
  ADDITIONAL_SRCS
    itkMeshTovtkPolyData.cxx
  TARGET_LIBRARIES PRIVATE
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
  )

# MeshMath
spharm_add_executable(
  NAME MeshMath
  ADDITIONAL_SRCS
    itkMeshTovtkPolyData.cxx
    vtkPolyDataToitkMesh.cxx
  TARGET_LIBRARIES PRIVATE
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
    ShapeAlgorithms
  )

# Meta -> VTK
spharm_add_executable(
  NAME Meta2VTK
  ADDITIONAL_SRCS
    itkMeshTovtkPolyData.cxx
  TARGET_LIBRARIES PRIVATE
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
  )

# VTK -> Meta
spharm_add_executable(
  NAME VTK2Meta
  ADDITIONAL_SRCS
    vtkPolyDataToitkMesh.cxx
  TARGET_LIBRARIES PRIVATE
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
  )

# asc -> VTK
spharm_add_executable(NAME asc2vtk)

# asc -> Meta
spharm_add_executable(NAME asc2meta)

# BYU -> VTK
spharm_add_executable(
  NAME BYU2VTK
  TARGET_LIBRARIES PRIVATE
    ${VTK_LIBRARIES}
  )

# STL -> VTK
spharm_add_executable(
  NAME STL2VTK
  TARGET_LIBRARIES PRIVATE
    ${VTK_LIBRARIES}
  )
