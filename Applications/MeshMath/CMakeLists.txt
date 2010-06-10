

# Shorthand to source and binary directories
SET(MY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Include directories
INCLUDE_DIRECTORIES( ${MY_SOURCE_DIR} ${MY_BINARY_DIR})

LINK_DIRECTORIES( ${CMAKE_BINARY_SOURCE_DIR})
LINK_LIBRARIES( 
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkHybrid
  vtkFiltering 
)



SET(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
ADD_EXECUTABLE(MeshMath ${SRC_MESHMATH_CXX})
TARGET_LINK_LIBRARIES(MeshMath ${ITK_LIBRARIES})

