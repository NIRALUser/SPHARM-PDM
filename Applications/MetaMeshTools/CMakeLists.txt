#project(MetaMeshTools)

# Shorthand to source and binary directories
set(MY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Include directories
include_directories( ${MY_SOURCE_DIR} ${MY_BINARY_DIR})
link_directories( ${CMAKE_BINARY_SOURCE_DIR})

set(SRC_STL2META_CXX STL2Meta.cxx vtkPolyDataToitkMesh.cxx)
add_executable(STL2Meta ${SRC_STL2META_CXX})
target_link_libraries(STL2Meta ${ITK_LIBRARIES} ${VTK_LIBRARIES})


set(SRC_META2STL_CXX Meta2STL.cxx itkMeshTovtkPolyData.cxx)
add_executable(Meta2STL ${SRC_META2STL_CXX})
target_link_libraries(Meta2STL ${ITK_LIBRARIES} ${VTK_LIBRARIES})

#set(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
#add_executable(MeshMath ${SRC_MESHMATH_CXX})
#target_link_libraries(MeshMath ${ITK_LIBRARIES} ${VTK_LIBRARIES})
SET(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
ADD_EXECUTABLE(MeshMath ${SRC_MESHMATH_CXX})
TARGET_LINK_LIBRARIES(MeshMath ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(Meta2VTK Meta2VTK.cxx itkMeshTovtkPolyData.cxx)
target_link_libraries(Meta2VTK ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(VTK2Meta VTK2Meta.cxx vtkPolyDataToitkMesh.cxx)
target_link_libraries(VTK2Meta ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

ADD_EXECUTABLE(asc2vtk asc2vtk.cxx)
ADD_EXECUTABLE(asc2meta asc2meta.cxx)

ADD_EXECUTABLE(BYU2VTK BYU2VTK.cxx)
target_link_libraries(BYU2VTK ${VTK_LIBRARIES} )

install(TARGETS
  STL2Meta
  Meta2STL
  MeshMath
  Meta2VTK
  VTK2Meta
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
