project(MetaMeshTools)

cmake_minimum_required(VERSION 2.6)

find_package(Slicer3 QUIET NO_DEFAULT_PATH)

if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

else(Slicer3_FOUND)
  find_package(ITK REQUIRED)
  if(ITK_FOUND)
    include(${USE_ITK_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  endif(ITK_FOUND)

  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    include(${USE_VTK_FILE})
  else(VTK_FOUND)
    message(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
  endif(VTK_FOUND)

endif(Slicer3_FOUND)






# Shorthand to source and binary directories
set(MY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(MY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Include directories
include_directories( ${MY_SOURCE_DIR} ${MY_BINARY_DIR})

link_directories( ${CMAKE_BINARY_SOURCE_DIR})
link_libraries(
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkHybrid
  vtkFiltering
)

set(SRC_STL2META_CXX STL2Meta.cxx vtkPolyDataToitkMesh.cxx)
add_executable(STL2Meta ${SRC_STL2META_CXX})
target_link_libraries(STL2Meta ${ITK_LIBRARIES})


set(SRC_META2STL_CXX Meta2STL.cxx itkMeshTovtkPolyData.cxx)
add_executable(Meta2STL ${SRC_META2STL_CXX})
target_link_libraries(Meta2STL ${ITK_LIBRARIES})

set(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
add_executable(MeshMath ${SRC_MESHMATH_CXX})
target_link_libraries(MeshMath ${ITK_LIBRARIES})


add_executable(Meta2VTK Meta2VTK.cxx itkMeshTovtkPolyData.cxx)
target_link_libraries(Meta2VTK ${ITK_LIBRARIES})

add_executable(VTK2Meta VTK2Meta.cxx vtkPolyDataToitkMesh.cxx)
target_link_libraries(VTK2Meta ${ITK_LIBRARIES} )

if(Slicer3_FOUND)
  slicer3_set_plugins_output_path(STL2Meta)
  slicer3_install_plugins(STL2Meta)
  slicer3_set_plugins_output_path(Meta2STL)
  slicer3_install_plugins(Meta2STL)
  slicer3_set_plugins_output_path(MeshMath)
  slicer3_install_plugins(MeshMath)
  slicer3_set_plugins_output_path(Meta2VTK)
  slicer3_install_plugins(Meta2VTK)
  slicer3_set_plugins_output_path(VTK2Meta)
  slicer3_install_plugins(VTK2Meta)
endif(Slicer3_FOUND)

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
