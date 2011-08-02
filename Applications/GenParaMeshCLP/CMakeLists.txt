project(GenParaMeshCLP)

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

  find_package(GenerateCLP REQUIRED)
  if(GenerateCLP_FOUND)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_FOUND)
    message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  endif(GenerateCLP_FOUND)

endif(Slicer3_FOUND)


set(GenParaMeshCLP_SOURCE GenParaMeshCLP.cxx itkMeshTovtkPolyData.cxx)

GENERATECLP(GenParaMeshCLP_SOURCE GenParaMeshCLP.xml)
add_executable(GenParaMeshCLP ${GenParaMeshCLP_SOURCE})
target_link_libraries( GenParaMeshCLP
  ${ITK_LIBRARIES}
  vtkIO ShapeAlgorithms SparseMatrixLib)


if(Slicer3_FOUND)
 slicer3_set_plugins_output_path(GenParaMeshCLP)
 slicer3_install_plugins(GenParaMeshCLP)
endif(Slicer3_FOUND)

install(TARGETS
  GenParaMeshCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
