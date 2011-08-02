project(ParaToSPHARMMeshCLP)

cmake_minimum_required(VERSION 2.6)

set(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)

GENERATECLP(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.xml)

add_executable(ParaToSPHARMMeshCLP ${ParaToSPHARMMeshCLP_SOURCE})

target_link_libraries( ParaToSPHARMMeshCLP
       ${ITK_LIBRARIES}
       vtkIO vtkHybrid vtkCommon vtkFiltering
       ShapeIO SpatialObject ShapeAlgorithms SpatialObject
       lapack
       )

if(Slicer3_FOUND)
 slicer3_set_plugins_output_path(ParaToSPHARMMeshCLP)
 slicer3_install_plugins(ParaToSPHARMMeshCLP)
endif(Slicer3_FOUND)

install(TARGETS
  ParaToSPHARMMeshCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
