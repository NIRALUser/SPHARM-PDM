#project(ParaToSPHARMMeshCLP)

set(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)

GENERATECLP(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.xml)

add_executable(ParaToSPHARMMeshCLP ${ParaToSPHARMMeshCLP_SOURCE})

target_link_libraries( ParaToSPHARMMeshCLP
       ${ITK_LIBRARIES}
       ${VTK_LIBRARIES}
       ShapeIO SpatialObject ShapeAlgorithms SpatialObject
       ${CLAPACK_LIBRARIES}
       )

install(TARGETS
  ParaToSPHARMMeshCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
