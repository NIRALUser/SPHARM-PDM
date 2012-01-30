#project(GenParaMeshCLP)

set(GenParaMeshCLP_SOURCE GenParaMeshCLP.cxx itkMeshTovtkPolyData.cxx)

GENERATECLP(GenParaMeshCLP_SOURCE GenParaMeshCLP.xml)
add_executable(GenParaMeshCLP ${GenParaMeshCLP_SOURCE})
target_link_libraries( GenParaMeshCLP
  ${ITK_LIBRARIES}
  vtkIO ShapeAlgorithms SparseMatrixLib)


install(TARGETS
  GenParaMeshCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
