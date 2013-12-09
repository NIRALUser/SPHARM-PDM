#project(GenParaMeshCLP)

set(GenParaMeshCLP_SOURCE GenParaMeshCLP.cxx itkMeshTovtkPolyData.cxx)

GENERATECLP(GenParaMeshCLP_SOURCE GenParaMeshCLP.xml)
add_executable(GenParaMeshCLP ${GenParaMeshCLP_SOURCE})
target_link_libraries( GenParaMeshCLP
  ${ITK_LIBRARIES}
  vtkIO ShapeAlgorithms SparseMatrixLib)

add_library(GenParaMeshCLPLib STATIC ${GenParaMeshCLP_SOURCE})
set_target_properties (GenParaMeshCLPLib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
target_link_libraries(GenParaMeshCLPLib
  ${ITK_LIBRARIES}
  vtkIO ShapeAlgorithms SparseMatrixLib
)

if( NOT EXTENSION_SUPERBUILD_BINARY_DIR )
 install(TARGETS
   GenParaMeshCLP
   BUNDLE DESTINATION  bin
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib/static  
 )
endif()
