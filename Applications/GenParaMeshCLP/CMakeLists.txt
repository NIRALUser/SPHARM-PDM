
SET (GenParaMeshCLP_SOURCE GenParaMeshCLP.cxx itkMeshTovtkPolyData.cxx)
GENERATECLP(GenParaMeshCLP_SOURCE GenParaMeshCLP.xml)
ADD_EXECUTABLE (GenParaMeshCLP ${GenParaMeshCLP_SOURCE})
TARGET_LINK_LIBRARIES ( GenParaMeshCLP
		   ${ITK_LIBRARIES}
		   vtkIO SparseMatrixLib)


