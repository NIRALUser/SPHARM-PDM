SET (ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
GENERATECLP(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.xml)
ADD_EXECUTABLE (ParaToSPHARMMeshCLP ${ParaToSPHARMMeshCLP_SOURCE})
TARGET_LINK_LIBRARIES ( ParaToSPHARMMeshCLP
       ${ITK_LIBRARIES}
       vtkIO vtkHybrid vtkCommon vtkFiltering 
	ShapeIO SpatialObject ShapeAlgorithms SpatialObject)



