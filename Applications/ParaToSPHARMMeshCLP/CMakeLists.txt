PROJECT(ParaToSPHARMMeshCLP)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

SET (ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)

GENERATECLP(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.xml)

ADD_EXECUTABLE (ParaToSPHARMMeshCLP ${ParaToSPHARMMeshCLP_SOURCE})

TARGET_LINK_LIBRARIES ( ParaToSPHARMMeshCLP
       ${ITK_LIBRARIES}
       vtkIO vtkHybrid vtkCommon vtkFiltering 
	ShapeIO SpatialObject ShapeAlgorithms SpatialObject)

IF (Slicer3_FOUND)
 slicer3_set_plugins_output_path(ParaToSPHARMMeshCLP)
 slicer3_install_plugins(ParaToSPHARMMeshCLP)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  ParaToSPHARMMeshCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
