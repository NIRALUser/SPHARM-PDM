PROJECT(ParaToSPHARMMeshCLP)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(Slicer3 QUIET NO_DEFAULT_PATH)

IF (Slicer3_FOUND)
  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

ELSE (Slicer3_FOUND)
  FIND_PACKAGE(ITK REQUIRED)
  IF (ITK_FOUND)
    INCLUDE(${USE_ITK_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  ENDIF (ITK_FOUND)

  FIND_PACKAGE(VTK REQUIRED)
  IF (VTK_FOUND)
    INCLUDE(${USE_VTK_FILE})
  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
  ENDIF (VTK_FOUND)

  FIND_PACKAGE(GenerateCLP REQUIRED)
  IF(GenerateCLP_FOUND)
    INCLUDE(${GenerateCLP_USE_FILE})
  ELSE(GenerateCLP_FOUND)
    MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  ENDIF(GenerateCLP_FOUND)

ENDIF (Slicer3_FOUND)



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
