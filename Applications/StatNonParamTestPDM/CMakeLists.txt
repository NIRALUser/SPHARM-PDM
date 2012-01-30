#project(shapeAnalysisMANCOVA)

   set(shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.cxx MANCOVA.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx  miscMANCOVA.cxx nonMANCOVA.cxx )
  GENERATECLP(shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.xml)

  add_executable(shapeAnalysisMANCOVA ${shapeAnalysisMANCOVA_SRCS})

  target_link_libraries(shapeAnalysisMANCOVA
    ${ITK_LIBRARIES}
    ${VTK_LIBRARIES}
    boost_date_time-mt
    boost_system-mt
    )

install(TARGETS
  shapeAnalysisMANCOVA
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

