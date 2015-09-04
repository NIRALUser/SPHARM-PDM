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

if( NOT SPHARM-PDM_BUILD_SLICER_EXTENSION )
  install(TARGETS
    ParaToSPHARMMeshCLP
    BUNDLE DESTINATION  bin
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static  
  )
else()
  install(TARGETS ParaToSPHARMMeshCLP DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} )
endif()
