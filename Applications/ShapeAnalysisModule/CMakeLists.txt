project(ShapeAnalysisModule)

set(ShapeAnalysisModule_SOURCE
    ShapeAnalysisModuleCLP.cxx
    Parameters.cxx
    Parameters.h
    ShapeAnalysisModuleComputation.cxx
    ShapeAnalysisModuleComputation.h 
)

GENERATECLP(ShapeAnalysisModule_SOURCE ShapeAnalysisModule.xml)
add_executable(ShapeAnalysisModule ${ShapeAnalysisModule_SOURCE})

target_link_libraries(ShapeAnalysisModule
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${BatchMake_LIBRARIES} )


if( NOT SPHARM-PDM_BUILD_SLICER_EXTENSION )
  install(TARGETS
    ShapeAnalysisModule
    BUNDLE DESTINATION  bin
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static  
  )
else()
  install(TARGETS ShapeAnalysisModule DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} )
endif()
