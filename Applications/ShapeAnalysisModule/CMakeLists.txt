#project(ShapeAnalysisModule)

set(ShapeAnalysisModule_SOURCE
    ShapeAnalysisModuleCLP.cxx
    Parameters.cxx
    Parameters.h
    ShapeAnalysisModuleComputation.cxx
    ShapeAnalysisModuleComputation.h   )

GENERATECLP(ShapeAnalysisModule_SOURCE ShapeAnalysisModule.xml)
add_executable(ShapeAnalysisModule ${ShapeAnalysisModule_SOURCE})

target_link_libraries(ShapeAnalysisModule
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${BatchMake_LIBRARIES} )


install(TARGETS
  ShapeAnalysisModule
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
