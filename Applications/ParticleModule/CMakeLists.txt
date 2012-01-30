#project(ParticleModule)

set(ParticleModule_SOURCE
    ParticleModuleCLP.cxx
    ParticleModuleParameters.cxx
    ParticleModuleParameters.h
    ParticleModuleComputation.cxx
    ParticleModuleComputation.h   )

GENERATECLP(ParticleModule_SOURCE ParticleModule.xml)
add_executable(ParticleModule ${ParticleModule_SOURCE})

target_link_libraries(ParticleModule
  itksys
  vtkIO
  ${ITK_LIBRARIES}
  ${BatchMake_LIBRARIES} )

install(TARGETS
  ParticleModule
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
