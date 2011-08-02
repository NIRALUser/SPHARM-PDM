project(ParticleModule)

cmake_minimum_required(VERSION 2.6)

if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
else(Slicer3_FOUND)

  find_package(ITK REQUIRED)
  if(ITK_FOUND)
    include(${USE_ITK_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  endif(ITK_FOUND)

  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    include(${USE_VTK_FILE})
  else(VTK_FOUND)
    message(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
  endif(VTK_FOUND)

  find_package(GenerateCLP REQUIRED)
  if(GenerateCLP_FOUND)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_FOUND)
    message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  endif(GenerateCLP_FOUND)

endif(Slicer3_FOUND)

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

if(Slicer3_FOUND)
  slicer3_set_plugins_output_path(ParticleModule)
  slicer3_install_plugins(ParticleModule)
endif(Slicer3_FOUND)

install(TARGETS
  ParticleModule
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
