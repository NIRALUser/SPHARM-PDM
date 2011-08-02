project(SegPostProcessCLP)

cmake_minimum_required(VERSION 2.6)

find_package(Slicer3 QUIET NO_DEFAULT_PATH)

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

  find_package(GenerateCLP QUIET)
  if(GenerateCLP_FOUND)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_FOUND)
    message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  endif(GenerateCLP_FOUND)
endif(Slicer3_FOUND)

set(SegPostProcessCLP_SOURCE SegPostProcessCLP.cxx )
GENERATECLP(SegPostProcessCLP_SOURCE SegPostProcessCLP.xml)
add_executable(SegPostProcessCLP ${SegPostProcessCLP_SOURCE})

target_link_libraries(SegPostProcessCLP ${ITK_LIBRARIES})

if(Slicer3_FOUND)
  slicer3_set_plugins_output_path(SegPostProcessCLP)
  slicer3_install_plugins(SegPostProcessCLP)
endif(Slicer3_FOUND)

install(TARGETS
  SegPostProcessCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
