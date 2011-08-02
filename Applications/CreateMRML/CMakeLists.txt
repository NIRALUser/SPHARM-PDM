project(CreateMRML)

cmake_minimum_required(VERSION 2.6)

find_package(Slicer3 REQUIRED)
if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
endif(Slicer3_FOUND)

set(CreateMRML_SOURCE
  CreateMRML.cxx
  MRMLColorableHelper.h
  MRMLColorableHelper.cxx
  MRMLFiducialHelper.cxx
  MRMLFiducialHelper.h
  MRMLNodeHelper.h
  MRMLModelHelper.h
  MRMLTransformHelper.h
  MRMLVolumeHelper.h
  CreateMRMLSceneHelper.h
  CreateMRMLSceneHelper.cxx
  MagicParser.h
  MagicParser.cxx
)

add_executable(CreateMRML ${CreateMRML_SOURCE})
target_link_libraries(CreateMRML ${Slicer3_Libs_LIBRARIES} )

if(Slicer3_FOUND)
  slicer3_set_plugins_output_path(CreateMRML)
  slicer3_install_plugins(CreateMRML)
endif(Slicer3_FOUND)

install(TARGETS
  CreateMRML
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
