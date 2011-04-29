cmake_minimum_required(VERSION 2.6)
project (CreateMRML)

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

  find_package(Slicer3 REQUIRED)
  include(${Slicer3_USE_FILE})

add_executable(CreateMRML ${CreateMRML_SOURCE})
target_link_libraries(CreateMRML ${Slicer3_Libs_LIBRARIES} )


