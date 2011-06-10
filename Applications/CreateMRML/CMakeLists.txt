PROJECT(CreateMRML)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(Slicer3 REQUIRED)
IF (Slicer3_FOUND)
  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
ENDIF(Slicer3_FOUND)

SET(CreateMRML_SOURCE
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

ADD_EXECUTABLE(CreateMRML ${CreateMRML_SOURCE})
TARGET_LINK_LIBRARIES(CreateMRML ${Slicer3_Libs_LIBRARIES} )

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(CreateMRML)
  slicer3_install_plugins(CreateMRML)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  CreateMRML
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
