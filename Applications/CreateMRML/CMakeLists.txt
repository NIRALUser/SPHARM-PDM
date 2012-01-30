#project(CreateMRML)

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

install(TARGETS
  CreateMRML
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
