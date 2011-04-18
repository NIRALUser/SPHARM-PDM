##project (CreateMRML)

#set(CreateMRML_SOURCE
 #CreateMRML.cxx
 #MRMLColorableHelper.h
# MRMLColorableHelper.cxx
 #MRMLFiducialHelper.cxx
# MRMLFiducialHelper.h
# MRMLNodeHelper.h
 #MRMLModelHelper.h
 #MRMLTransformHelper.h
# MRMLVolumeHelper.h
# CreateMRMLSceneHelper.h
 #CreateMRMLSceneHelper.cxx
# MagicParser.h
# MagicParser.cxx
#)

 # find_package(Slicer3 REQUIRED)
 # include(${Slicer3_USE_FILE})
#
#add_executable(CreateMRML ${CreateMRML_SOURCE})
#target_link_libraries(CreateMRML ${Slicer3_Libs_LIBRARIES} )

#PROJECT(Mrml creator)

cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(ITK REQUIRED)
IF (ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
	MESSAGE(FATAL_ERROR Please set ITK_DIR.)
ENDIF (ITK_FOUND)

FIND_PACKAGE(VTK REQUIRED)
IF (VTK_FOUND)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
	MESSAGE(FATAL_ERROR Please set VTK_DIR.)
ENDIF (VTK_FOUND)

include_directories(/biomed-resimg/home/budin/LocalWork/devel/tempMRML/MRML-linux64-static /biomed-resimg/home/budin/LocalWork/devel/tempMRML/MRML)
link_directories(/biomed-resimg/home/budin/LocalWork/devel/tempMRML/MRML-linux64-static)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF (GenerateCLP_FOUND)
   INCLUDE(${GenerateCLP_USE_FILE})
ENDIF(GenerateCLP_FOUND)

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
ADD_EXECUTABLE (CreateMRML_SOURCE ${CreateMRML_SOURCE})

