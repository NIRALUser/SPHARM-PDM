PROJECT(ShapeAnalysisModule)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

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

FIND_PACKAGE(GenerateCLP REQUIRED)
IF (GenerateCLP_FOUND)
	INCLUDE(${GenerateCLP_USE_FILE})
ENDIF (GenerateCLP_FOUND)

FIND_PACKAGE(BatchMake REQUIRED)
IF (BatchMake_FOUND)
	INCLUDE(${BatchMake_USE_FILE})
ENDIF (BatchMake_FOUND)

#FIND_PACKAGE(VTK REQUIRED)
#IF(VTK_FOUND)
	#INCLUDE(${VTK_USE_FILE})
#ENDIF(VTK_FOUND)
 
SET (ShapeAnalysisModule_SOURCE 
	ShapeAnalysisModuleCLP.cxx
	Parameters.cxx 
	Parameters.h
	ShapeAnalysisModuleComputation.cxx
	ShapeAnalysisModuleComputation.h
	
) 

GENERATECLP(ShapeAnalysisModule_SOURCE ShapeAnalysisModule.xml)
ADD_EXECUTABLE (ShapeAnalysisModule ${ShapeAnalysisModule_SOURCE})
#INCLUDE_DIRECTORIES("/opt/local/ITK_3.X/Insight_3.16.0/ITK_3.16.0_linux64_dyn/Utilities/itksys ")
# 
TARGET_LINK_LIBRARIES(ShapeAnalysisModule
  itksys
  vtkIO
  ${ITK_LIBRARIES}
  ${BatchMake_LIBRARIES}

)








