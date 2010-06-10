

FIND_PACKAGE(GenerateCLP REQUIRED)
IF (GenerateCLP_FOUND)
	INCLUDE(${GenerateCLP_USE_FILE})
ENDIF (GenerateCLP_FOUND)

FIND_PACKAGE(BatchMake REQUIRED)
IF (BatchMake_FOUND)
	INCLUDE(${BatchMake_USE_FILE})
ENDIF (BatchMake_FOUND)

 
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








