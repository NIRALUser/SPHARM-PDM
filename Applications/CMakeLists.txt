#LINK_LIBRARIES(ShapeAlgorithms)

OPTION(COMPILE_ImageMath "Compile ImageMath." ON)
IF(COMPILE_ImageMath)
  SUBDIRS(ImageMath)
ENDIF(COMPILE_ImageMath)

OPTION(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
IF(COMPILE_MetaMeshTools)
  SUBDIRS(MetaMeshTools)
ENDIF(COMPILE_MetaMeshTools)

OPTION(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
IF(COMPILE_SegPostProcessCLP)
  SUBDIRS(SegPostProcessCLP)
ENDIF(COMPILE_SegPostProcessCLP)

OPTION(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
IF(COMPILE_GenParaMeshCLP)
  SUBDIRS(GenParaMeshCLP)
ENDIF(COMPILE_GenParaMeshCLP)

OPTION(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
IF(COMPILE_ParaToSPHARMMeshCLP)
  SUBDIRS(ParaToSPHARMMeshCLP)
ENDIF(COMPILE_ParaToSPHARMMeshCLP)

OPTION(COMPILE_ShapeAnalysisModule "Compile ShapeAnalysisModule." ON)
IF(COMPILE_ShapeAnalysisModule)
  SUBDIRS(ShapeAnalysisModule)
ENDIF(COMPILE_ShapeAnalysisModule)

OPTION(COMPILE_StatNonParamTestPDM "Compile StatNonParamTestPDM." ON)
IF(COMPILE_StatNonParamTestPDM)
  SUBDIRS(StatNonParamTestPDM)
ENDIF(COMPILE_StatNonParamTestPDM)

OPTION(COMPILE_ParticleModule "Compile ParticleModule." ON)
IF(COMPILE_ParticleModule)
  SUBDIRS(ParticleModule)
ENDIF(COMPILE_ParticleModule)

OPTION(COMPILE_shapeAnalysisMANCOVA_Wizard "Compile shapeAnalysisMANCOVA_Wizard." ON)
IF(COMPILE_shapeAnalysisMANCOVA_Wizard)
  SUBDIRS(shapeAnalysisMANCOVA_Wizard)
ENDIF(COMPILE_shapeAnalysisMANCOVA_Wizard)


OPTION(COMPILE_CreateMRML "Compile CreateMRML." OFF)

IF (Slicer3_FOUND)

#OPTION(COMPILE_CreateMRML "Compile CreateMRML." ON)



ENDIF (Slicer3_FOUND)
IF(COMPILE_CreateMRML)
  SUBDIRS(CreateMRML)
ENDIF(COMPILE_CreateMRML)
