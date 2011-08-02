#link_libraries(ShapeAlgorithms)

option(COMPILE_ImageMath "Compile ImageMath." ON)
if(COMPILE_ImageMath)
  subdirs(ImageMath)
endif(COMPILE_ImageMath)

option(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
if(COMPILE_MetaMeshTools)
  subdirs(MetaMeshTools)
endif(COMPILE_MetaMeshTools)

option(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
if(COMPILE_SegPostProcessCLP)
  subdirs(SegPostProcessCLP)
endif(COMPILE_SegPostProcessCLP)

option(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
if(COMPILE_GenParaMeshCLP)
  subdirs(GenParaMeshCLP)
endif(COMPILE_GenParaMeshCLP)

option(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
if(COMPILE_ParaToSPHARMMeshCLP)
  subdirs(ParaToSPHARMMeshCLP)
endif(COMPILE_ParaToSPHARMMeshCLP)

option(COMPILE_ShapeAnalysisModule "Compile ShapeAnalysisModule." ON)
if(COMPILE_ShapeAnalysisModule)
  subdirs(ShapeAnalysisModule)
endif(COMPILE_ShapeAnalysisModule)

option(COMPILE_StatNonParamTestPDM "Compile StatNonParamTestPDM." ON)
if(COMPILE_StatNonParamTestPDM)
  subdirs(StatNonParamTestPDM)
endif(COMPILE_StatNonParamTestPDM)

option(COMPILE_ParticleModule "Compile ParticleModule." ON)
if(COMPILE_ParticleModule)
  subdirs(ParticleModule)
endif(COMPILE_ParticleModule)

option(COMPILE_shapeAnalysisMANCOVA_Wizard "Compile shapeAnalysisMANCOVA_Wizard." ON)
if(COMPILE_shapeAnalysisMANCOVA_Wizard)
  subdirs(shapeAnalysisMANCOVA_Wizard)
endif(COMPILE_shapeAnalysisMANCOVA_Wizard)


option(COMPILE_CreateMRML "Compile CreateMRML." OFF)

if(Slicer3_FOUND)

#option(COMPILE_CreateMRML "Compile CreateMRML." ON)



endif(Slicer3_FOUND)
if(COMPILE_CreateMRML)
  subdirs(CreateMRML)
endif(COMPILE_CreateMRML)
