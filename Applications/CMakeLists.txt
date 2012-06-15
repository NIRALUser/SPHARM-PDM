#link_libraries(ShapeAlgorithms)

option(COMPILE_ImageMath "Compile ImageMath." ON)
if(COMPILE_ImageMath)
  add_subdirectory(ImageMath)
endif(COMPILE_ImageMath)

option(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
if(COMPILE_MetaMeshTools)
  add_subdirectory(MetaMeshTools)
endif(COMPILE_MetaMeshTools)

option(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
if(COMPILE_SegPostProcessCLP)
  add_subdirectory(SegPostProcessCLP)
endif(COMPILE_SegPostProcessCLP)

option(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
if(COMPILE_GenParaMeshCLP)
  add_subdirectory(GenParaMeshCLP)
endif(COMPILE_GenParaMeshCLP)

option(COMPILE_RadiusToMesh "Compile GenParaMeshCLP." ON)
if(COMPILE_RadiusToMesh)
  add_subdirectory(RadiusToMesh)
endif(COMPILE_RadiusToMesh)

option(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
if(COMPILE_ParaToSPHARMMeshCLP)
  add_subdirectory(ParaToSPHARMMeshCLP)
endif(COMPILE_ParaToSPHARMMeshCLP)

option(COMPILE_ShapeAnalysisModule "Compile ShapeAnalysisModule." ON)
if(COMPILE_ShapeAnalysisModule)
  add_subdirectory(ShapeAnalysisModule)
endif(COMPILE_ShapeAnalysisModule)

option(COMPILE_StatNonParamTestPDM "Compile StatNonParamTestPDM." ON)
if(COMPILE_StatNonParamTestPDM)
  add_subdirectory(StatNonParamTestPDM)
endif(COMPILE_StatNonParamTestPDM)

option(COMPILE_ParticleModule "Compile ParticleModule." ON)
if(COMPILE_ParticleModule)
  add_subdirectory(ParticleModule)
endif(COMPILE_ParticleModule)

if(Qt_FOUND)
option(COMPILE_shapeAnalysisMANCOVA_Wizard "Compile shapeAnalysisMANCOVA_Wizard." ON)
if(COMPILE_shapeAnalysisMANCOVA_Wizard)
  add_subdirectory(shapeAnalysisMANCOVA_Wizard)
endif(COMPILE_shapeAnalysisMANCOVA_Wizard)
endif()


option(COMPILE_CreateMRML "Compile CreateMRML." OFF)
if(COMPILE_CreateMRML)
  add_subdirectory(CreateMRML)
endif(COMPILE_CreateMRML)
