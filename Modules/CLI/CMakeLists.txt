option(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
if(COMPILE_MetaMeshTools)
  add_subdirectory(MetaMeshTools)
endif()

option(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
if(COMPILE_SegPostProcessCLP)
  add_subdirectory(SegPostProcessCLP)
endif()

option(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
if(COMPILE_GenParaMeshCLP)
  add_subdirectory(GenParaMeshCLP)
endif()

option(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
if(COMPILE_ParaToSPHARMMeshCLP)
  add_subdirectory(ParaToSPHARMMeshCLP)
endif()

option(COMPILE_SpharmTool "Compile SpharmTool." ON)
if(COMPILE_SpharmTool)
  add_subdirectory(SpharmTool)
endif()
