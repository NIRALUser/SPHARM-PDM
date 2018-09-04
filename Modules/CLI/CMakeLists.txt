if(COMPILE_MetaMeshTools)
  add_subdirectory(MetaMeshTools)
endif()

if(COMPILE_SegPostProcessCLP)
  add_subdirectory(SegPostProcessCLP)
endif()

if(COMPILE_GenParaMeshCLP)
  add_subdirectory(GenParaMeshCLP)
endif()

if(COMPILE_ParaToSPHARMMeshCLP)
  add_subdirectory(ParaToSPHARMMeshCLP)
endif()

if(COMPILE_SpharmTool)
  add_subdirectory(SpharmTool)
endif()
