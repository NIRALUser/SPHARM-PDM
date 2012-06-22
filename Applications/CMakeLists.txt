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


find_package(Slicer3 QUIET )
if( Slicer3_FOUND )
  option(COMPILE_CreateMRML "Compile CreateMRML." ON)
  if(COMPILE_CreateMRML)
    include(ExternalProject)
    # Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()
    set( GIT_PROTOCOL "git" )
    if(DEFINED ${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL)
      if( NOT ${${CMAKE_PROJECT_NAME}_USE_GIT_PROTOCOL} )
        set(GIT_PROTOCOL "http")
      endif()
    endif()
    set(CreateMRML_REPOSITORY ${GIT_PROTOCOL}://github.com/booboo69500/CreateMRML.git CACHE STRING "" FORCE)
    ExternalProject_Add(CreateMRML
      GIT_REPOSITORY ${CreateMRML_REPOSITORY}
      GIT_TAG "1e17d6dc66dca98b5622918408825d888ebd5a42"
      UPDATE_COMMAND ""
      SOURCE_DIR CreateMRML
      BINARY_DIR CreateMRML-build
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
        -DSlicer3_DIR:PATH=${Slicer3_DIR}
        INSTALL_COMMAND ""
        #    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}
    )
  endif()
endif()
