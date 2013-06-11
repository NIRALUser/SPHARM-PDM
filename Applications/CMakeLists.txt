#link_libraries(ShapeAlgorithms)
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


set(external_projects_vars
  -DMAKECOMMAND:STRING=${MAKECOMMAND}
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
  -DCMAKE_CXX_COMPILER:PATH=${CMAKE_CXX_COMPILER}
  -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_COMPILER:PATH=${CMAKE_C_COMPILER}
  -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
  -DCMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
  -DCMAKE_EXTRA_GENERATOR:STRING=${CMAKE_EXTRA_GENERATOR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  -DCMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH=${CMAKE_BUNDLE_OUTPUT_DIRECTORY}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
  )


option(COMPILE_ImageMath "Compile ImageMath." ON)
if(COMPILE_ImageMath)
  set(ImageMath_REPOSITORY https://www.nitrc.org/svn/niral_utilities/trunk/ImageMath CACHE STRING "" FORCE)
  ExternalProject_Add(ImageMath
    SVN_REPOSITORY ${ImageMath_REPOSITORY}
    SVN_REVISION -r 36
    SVN_USERNAME slicerbot
    SVN_PASSWORD slicer
    UPDATE_COMMAND ""
    SOURCE_DIR ImageMath
    BINARY_DIR ImageMath-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DITK_DIR:PATH=${ITK_DIR}
      ${external_projects_vars}
      INSTALL_COMMAND ""
  )
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

option(COMPILE_shapeAnalysisMANCOVA_Wizard "Compile shapeAnalysisMANCOVA_Wizard." OFF)
option(COMPILE_StatNonParamTestPDM "Compile StatNonParamTestPDM." ON)
if(COMPILE_StatNonParamTestPDM)
  ###add download Mancova and configure mancova project. Qt for wizard....
  set(MANCOVA_REPOSITORY https://www.nitrc.org/svn/shape_mancova CACHE STRING "" FORCE)
  ExternalProject_Add(MANCOVA
    SVN_REPOSITORY ${MANCOVA_REPOSITORY}
    SVN_REVISION -r 4
    SVN_USERNAME slicerbot
    SVN_PASSWORD slicer
    UPDATE_COMMAND ""
    SOURCE_DIR MANCOVA
    BINARY_DIR MANCOVA-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DCOMPILE_shapeAnalysisMANCOVA_Wizard:BOOL={COMPILE_shapeAnalysisMANCOVA_Wizard}
      -DITK_DIR:PATH=${ITK_DIR}
      -DVTK_DIR:PATH=${VTK_DIR}
      -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
      -DBoost_DIR:PATH=${Boost_DIR}/share/boost-${boost_version}/cmake
      ${external_projects_vars}
      INSTALL_COMMAND ""
  )
endif(COMPILE_StatNonParamTestPDM)

option(COMPILE_shapeworks "Compile shapeworks." ON)
if(COMPILE_shapeworks)
  set(shapeworks_REPOSITORY  https://www.nitrc.org/svn/shapeworks/code CACHE STRING "" FORCE)
  ExternalProject_Add(shapeworks
    SVN_REPOSITORY ${shapeworks_REPOSITORY}
    SVN_REVISION -r 78
    SVN_USERNAME slicerbot
    SVN_PASSWORD slicer
    UPDATE_COMMAND ""
    SOURCE_DIR shapeworks
    BINARY_DIR shapeworks-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      -DITK_DIR:PATH=${ITK_DIR}
      ${external_projects_vars}
      INSTALL_COMMAND ""
  )
endif(COMPILE_shapeworks)

option(COMPILE_ParticleModule "Compile ParticleModule." ON)
if(COMPILE_ParticleModule)
  add_subdirectory(ParticleModule)
endif(COMPILE_ParticleModule)



