cmake_minimum_required(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

# Testing for 

macro( BULD_AND_INSTALL_TESTS NAME TEST_FILE)
  add_executable(${NAME}Test ${TEST_FILE})
  target_link_libraries(${NAME}Test ${NAME}Lib)
  if( ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
    install( TARGETS ${NAME}Test DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} )
  endif()
  set(${NAME}TestPATH $<TARGET_FILE:${NAME}Test>)
endmacro()

if( NOT ${LOCAL_PROJECT_NAME}_SUPERBUILD )
  BULD_AND_INSTALL_TESTS( SegPostProcessCLP itkTest.cxx)
  BULD_AND_INSTALL_TESTS( GenParaMeshCLP vtkTest.cxx)
  BULD_AND_INSTALL_TESTS( ParaToSPHARMMeshCLP vtkTest.cxx)
endif()



set(TestDataDirectory ${CMAKE_CURRENT_BINARY_DIR}/Data)
set(TestingBinarydirectory ${CMAKE_CURRENT_BINARY_DIR}/Temporary)

## Test 1 : Segpostprocess
add_test(NAME SegPostProcessCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

## Test 2
set(InputImageForSeg ${TestDataDirectory}/InputImageForSegPostProcessCLP.nrrd)
set(ReferenceImageSeg ${TestDataDirectory}/OutputImageToCompareSegPostProcessCLP.nrrd)
set(OutputImageSeg ${TestingBinarydirectory}/OutputImageSegPostProcessCLP.nrrd)

#       Downloading of data on Girder
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f4568d777f16d01e0634/download" ${InputImageForSeg})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f4578d777f16d01e0637/download" ${ReferenceImageSeg})

add_test(
    NAME SegPostProcessCLPTest
    COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLPTest>
  --compare
    ${ReferenceImageSeg}
    ${OutputImageSeg}
  ModuleEntryPoint
    ${InputImageForSeg}
    ${OutputImageSeg}
    --UseUnitSpacing
  )

## Test 3 : GenParaMeshCLP
add_test(NAME GenParaMeshCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:GenParaMeshCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

## Test 4

set(InputImageForGen ${TestDataDirectory}/InputImageForGenParaMeshCLP.nrrd)
set(ReferenceImageGen ${TestDataDirectory}/OutputImageToCompareGenParaMeshCLP_surf.vtk)
set(OutputImagePara ${TestingBinarydirectory}/OutputImageGenParaMeshCLP_para.vtk)
set(OutputImageSurf ${TestingBinarydirectory}/OutputImageGenParaMeshCLP_surf.vtk)

#       Downloading of data on Girder
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f0ad8d777f16d01e0616/download" ${InputImageForGen})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f0b48d777f16d01e0619/download" ${ReferenceImageGen})

 add_test(
    NAME GenParaMeshCLPTest
    COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:GenParaMeshCLPTest>
  --compare
    ${ReferenceImageGen}
    ${OutputImageSurf}
  ModuleEntryPoint
    ${InputImageForGen}
    ${OutputImagePara}
    ${OutputImageSurf}
    --iter 100
  )

## Test 5 : ParaToSPHARMMeshCLP
add_test(NAME ParaToSPHARMMeshCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:ParaToSPHARMMeshCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

## Test 6

set(InputImageForSPHARM_Para ${TestDataDirectory}/InputImageTestForSPHARM_para.vtk)
set(InputImageForSPHARM_Surf ${TestDataDirectory}/InputImageTestForSPHARM_surf.vtk)
set(RegTemplateForSPHARM ${TestDataDirectory}/RegistrationTemplateForParaToSPHARMMeshCLP.vtk)

set(ReferenceImageSPHARM ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMeshCLP_SPHARM.vtk)
set(ReferenceImageSPHARM_MedialMesh ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMeshCLP_MedialMesh.vtk)
set(ReferenceImageSPHARM_Ellalign ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMeshCLP_SPHARM_ellalign.vtk)
set(ReferenceImageSPHARM_Procalign ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMeshCLP_SPHARM_procalign.vtk)

set(OutputDirectorySPHARM ${TestingBinarydirectory}/OutputImageParaToSPHARMMeshCLP)
set(OutputImageSPHARM ${TestingBinarydirectory}/OutputImageParaToSPHARMMeshCLP_SPHARM.vtk)
set(OutputImageSPHARM_medialMesh ${TestingBinarydirectory}/OutputImageParaToSPHARMMeshCLP_SPHARMMedialMesh.vtk)
set(OutputImageSPHARM_Ellalign ${TestingBinarydirectory}/OutputImageParaToSPHARMMeshCLP_SPHARM_ellalign.vtk)
set(OutputImageSPHARM_Procalign ${TestingBinarydirectory}/OutputImageParaToSPHARMMeshCLP_SPHARM_procalign.vtk)

#       Downloading of data on Girder
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f3038d777f16d01e061f/download" ${InputImageForSPHARM_Para})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f30a8d777f16d01e0622/download" ${InputImageForSPHARM_Surf})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f39e8d777f16d01e0631/download" ${RegTemplateForSPHARM})

file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f3188d777f16d01e0628/download" ${ReferenceImageSPHARM})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f3188d777f16d01e0625/download" ${ReferenceImageSPHARM_MedialMesh})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f3198d777f16d01e062b/download" ${ReferenceImageSPHARM_Ellalign})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/5909f31a8d777f16d01e062e/download" ${ReferenceImageSPHARM_Procalign})

add_test(
    NAME ParaToSPHARMMeshCLPTest
    COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:ParaToSPHARMMeshCLPTest>
  --compare
    ${ReferenceImageSPHARM}
    ${OutputImageSPHARM}
  --compare
    ${ReferenceImageSPHARM_MedialMesh}
    ${OutputImageSPHARM_medialMesh}
  --compare
    ${ReferenceImageSPHARM_Ellalign}
    ${OutputImageSPHARM_Ellalign}
  --compare
    ${ReferenceImageSPHARM_Procalign}
    ${OutputImageSPHARM_Procalign}
  ModuleEntryPoint
    ${InputImageForSPHARM_Para}
    ${InputImageForSPHARM_Surf}
    ${OutputDirectorySPHARM}
    --subdivLevel 10
    --spharmDegree 15
    --medialMesh
    --FinalFlip 3
    --regTemplateFileOn
    --regTemplate ${RegTemplateForSPHARM}
 )
