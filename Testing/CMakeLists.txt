function(spharm_add_test_executable NAME TEST_FILE)
  add_executable(${NAME}Test ${TEST_FILE})
  target_link_libraries(${NAME}Test ${NAME}Lib)
endfunction()

spharm_add_test_executable( SegPostProcessCLP itkTest.cxx)
spharm_add_test_executable( GenParaMeshCLP vtkTest.cxx)
spharm_add_test_executable( ParaToSPHARMMeshCLP vtkTest.cxx)

set(TestDataDirectory ${CMAKE_CURRENT_BINARY_DIR}/Data)
set(TestingBinarydirectory ${CMAKE_CURRENT_BINARY_DIR}/Temporary)

#----------------------------------------------------------------------------------------
## Test 1 : Segpostprocess
add_test(NAME SegPostProcessCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

#----------------------------------------------------------------------------------------
## Test 2
set(InputImageForSeg ${TestDataDirectory}/InputImageForSegPostProcess.gipl)
set(ReferenceImageSeg ${TestDataDirectory}/OutputImageToCompareSegPostProcess.nrrd)
set(OutputImageSeg ${TestingBinarydirectory}/OutputImageSegPostProcess.nrrd)

#       Downloading of data on Girder
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59945eb38d777f7d33e9c3c1/download" ${InputImageForSeg})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59945ee08d777f7d33e9c3d3/download" ${ReferenceImageSeg})

add_test(
    NAME SegPostProcessCLPTest
    COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLPTest>
  --compare
    ${ReferenceImageSeg}
    ${OutputImageSeg}
  ModuleEntryPoint
    ${InputImageForSeg}
    ${OutputImageSeg}
    --rescale 
    --space 0.5,0.5,0.5
  )

#----------------------------------------------------------------------------------------
## Test 3 : GenParaMeshCLP
add_test(NAME GenParaMeshCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:GenParaMeshCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

#----------------------------------------------------------------------------------------
## Test 4

set(InputImageForGen ${TestDataDirectory}/InputImageForGenParaMesh.nrrd)
set(ReferenceImageGen ${TestDataDirectory}/OutputImageToCompareGenParaMesh_surf.vtk)
set(OutputImagePara ${TestingBinarydirectory}/OutputImageGenParaMesh_para.vtk)
set(OutputImageSurf ${TestingBinarydirectory}/OutputImageGenParaMesh_surf.vtk)

#       Downloading of data on Girder
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59945ee08d777f7d33e9c3d6/download" ${InputImageForGen})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59945ece8d777f7d33e9c3c7/download" ${ReferenceImageGen})

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
    --iter 25
  )

#----------------------------------------------------------------------------------------
## Test 5 : ParaToSPHARMMeshCLP
add_test(NAME ParaToSPHARMMeshCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:ParaToSPHARMMeshCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

#----------------------------------------------------------------------------------------
## Test 6

set(InputImageForSPHARM_Para ${TestDataDirectory}/InputImageTestForSPHARM_para.vtk)
set(InputImageForSPHARM_Surf ${TestDataDirectory}/InputImageTestForSPHARM_surf.vtk)
set(RegTemplateForSPHARM ${TestDataDirectory}/RegistrationTemplateForParaToSPHARMMesh.vtk)

set(ReferenceImageSPHARM ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMesh_SPHARM.vtk)
set(ReferenceImageSPHARM_MedialMesh ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMesh_MedialMesh.vtk)
set(ReferenceImageSPHARM_Ellalign ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMesh_SPHARM_ellalign.vtk)
set(ReferenceImageSPHARM_Procalign ${TestDataDirectory}/OutputImageToCompareParaToSPHARMMesh_SPHARM_procalign.vtk)

set(OutputDirectorySPHARM ${TestingBinarydirectory}/OutputImageParaToSPHARMMesh)
set(OutputImageSPHARM ${TestingBinarydirectory}/OutputImageParaToSPHARMMesh_SPHARM.vtk)
set(OutputImageSPHARM_medialMesh ${TestingBinarydirectory}/OutputImageParaToSPHARMMesh_SPHARMMedialMesh.vtk)
set(OutputImageSPHARM_Ellalign ${TestingBinarydirectory}/OutputImageParaToSPHARMMesh_SPHARM_ellalign.vtk)
set(OutputImageSPHARM_Procalign ${TestingBinarydirectory}/OutputImageParaToSPHARMMesh_SPHARM_procalign.vtk)

#       Downloading of data on Girder
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59af08798d777f7d33e9cf8e/download" ${InputImageForSPHARM_Para})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59945ed08d777f7d33e9c3cd/download" ${InputImageForSPHARM_Surf})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/599462f78d777f7d33e9c3e6/download" ${RegTemplateForSPHARM})

file(DOWNLOAD "https://data.kitware.com/api/v1/file/59af09028d777f7d33e9cf9a/download" ${ReferenceImageSPHARM})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59af09018d777f7d33e9cf94/download" ${ReferenceImageSPHARM_MedialMesh})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59af09018d777f7d33e9cf91/download" ${ReferenceImageSPHARM_Ellalign})
file(DOWNLOAD "https://data.kitware.com/api/v1/file/59af09028d777f7d33e9cf97/download" ${ReferenceImageSPHARM_Procalign})

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
    --paraOut
    --regTemplateFileOn
    --regTemplate ${RegTemplateForSPHARM}
 )
