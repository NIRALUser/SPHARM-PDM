cmake_minimum_required(VERSION 2.8)
CMAKE_POLICY(VERSION 2.8)

# Testing for 

set(TestingBinarydirectory ${CMAKE_CURRENT_BINARY_DIR}/Temporary)
set(TestDataDirectory ${CMAKE_CURRENT_SOURCE_DIR}/Data)


## Test 1 : Segpostprocess
#add_test(NAME SegPostProcessCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

add_executable(SegPostProcessCLPTest itkTest.cxx)
target_link_libraries(SegPostProcessCLPTest SegPostProcessCLPLib)

set(InputImage ${TestDataDirectory}/InputImageTest.nrrd )
set(OutputImageSeg ${TestingBinarydirectory}/OutputImageSegPostProcessCLP.nrrd )
set(ReferenceImage ${TestDataDirectory}/OutputImageTestToCompare.nrrd)

## Test 2 
add_test(NAME SegPostProcessCLPTest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLPTest>
  --compare
    ${ReferenceImage}
    ${OutputImageSeg}
  ModuleEntryPoint
    ${InputImage}
    ${OutputImageSeg}
    --UseUnitSpacing
  )

## Test 3 : GenParaMeshCLP
add_test(NAME GenParaMeshCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:GenParaMeshCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

## Test 4 
add_executable(GenParaMeshCLPTest vtkTest.cxx)
target_link_libraries(GenParaMeshCLPTest GenParaMeshCLPLib)

set(ReferenceImageGen ${TestDataDirectory}/OutputImageToCompareGenParaMeshCLP_surf.vtk)
set(OutputImagePara ${TestingBinarydirectory}/OutputImageGenParaMeshCLP_para.vtk)
set(OutputImageSurf ${TestingBinarydirectory}/OutputImageGenParaMeshCLP_surf.vtk)
set(InputImageForGen ${TestDataDirectory}/InputImageTestForGen.nrrd)

if( NOT VTK_DEBUG_LEAKS )
add_test(NAME GenParaMeshCLPTest COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:GenParaMeshCLPTest>
  --compare
    ${ReferenceImageGen}
    ${OutputImageSurf}
  ModuleEntryPoint
    ${InputImageForGen}
    ${OutputImagePara}
    ${OutputImageSurf}
    --iter 100
  )
endif()
