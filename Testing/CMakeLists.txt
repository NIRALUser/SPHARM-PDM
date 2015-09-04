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

if( ${LOCAL_PROJECT_NAME}_SUPERBUILD )
  set(GenParaMeshCLPTestPATH ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}/GenParaMeshCLPTest)
  set(SegPostProcessCLPTestPATH ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}/SegPostProcessCLPTest)
  set(GenParaMeshCLPPATH ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}/GenParaMeshCLP)
else()
  BULD_AND_INSTALL_TESTS( SegPostProcessCLP itkTest.cxx)
  BULD_AND_INSTALL_TESTS( GenParaMeshCLP vtkTest.cxx)
  set(GenParaMeshCLPPATH $<TARGET_FILE:GenParaMeshCLP>)
endif()



set(TestingBinarydirectory ${CMAKE_CURRENT_BINARY_DIR}/Temporary)
set(TestDataDirectory ${CMAKE_CURRENT_SOURCE_DIR}/Data)


## Test 1 : Segpostprocess
#add_test(NAME SegPostProcessCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} $<TARGET_FILE:SegPostProcessCLP> --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

set(InputImage ${TestDataDirectory}/InputImageTest.nrrd )
set(OutputImageSeg ${TestingBinarydirectory}/OutputImageSegPostProcessCLP.nrrd )
set(ReferenceImage ${TestDataDirectory}/OutputImageTestToCompare.nrrd)

## Test 2 
add_test(NAME SegPostProcessCLPTest COMMAND ${Slicer_LAUNCH_COMMAND} ${SegPostProcessCLPTestPATH}
  --compare
    ${ReferenceImage}
    ${OutputImageSeg}
  ModuleEntryPoint
    ${InputImage}
    ${OutputImageSeg}
    --UseUnitSpacing
  )

## Test 3 : GenParaMeshCLP
add_test(NAME GenParaMeshCLPHelp COMMAND ${Slicer_LAUNCH_COMMAND} ${GenParaMeshCLPPATH} --help) # Should just display the Usage : ${CMAKE_CURRENT_BINARY_DIR} = DTIAtlasBuilder-build/Testing -> need to ../

## Test 4

set(ReferenceImageGen ${TestDataDirectory}/OutputImageToCompareGenParaMeshCLP_surf.vtk)
set(OutputImagePara ${TestingBinarydirectory}/OutputImageGenParaMeshCLP_para.vtk)
set(OutputImageSurf ${TestingBinarydirectory}/OutputImageGenParaMeshCLP_surf.vtk)
set(InputImageForGen ${TestDataDirectory}/InputImageTestForGen.nrrd)

if( NOT VTK_DEBUG_LEAKS )
  add_test(NAME GenParaMeshCLPTest COMMAND ${Slicer_LAUNCH_COMMAND} ${GenParaMeshCLPTestPATH}
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
