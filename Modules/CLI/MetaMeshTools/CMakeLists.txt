
set(SRC_STL2META_CXX STL2Meta.cxx vtkPolyDataToitkMesh.cxx)
add_executable(STL2Meta ${SRC_STL2META_CXX})
target_link_libraries(STL2Meta PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES})

set(SRC_META2STL_CXX Meta2STL.cxx itkMeshTovtkPolyData.cxx)
add_executable(Meta2STL ${SRC_META2STL_CXX})
target_link_libraries(Meta2STL PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES})

set(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
add_executable(MeshMath ${SRC_MESHMATH_CXX})
target_link_libraries(MeshMath PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES} ShapeAlgorithms)

add_executable(Meta2VTK Meta2VTK.cxx itkMeshTovtkPolyData.cxx)
target_link_libraries(Meta2VTK PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES})

add_executable(VTK2Meta VTK2Meta.cxx vtkPolyDataToitkMesh.cxx)
target_link_libraries(VTK2Meta PRIVATE ${ITK_LIBRARIES} ${VTK_LIBRARIES} )

add_executable(asc2vtk asc2vtk.cxx)
add_executable(asc2meta asc2meta.cxx)

add_executable(BYU2VTK BYU2VTK.cxx)
target_link_libraries(BYU2VTK PRIVATE ${VTK_LIBRARIES} )

add_executable(STL2VTK STL2VTK.cxx)
target_link_libraries(STL2VTK PRIVATE ${VTK_LIBRARIES} )

macro(INSTALL_EXECUTABLE)
  set(options )
  set( oneValueArgs OUTPUT_DIR )
  set(multiValueArgs LIST_EXEC )
  cmake_parse_arguments(LOCAL
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  foreach( tool ${LOCAL_LIST_EXEC})
    install(TARGETS ${tool} DESTINATION ${LOCAL_OUTPUT_DIR} )
  endforeach()
endmacro()

if( NOT SPHARM-PDM_BUILD_SLICER_EXTENSION )
  install(TARGETS
   STL2Meta
   Meta2STL
   MeshMath
   Meta2VTK
   BYU2VTK
   VTK2Meta
   STL2VTK
   asc2vtk
   asc2meta
   BUNDLE DESTINATION  bin
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib/static  
   )
else()
  set( NotCLIToolsList
    MeshMath STL2Meta Meta2STL Meta2VTK VTK2Meta asc2vtk asc2meta BYU2VTK STL2VTK
     )
  INSTALL_EXECUTABLE( OUTPUT_DIR ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} LIST_EXEC ${NotCLIToolsList} )
endif()
