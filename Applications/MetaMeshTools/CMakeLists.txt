PROJECT(MetaMeshTools)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(Slicer3 QUIET NO_DEFAULT_PATH)

IF (Slicer3_FOUND)
  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

ELSE (Slicer3_FOUND)
  FIND_PACKAGE(ITK REQUIRED)
  IF (ITK_FOUND)
    INCLUDE(${USE_ITK_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  ENDIF (ITK_FOUND)

  FIND_PACKAGE(VTK REQUIRED)
  IF (VTK_FOUND)
    INCLUDE(${USE_VTK_FILE})
  ELSE(VTK_FOUND)
    MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
  ENDIF (VTK_FOUND)

ENDIF (Slicer3_FOUND)






# Shorthand to source and binary directories
SET(MY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MY_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Include directories
INCLUDE_DIRECTORIES( ${MY_SOURCE_DIR} ${MY_BINARY_DIR})

LINK_DIRECTORIES( ${CMAKE_BINARY_SOURCE_DIR})
LINK_LIBRARIES( 
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkHybrid
  vtkFiltering 
)

SET(SRC_STL2META_CXX STL2Meta.cxx vtkPolyDataToitkMesh.cxx)
ADD_EXECUTABLE(STL2Meta ${SRC_STL2META_CXX})
TARGET_LINK_LIBRARIES(STL2Meta ${ITK_LIBRARIES})


SET(SRC_META2STL_CXX Meta2STL.cxx itkMeshTovtkPolyData.cxx)
ADD_EXECUTABLE(Meta2STL ${SRC_META2STL_CXX})
TARGET_LINK_LIBRARIES(Meta2STL ${ITK_LIBRARIES})

SET(SRC_MESHMATH_CXX MeshMath.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
ADD_EXECUTABLE(MeshMath ${SRC_MESHMATH_CXX})
TARGET_LINK_LIBRARIES(MeshMath ${ITK_LIBRARIES})


ADD_EXECUTABLE(Meta2VTK Meta2VTK.cxx itkMeshTovtkPolyData.cxx)
TARGET_LINK_LIBRARIES(Meta2VTK ${ITK_LIBRARIES})

ADD_EXECUTABLE(VTK2Meta VTK2Meta.cxx vtkPolyDataToitkMesh.cxx)
TARGET_LINK_LIBRARIES(VTK2Meta ${ITK_LIBRARIES} )

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(STL2Meta)
  slicer3_install_plugins(STL2Meta)
 slicer3_set_plugins_output_path(Meta2STL)
  slicer3_install_plugins(Meta2STL)
 slicer3_set_plugins_output_path(MeshMath)
  slicer3_install_plugins(MeshMath)
 slicer3_set_plugins_output_path(Meta2VTK)
  slicer3_install_plugins(Meta2VTK)
 slicer3_set_plugins_output_path(VTK2Meta)
  slicer3_install_plugins(VTK2Meta)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  STL2Meta
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
INSTALL(TARGETS
  Meta2STL
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
INSTALL(TARGETS
  MeshMath
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
INSTALL(TARGETS
 Meta2VTK
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
INSTALL(TARGETS
  VTK2Meta
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

