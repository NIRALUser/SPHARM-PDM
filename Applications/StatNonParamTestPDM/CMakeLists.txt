PROJECT(shapeAnalysisMANCOVA)

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

  FIND_PACKAGE(GenerateCLP REQUIRED)
  IF(GenerateCLP_FOUND)
    INCLUDE(${GenerateCLP_USE_FILE})
  ELSE(GenerateCLP_FOUND)
    MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  ENDIF(GenerateCLP_FOUND)

ENDIF (Slicer3_FOUND)

IF(NOT BOOST_FOUND)
  FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp
  /opt/local/Boost/boost_1_39_0
  /usr/include
  /usr/local/include
  )
  INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})

#  FIND_PATH(BOOST_LIBRARY_DIR boost_date_time-vc71-mt.lib)
  FIND_PATH(BOOST_LIBRARY_DIR libboost_filesystem.a
    /opt/local/Boost/boost_1_39_0_linux64/lib
    /usr/lib64 /usr/lib /opt/local/lib /usr/local/lib /tools/lib  
    DOC " directory for boost, where libboost can be found"
  )

  LINK_DIRECTORIES(${BOOST_LIBRARY_DIR})

  # Win32 has automatic linking of boost libraries so only add boost
  # libaries on unix
  IF(NOT WIN32 OR CYGWIN)
    SET(BOOST_LIBRARIES boost_program_options boost_system boost_filesystem)
  ENDIF(NOT WIN32 OR CYGWIN)
ENDIF(NOT BOOST_FOUND)


   SET (shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.cxx MANCOVA.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx  miscMANCOVA.cxx nonMANCOVA.cxx )
  GENERATECLP(shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.xml)

  ADD_EXECUTABLE(shapeAnalysisMANCOVA ${shapeAnalysisMANCOVA_SRCS})

  TARGET_LINK_LIBRARIES(shapeAnalysisMANCOVA
    ITKCommon 
    ITKAlgorithms  
    vtkGraphics
    vtkIO
    ITKIO
    vtkCommon
    vtkFiltering
    ${BOOST_LIBRARIES}	
    )

IF (Slicer3_FOUND)
  slicer3_set_plugins_output_path(shapeAnalysisMANCOVA)
  slicer3_install_plugins(shapeAnalysisMANCOVA)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  shapeAnalysisMANCOVA
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )





