
SET(StatNonParamTestPDMSRCS 
  createPerm.cxx 
  main.cxx
  computeStatPDM.cxx 
  computeStatPDM_vector.cxx 
  itkMeshTovtkPolyData.cxx 
  vtkPolyDataToitkMesh.cxx
)

ADD_EXECUTABLE( StatNonParamTestPDM ${StatNonParamTestPDMSRCS} )

TARGET_LINK_LIBRARIES ( StatNonParamTestPDM 
  ${ITK_LIBRARIES}  
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkCommon
  vtkHybrid
  vtkFiltering 
)

FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp)

INCLUDE_DIRECTORIES(
${BOOST_INCLUDE_DIR}
)

FIND_PATH(BOOST_LIBRARY_DIR boost_date_time-vc71-mt.lib)
LINK_DIRECTORIES(${BOOST_LIBRARY_DIR})

# Win32 has automatic linking of boost libraries so only add boost
# libaries on unix
IF(NOT WIN32 OR CYGWIN)
  SET(BOOST_LIBRARIES boost_program_options boost_system boost_filesystem)
ENDIF(NOT WIN32 OR CYGWIN)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF(GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
  SET (shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.cxx MANCOVA.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx  miscMANCOVA.cxx nonMANCOVA.cxx)
  GENERATECLP(shapeAnalysisMANCOVA.cxx shapeAnalysisMANCOVA.xml)

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
ENDIF(GenerateCLP_FOUND)


