PROJECT(shapeAnalysisMANCOVA_Wizard)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(Slicer3 QUIET NO_DEFAULT_PATH)
IF (Slicer3_FOUND)
  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

 FIND_PACKAGE(BatchMake REQUIRED)
 INCLUDE(${BatchMake_USE_FILE})

ELSE (Slicer3_FOUND)
  FIND_PACKAGE(ITK REQUIRED)
  IF (ITK_FOUND)
    INCLUDE(${USE_ITK_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  ENDIF (ITK_FOUND)

  FIND_PACKAGE(GenerateCLP REQUIRED)
  IF(GenerateCLP_FOUND)
    INCLUDE(${GenerateCLP_USE_FILE})
  ELSE(GenerateCLP_FOUND)
    MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  ENDIF(GenerateCLP_FOUND)

  FIND_PACKAGE(BatchMake REQUIRED)
  IF (BatchMake_FOUND)
    INCLUDE(${BatchMake_USE_FILE})
  ENDIF (BatchMake_FOUND)
ENDIF (Slicer3_FOUND)

FIND_PACKAGE(Qt4 REQUIRED)
IF(QT_USE_FILE)
  INCLUDE_DIRECTORIES(${QT_INCLUDE_DIR})
  INCLUDE(${QT_USE_FILE})
ELSE(QT_USE_FILE)
   MESSAGE(FATAL_ERROR, "QT not found. Please set QT_DIR.")
ENDIF(QT_USE_FILE)

SET (shapeAnalysisMANCOVA_Wizard_SRCS shapeAnalysisMANCOVA_Wizard.cxx main.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx )
GENERATECLP(main.cxx shapeAnalysisMANCOVA_Wizard.xml)

SET (GUI_HDRS shapeAnalysisMANCOVA_Wizard.ui)
QT4_WRAP_UI(GUI_files ${GUI_HDRS})
SET (MANCOVA_HDRS shapeAnalysisMANCOVA_Wizard.h)
QT4_WRAP_CPP(MOC_files ${MANCOVA_HDRS})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(CXX_FILES 
  shapeAnalysisMANCOVA_Wizard.cxx
  main.cxx
)

ADD_DEFINITIONS(-DQT_CORE_LIB -DQT3_SUPPORT)
ADD_EXECUTABLE(shapeAnalysisMANCOVA_Wizard ${CXX_FILES} ${MOC_files} ${GUI_files})

TARGET_LINK_LIBRARIES(shapeAnalysisMANCOVA_Wizard
  ${ITK_LIBRARIES}
  ${QT_LIBRARIES}
)


IF (Slicer3_FOUND)
 slicer3_set_plugins_output_path(shapeAnalysisMANCOVA_Wizard)
 slicer3_install_plugins(shapeAnalysisMANCOVA_Wizard)
ENDIF (Slicer3_FOUND)

INSTALL(TARGETS
  shapeAnalysisMANCOVA_Wizard
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

