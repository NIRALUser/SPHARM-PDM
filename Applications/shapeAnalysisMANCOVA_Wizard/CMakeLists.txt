#project(shapeAnalysisMANCOVA_Wizard)

find_package(Qt4 REQUIRED)
if(QT_USE_FILE)
  include_directories(${QT_INCLUDE_DIR})
  include(${QT_USE_FILE})
else(QT_USE_FILE)
   message(FATAL_ERROR, "QT not found. Please set QT_DIR.")
endif(QT_USE_FILE)

set(shapeAnalysisMANCOVA_Wizard_SRCS shapeAnalysisMANCOVA_Wizard.cxx main.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx )
GENERATECLP(main.cxx shapeAnalysisMANCOVA_Wizard.xml)

set(GUI_HDRS shapeAnalysisMANCOVA_Wizard.ui)
QT4_WRAP_UI(GUI_files ${GUI_HDRS})
set(MANCOVA_HDRS shapeAnalysisMANCOVA_Wizard.h)
QT4_WRAP_CPP(MOC_files ${MANCOVA_HDRS})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(CXX_FILES
  shapeAnalysisMANCOVA_Wizard.cxx
  main.cxx
)

add_definitions(-DQT_CORE_LIB -DQT3_SUPPORT)
add_executable(shapeAnalysisMANCOVA_Wizard ${CXX_FILES} ${MOC_files} ${GUI_files})

target_link_libraries(shapeAnalysisMANCOVA_Wizard
  ${ITK_LIBRARIES}
  ${QT_LIBRARIES}
)


install(TARGETS
  shapeAnalysisMANCOVA_Wizard
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

