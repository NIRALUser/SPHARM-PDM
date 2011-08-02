project(shapeAnalysisMANCOVA_Wizard)

cmake_minimum_required(VERSION 2.6)

find_package(Slicer3 QUIET NO_DEFAULT_PATH)
if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

 find_package(BatchMake REQUIRED)
 include(${BatchMake_USE_FILE})

else(Slicer3_FOUND)
  find_package(ITK REQUIRED)
  if(ITK_FOUND)
    include(${USE_ITK_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  endif(ITK_FOUND)

  find_package(GenerateCLP REQUIRED)
  if(GenerateCLP_FOUND)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_FOUND)
    message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  endif(GenerateCLP_FOUND)

  find_package(BatchMake REQUIRED)
  if(BatchMake_FOUND)
    include(${BatchMake_USE_FILE})
  endif(BatchMake_FOUND)
endif(Slicer3_FOUND)

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


if(Slicer3_FOUND)
 slicer3_set_plugins_output_path(shapeAnalysisMANCOVA_Wizard)
 slicer3_install_plugins(shapeAnalysisMANCOVA_Wizard)
endif(Slicer3_FOUND)

install(TARGETS
  shapeAnalysisMANCOVA_Wizard
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

