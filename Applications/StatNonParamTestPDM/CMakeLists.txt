project(shapeAnalysisMANCOVA)

cmake_minimum_required(VERSION 2.6)

find_package(Slicer3 QUIET NO_DEFAULT_PATH)

set(boost_version 1.41.0)
find_package(Boost ${boost_version} REQUIRED COMPONENTS date_time filesystem system)
set(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIR} PARENT_SCOPE)
set(BOOST_LIBRARYDIR "${Boost_LIBRARY_DIRS}" PARENT_SCOPE)
set(BOOST_LIBRARYDIR ${BOOST_LIBRARY_DIR})
set(BOOST_INCLUDEDIR ${Boost_INCLUDE_DIR})

if(NOT BOOST_FOUND)
  find_path(Boost_INCLUDE_DIR boost/config.hpp
  /opt/local/Boost/boost_1_39_0
  /usr/include
  /usr/local/include
  )
  include_directories(${Boost_INCLUDE_DIR})

#  find_path(BOOST_LIBRARY_DIR boost_date_time-vc71-mt.lib)
  find_path(BOOST_LIBRARY_DIR libboost_filesystem.a
    /opt/local/Boost/boost_1_39_0_linux64/lib
    /usr/lib64 /usr/lib /opt/local/lib /usr/local/lib /tools/lib
    DOC " directory for boost, where libboost can be found"
  )

  link_directories(${BOOST_LIBRARY_DIR})

  # Win32 has automatic linking of boost libraries so only add boost
  # libaries on unix
  if(NOT WIN32 OR CYGWIN)
    #set(BOOST_LIBRARIES boost_program_options boost_system boost_filesystem)
    #set(BOOST_LIBRARIES boost_system boost_filesystem)
    set(BOOST_LIBRARIES ${Boost_LIBRARIES})
  endif(NOT WIN32 OR CYGWIN)
endif(NOT BOOST_FOUND)

if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()
else(Slicer3_FOUND)

  find_package(ITK REQUIRED)
  if(ITK_FOUND)
    include(${USE_ITK_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  endif(ITK_FOUND)

  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    include(${USE_VTK_FILE})
  else(VTK_FOUND)
    message(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
  endif(VTK_FOUND)

  find_package(GenerateCLP REQUIRED)
  if(GenerateCLP_FOUND)
    include(${GenerateCLP_USE_FILE})
  else(GenerateCLP_FOUND)
    message(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
  endif(GenerateCLP_FOUND)

endif(Slicer3_FOUND)

   set(shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.cxx MANCOVA.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx  miscMANCOVA.cxx nonMANCOVA.cxx )
  GENERATECLP(shapeAnalysisMANCOVA_SRCS shapeAnalysisMANCOVA.xml)

  add_executable(shapeAnalysisMANCOVA ${shapeAnalysisMANCOVA_SRCS})

  target_link_libraries(shapeAnalysisMANCOVA
    ${ITK_LIBRARIES}
    vtkGraphics
    vtkIO
    vtkCommon
    vtkFiltering
    ${Boost_LIBRARIES}
    )

if(Slicer3_FOUND)
  slicer3_set_plugins_output_path(shapeAnalysisMANCOVA)
  slicer3_install_plugins(shapeAnalysisMANCOVA)
endif(Slicer3_FOUND)

install(TARGETS
  shapeAnalysisMANCOVA
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

