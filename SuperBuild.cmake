#-----------------------------------------------------------------------------
#If it is build as an extension
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  if( APPLE )
    set( CMAKE_EXE_LINKER_FLAGS -Wl,-rpath,@loader_path/../../../../../ )
  endif()
endif()

#-----------------------------------------------------------------------------
# Superbuild option(s)
#-----------------------------------------------------------------------------
option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
set(Slicer_USE_SYSTE_ITKv4 ${USE_SYSTEM_ITK})

option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
set(Slicer_USE_SYSTEM_SlicerExecutionModel ${USE_SYSTEM_SlicerExecutionModel})

option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
set(Slicer_USE_SYSTEM_VTK ${USE_SYSTEM_VTK})

option(USE_SYSTEM_CLAPACK "Build using an externally defined version of CLAPACK" OFF)
set(Slicer_USE_SYSTEM_CLAPACK ${USE_SYSTEM_CLAPACK})

#-----------------------------------------------------------------------------
# Top-level "external" project
#-----------------------------------------------------------------------------

foreach(dep ${EXTENSION_DEPENDS})
  mark_as_superbuild(${dep}_DIR)
endforeach()

# Project dependencies
set(${LOCAL_PROJECT_NAME}_DEPENDS
  CLAPACK
  )
if(NOT SPHARM-PDM_BUILD_SLICER_EXTENSION)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDS
    ITKv4
    SlicerExecutionModel
    VTK
    )
endif()

set(proj ${SUPERBUILD_TOPLEVEL_PROJECT})

ExternalProject_Include_Dependencies(${proj}
  PROJECT_VAR proj
  SUPERBUILD_VAR ${PROJECT_NAME}_SUPERBUILD
  )

ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  DOWNLOAD_COMMAND ""
  INSTALL_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${EXTENSION_BUILD_SUBDIRECTORY}
  CMAKE_CACHE_ARGS
    # Compiler settings
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
    # Packaging
    -DMIDAS_PACKAGE_EMAIL:STRING=${MIDAS_PACKAGE_EMAIL}
    -DMIDAS_PACKAGE_API_KEY:STRING=${MIDAS_PACKAGE_API_KEY}
    # Superbuild
    -D${EXTENSION_NAME}_SUPERBUILD:BOOL=OFF
    -DEXTENSION_SUPERBUILD_BINARY_DIR:PATH=${${EXTENSION_NAME}_BINARY_DIR}
  INSTALL_COMMAND ""
  DEPENDS
    ${${LOCAL_PROJECT_NAME}_DEPENDS}
  )

ExternalProject_AlwaysConfigure(${proj})
