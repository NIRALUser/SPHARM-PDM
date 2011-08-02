project(SPHARMPDM)

project(ShapeAnalysisModule)

cmake_minimum_required(VERSION 2.6)

find_package(Slicer3 QUIET NO_DEFAULT_PATH)

if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

else(Slicer3_FOUND)
  set(LIBRARY_OUTPUT_PATH ${SPHARMPDM_BINARY_DIR}/bin CACHE PATH "Single output directory for building all libraries.")
 set(EXECUTABLE_OUTPUT_PATH ${SPHARMPDM_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")

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

include(${SPHARMPDM_SOURCE_DIR}/CMake/lapackSetup.cmake)
include(${SPHARMPDM_SOURCE_DIR}/CMake/boostSetup.cmake)

include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/Shape/IO)
include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/Shape/Algorithms)
include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/Shape/SpatialObject)
include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/Shape/Numerics)
include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/Shape/Statistics)

include_directories(${SPHARMPDM_SOURCE_DIR}/Libraries/SparseLibMVIml)

subdirs(Libraries)
subdirs(Applications )

