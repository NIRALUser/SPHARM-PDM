if(Slicer3_DIR)
  find_package(Slicer3 REQUIRED NO_DEFAULT_PATH)
  include(${Slicer3_USE_FILE})
  get_filename_component(Slicer3_BUILD_DIR ${Slicer3_USE_FILE} PATH CACHE)
  include(${Slicer3_BUILD_DIR}/Slicer3Config.cmake)
  slicer3_set_default_install_prefix_for_external_projects()
endif()

find_package(ITK REQUIRED)
message(STATUS "###################>${ITK_FOUND}<####")
if(ITK_FOUND)
  message(STATUS "###################>${USE_ITK_FILE}<####")
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

set(CLAPACK_LIBRARY_DIRECTORIES
  ${CLAPACK_DIR}/SRC
  ${CLAPACK_DIR}/BLAS/SRC
  ${CLAPACK_DIR}/F2CLIBS/libf2c)
set(CLAPACK_LIBRARIES lapack blas f2c)

link_directories(${CLAPACK_LIBRARY_DIRECTORIES})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/IO)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/SpatialObject)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/Algorithms)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/Numerics)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/Statistics)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/SparseLibMVIml)

message("Boost_DIR=${Boost_DIR}")
include_directories(${Boost_DIR}/include)
link_directories(${Boost_DIR}/lib)
add_subdirectory(Libraries)
add_subdirectory(Applications)
