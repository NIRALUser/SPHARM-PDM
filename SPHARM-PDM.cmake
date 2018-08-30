project(SPHARM-PDM)

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

#-----------------------------------------------------------------------------
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#-----------------------------------------------------------------------------
set(CLAPACK_LIBRARY_DIRECTORIES
  ${CLAPACK_DIR}/F2CLIBS/libf2c
  ${CLAPACK_DIR}/BLAS/SRC
  ${CLAPACK_DIR}/SRC
   )
if(WIN32)
  set(CLAPACK_LIBRARIES lapack blas libf2c)
else()
  set(CLAPACK_LIBRARIES lapack blas f2c)
endif()

link_directories( ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} )
link_directories(${CLAPACK_LIBRARY_DIRECTORIES})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/IO)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/SpatialObject)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/Algorithms)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/Numerics)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/Shape/Statistics)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Libraries/SparseLibMVIml)

option(BUILD_LIBRARIES "Build libraries" ON)
mark_as_advanced(BUILD_LIBRARIES)
if(BUILD_LIBRARIES)
  add_subdirectory(Libraries)
endif()

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
option(BUILD_TESTING "Build testing" OFF)
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(Testing)
endif(BUILD_TESTING)
