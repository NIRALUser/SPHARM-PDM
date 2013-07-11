project(spharmpdm)

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)
#-----------------------------------------------------------------------------
find_package(ITK REQUIRED)
include(${USE_ITK_FILE})

#-----------------------------------------------------------------------------
find_package(VTK REQUIRED)
include(${USE_VTK_FILE})

#-----------------------------------------------------------------------------
find_package(SlicerExecutionModel REQUIRED GenerateCLP)
include(${GenerateCLP_USE_FILE})
include(${SlicerExecutionModel_USE_FILE})
include(${SlicerExecutionModel_CMAKE_DIR}/SEMMacroBuildCLI.cmake)

#-----------------------------------------------------------------------------
find_package(BatchMake REQUIRED)
include(${BatchMake_USE_FILE})

#-----------------------------------------------------------------------------

set(CLAPACK_LIBRARY_DIRECTORIES
  #${CLAPACK_DIR}/SRC
  #${CLAPACK_DIR}/BLAS/SRC
  ${CLAPACK_DIR}/F2CLIBS/libf2c)
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
include_directories(${Boost_DIR}/include)
include_directories(${Boost_DIR}/include)

link_directories(${Boost_DIR}/lib)
add_subdirectory(Libraries)
add_subdirectory(Applications)

FILE( GLOB list_files ${CMAKE_CURRENT_SOURCE_DIR}/bmm/*.* )
FILE( COPY ${list_files} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

#======================================================================================
  # Create Tests
 include(CTest)
 add_subdirectory(Testing) # contains a CMakeLists.txt
