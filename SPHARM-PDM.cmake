project(SPHARM-PDM)

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)

#-----------------------------------------------------------------------------
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(CLAPACK NO_MODULE REQUIRED)
# Workaround incomplete lapack target
set_target_properties(lapack PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${CLAPACK_DIR}/../CLAPACK/INCLUDE"
  )
set(CLAPACK_LIBRARIES lapack blas f2c)

#-----------------------------------------------------------------------------
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

add_subdirectory(Libraries)
add_subdirectory(Modules)
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
 add_subdirectory(CommandLineTool)
endif()

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
option(BUILD_TESTING "Build testing" OFF)
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(Testing)
endif(BUILD_TESTING)

#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  # CPack
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()
