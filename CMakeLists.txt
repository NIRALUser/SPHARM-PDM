cmake_minimum_required(VERSION 3.5.0)

project(SPHARM-PDM)
set(LOCAL_PROJECT_NAME SPHARM-PDM)

set(${LOCAL_PROJECT_NAME}_MAJOR_VERSION "2")
set(${LOCAL_PROJECT_NAME}_MINOR_VERSION "0")
set(${LOCAL_PROJECT_NAME}_PATCH_VERSION "0")
set(${LOCAL_PROJECT_NAME}_VERSION
  ${${LOCAL_PROJECT_NAME}_MAJOR_VERSION}.${${LOCAL_PROJECT_NAME}_MINOR_VERSION}.${${LOCAL_PROJECT_NAME}_PATCH_VERSION})

## NOTE THERE SHOULD BE NO PROJECT STATEMENT HERE!
## This file acts as a simple switch to initiate
## two completely independant CMake build environments.

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#                   Phase I:  ${LOCAL_PROJECT_NAME}_SUPERBUILD is set to ON, and the
#                             supporting packages defined in "SuperBuild.cmake"
#                             are built.  The last package in "SuperBuild.cmake"
#                             to be built is a recursive call to this
#                             file with ${LOCAL_PROJECT_NAME}_SUPERBUILD explicitly
#                             set to "OFF" to initiate Phase II
#
#                   Phase II: Build the ${LOCAL_PROJECT_NAME}, referencing the support
#                             packages built in Phase I.
#-----------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
   ${CMAKE_SOURCE_DIR}/CMake
   ${CMAKE_SOURCE_DIR}/SuperBuild
   ${CMAKE_MODULE_PATH}
)

option(${LOCAL_PROJECT_NAME}_SUPERBUILD "Build ${LOCAL_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${LOCAL_PROJECT_NAME}_SUPERBUILD)

#----------------------------------------------------------------------------------------
# Slicer
#----------------------------------------------------------------------------------------
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  set(extension_name_DIR "MeshToLabelMap")

  set(EXTENSION_NAME SPHARM-PDM)
  set(EXTENSION_HOMEPAGE "http://www.nitrc.org/projects/spharm-pdm")
  set(EXTENSION_CATEGORY "SPHARM")
  set(EXTENSION_CONTRIBUTORS "Beatriz Paniagua (UNC), Francois Budin (UNC), Martin Styner (UNC), Laura Pascal (Kitware), Hina Shah (Kitware)")
  set(EXTENSION_DESCRIPTION "SPHARM-PDM is a tool that computes point-based models using a parametric boundary description for the computing of Shape Analysis.")
  set(EXTENSION_ICONURL "http://www.na-mic.org/Wiki/images/a/ad/Spharm-pdm-icon.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.na-mic.org/Wiki/images/3/34/Spharm-pdm-snapshot.png")
  set(EXTENSION_DEPENDS ${extension_name_DIR})
  set(EXTENSION_BUILD_SUBDIRECTORY ${LOCAL_PROJECT_NAME}-build)
  set(SUPERBUILD_TOPLEVEL_PROJECT ${LOCAL_PROJECT_NAME})
  set(MODULE_NAME SPHARM-PDM)

  # Slicer
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  mark_as_superbuild(Slicer_DIR)

  # Slicer Execution Model
  find_package(SlicerExecutionModel REQUIRED)
  include(${SlicerExecutionModel_USE_FILE})
  include(${GenerateCLP_USE_FILE})

endif()

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  project(SuperBuild_${LOCAL_PROJECT_NAME})
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  project(${LOCAL_PROJECT_NAME}) 
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake")
#  return()
endif()

#-----------------------------------------------------------------------------
# Clear ${PROJECT_NAME}Targets.cmake
file(WRITE ${CMAKE_BINARY_DIR}/${Slicer_SHARE_DIR}/${PROJECT_NAME}Targets.cmake "")

#-----------------------------------------------------------------------------
add_subdirectory(Modules)
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
 add_subdirectory(CommandLineTool)
endif()

#-----------------------------------------------------------------------------
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  # CPack
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  include(${Slicer_EXTENSION_CPACK})
endif()
