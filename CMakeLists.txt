cmake_minimum_required(VERSION 3.5.0)

# C++ standard
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Update CMake module path
set(SPHARM-PDM_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
list(INSERT CMAKE_MODULE_PATH 0 ${SPHARM-PDM_CMAKE_DIR})

# Sanity checks
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)

project(SPHARM-PDM VERSION "2.0.0")
set(LOCAL_PROJECT_NAME SPHARM-PDM)

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
option(${LOCAL_PROJECT_NAME}_SUPERBUILD "Build ${LOCAL_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${LOCAL_PROJECT_NAME}_SUPERBUILD)

#-----------------------------------------------------------------------------
# Standalone vs Slicer extension option
#-----------------------------------------------------------------------------

# This option should be named after the project name, it corresponds to the
# option set to ON when the project is build by the Slicer Extension build
# system.

set(_default OFF)
set(_reason "${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION is ON")
if(NOT DEFINED ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION AND DEFINED Slicer_DIR)
  set(_default ON)
  set(_reason "Slicer_DIR is SET")
endif()

option(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION "Build as a Slicer Extension" ${_default})

set(_msg "Checking if building as a Slicer extension")
message(STATUS ${_msg})
if(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  message(STATUS "${_msg} - yes (${_reason})")
else()
  message(STATUS "${_msg} - no (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION is OFF)")
endif()

#----------------------------------------------------------------------------------------
# Slicer
#----------------------------------------------------------------------------------------
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  set(EXTENSION_NAME SPHARM-PDM)
  set(EXTENSION_HOMEPAGE "http://www.nitrc.org/projects/spharm-pdm")
  set(EXTENSION_CATEGORY "SPHARM")
  set(EXTENSION_CONTRIBUTORS "Beatriz Paniagua (UNC), Francois Budin (UNC), Martin Styner (UNC), Laura Pascal (Kitware), Hina Shah (Kitware)")
  set(EXTENSION_DESCRIPTION "SPHARM-PDM is a tool that computes point-based models using a parametric boundary description for the computing of Shape Analysis.")
  set(EXTENSION_ICONURL "http://www.na-mic.org/Wiki/images/a/ad/Spharm-pdm-icon.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.na-mic.org/Wiki/images/3/34/Spharm-pdm-snapshot.png")
  set(EXTENSION_DEPENDS MeshToLabelMap) # Specified as a space separated string, a list or 'NA' if any
  set(EXTENSION_BUILD_SUBDIRECTORY ${LOCAL_PROJECT_NAME}-build)
  set(SUPERBUILD_TOPLEVEL_PROJECT ${LOCAL_PROJECT_NAME})
  set(MODULE_NAME SPHARM-PDM)

  # Slicer
  find_package(Slicer REQUIRED)
  include(${Slicer_USE_FILE})
  mark_as_superbuild(Slicer_DIR)
endif()

#-----------------------------------------------------------------------------
# Common build features for both the superbuild and the main build
#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake) #<-- All feature options for top superbuild and inner product build

#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake") #<-- Harness for managing top superbuild, and finally delegate inner product build
  return()
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake") #<-- inner product build (a separate build environment from top level build)
  return()
endif()
