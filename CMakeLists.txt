cmake_minimum_required(VERSION 3.5.0)
if(POLICY CMP0083)
  cmake_policy(SET CMP0083 NEW)
endif()

# C++ standard
if(NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 11)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Update CMake module path
set(SPHARM-PDM_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
list(INSERT CMAKE_MODULE_PATH 0 ${SPHARM-PDM_CMAKE_DIR})
add_subdirectory(${SPHARM-PDM_CMAKE_DIR})
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_BINARY_DIR}/CMake")

# Sanity checks
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)

# Superbuild
set(SUPERBUILD_TOPLEVEL_PROJECT SPHARM-PDM)
set(EXTERNAL_PROJECT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild)
set(EXTERNAL_PROJECT_FILE_PREFIX "External_")
include(ExternalProject)
include(ExternalProjectDependency)

project(SPHARM-PDM VERSION "2.0.0")
set(LOCAL_PROJECT_NAME SPHARM-PDM)
set(EXTENSION_NAME "${LOCAL_PROJECT_NAME}")

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

include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake) #<-- All feature options for top superbuild and inner product build

if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake") #<-- Harness for managing top superbuild, and finally delegate inner product build
  return()
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/${LOCAL_PROJECT_NAME}.cmake") #<-- inner product build (a separate build environment from top level build)
  return()
endif()
