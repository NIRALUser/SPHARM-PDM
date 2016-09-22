cmake_minimum_required(VERSION 2.8.3)
cmake_policy(VERSION 2.8.3)

set(LOCAL_PROJECT_NAME SPHARM-PDM)

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
option(${LOCAL_PROJECT_NAME}_SUPERBUILD "Build ${LOCAL_PROJECT_NAME} and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(${LOCAL_PROJECT_NAME}_SUPERBUILD)


#----------------------------------------------------------------------------------------
# Common SPHARM-PDM extension variables that need to be declared in top level CMakeLists
#----------------------------------------------------------------------------------------
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  set(EXTENSION_NAME SPHARM-PDM)
  set(EXTENSION_HOMEPAGE "http://www.nitrc.org/projects/spharm-pdm")
  set(EXTENSION_CATEGORY "Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "Beatriz Paniagua (UNC), Francois Budin (UNC), Martin Styner (UNC)")
  set(EXTENSION_DESCRIPTION "SPHARM-PDM is a tool that computes point-based models using a parametric boundary description for the computing of Shape Analysis.")
  set(EXTENSION_ICONURL "http://www.na-mic.org/Wiki/images/a/ad/Spharm-pdm-icon.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.na-mic.org/Wiki/images/3/34/Spharm-pdm-snapshot.png")
  set(EXTENSION_BUILD_SUBDIRECTORY . )
  set(MODULE_NAME SPHARM-PDM)
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
  return()
endif()

message(FATAL_ERROR "You should never reach this point !")

