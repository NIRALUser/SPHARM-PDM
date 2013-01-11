
cmake_minimum_required(VERSION 2.8.3)
cmake_policy(VERSION 2.8.3)

set (LOCAL_PROJECT_NAME SPHARM-PDM)
set(Slicer_SKIP_PROJECT_COMMAND ON) # Useful for Slicer Extensions

option( SPHARM-PDM_BUILD_SLICER_EXTENSION  "Build SPHARM-PDM as a Slicer extension" OFF )

if (SPHARM-PDM_BUILD_SLICER_EXTENSION)
  set(EXTENSION_NAME SPHARM-PDM)
  set(EXTENSION_HOMEPAGE "http://www.nitrc.org/projects/spharm-pdm")
  set(EXTENSION_CATEGORY "Shape Analysis")
  set(EXTENSION_CONTRIBUTORS "Beatriz Paniagua (UNC), Francois Budin(UNC), Martin Styner (UNC)")
  set(EXTENSION_DESCRIPTION " SPHARM-PDM is a tool that computes point-based models using a parametric boundary description for the computing of Shape Analysis.")
  set(EXTENSION_ICONURL "http://www.na-mic.org/Wiki/images/a/ad/Spharm-pdm-icon.png")
  set(EXTENSION_SCREENSHOTURLS "http://www.na-mic.org/Wiki/images/3/34/Spharm-pdm-snapshot.png")
  set(EXTENSION_BUILD_SUBDIRECTORY SPHARM-PDM-build)

  set(MODULE_NAME SPHARM-PDM)

  # SlicerExecutionModel_DEFAULT_CLI_RUNTIME_OUTPUT_DIRECTORY and SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION defined in Slicer_USE_FILE
  # If build as Slicer Extension, CMAKE_INSTALL_PREFIX is set to Extensions/SPHARM-PDM
  # SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION is
  # ./lib/Slicer4.2/cli_module -> not install not cli modules in it
  set( USE_SYSTEM_BatchMake ON )
  set( USE_SYSTEM_VTK ON )

else()

   set( USE_SYSTEM_BatchMake OFF )
   set( USE_SYSTEM_VTK OFF )

endif()

#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#-----------------------------------------------------------------------------
option( SPHARM-PDM_USE_SUPERBUILD "Build SPHARM-PDM and the projects it depends on via SuperBuild.cmake." ON)
#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(SPHARM-PDM_USE_SUPERBUILD)
  project(SuperBuild_SPHARM-PDM)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  project(SPHARM-PDM)
  include("${CMAKE_CURRENT_SOURCE_DIR}/spharm-pdm.cmake")
endif()
