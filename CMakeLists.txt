project(Spharm_PDM_WRAPPER)

cmake_minimum_required(VERSION 2.8)


#-----------------------------------------------------------------------------
# Superbuild Option - Enabled by default
#-----------------------------------------------------------------------------
option(spharmpdm_SUPERBUILD "Build sparmpdm and the projects it depends on via SuperBuild.cmake." ON)
mark_as_advanced(spharmpdm_SUPERBUILD)
#-----------------------------------------------------------------------------
# Superbuild script
#-----------------------------------------------------------------------------
if(spharmpdm_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
else()
  include("${CMAKE_CURRENT_SOURCE_DIR}/spharm-pdm.cmake")
endif()

