
#-----------------------------------------------------------------------------
# Set a default build type if none was specified
#-----------------------------------------------------------------------------
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
# Pass variables to dependent projects
if(NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Configuring with build type '${CMAKE_BUILD_TYPE}'")
  mark_as_superbuild(VARS CMAKE_BUILD_TYPE ALL_PROJECTS)
else()
  mark_as_superbuild(VARS CMAKE_CONFIGURATION_TYPES ALL_PROJECTS)
endif()

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
mark_as_superbuild(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION:BOOL)

#-----------------------------------------------------------------------------
# Extension meta-information
#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.nitrc.org/projects/spharm-pdm")
set(EXTENSION_CATEGORY "SPHARM")
set(EXTENSION_CONTRIBUTORS "Beatriz Paniagua (UNC), Francois Budin (UNC), Martin Styner (UNC), Laura Pascal (Kitware), Hina Shah (Kitware)")
set(EXTENSION_DESCRIPTION "SPHARM-PDM is a tool that computes point-based models using a parametric boundary description for the computing of Shape Analysis.")
set(EXTENSION_ICONURL "http://www.na-mic.org/Wiki/images/a/ad/Spharm-pdm-icon.png")
set(EXTENSION_SCREENSHOTURLS "http://www.na-mic.org/Wiki/images/3/34/Spharm-pdm-snapshot.png")
set(EXTENSION_DEPENDS MeshToLabelMap) # Specified as a space separated string, a list or 'NA' if any

set(EXTENSION_BUILD_SUBDIRECTORY ${LOCAL_PROJECT_NAME}-build)

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------

find_package(Git REQUIRED)
mark_as_superbuild(GIT_EXECUTABLE)

set(${LOCAL_PROJECT_NAME}_CLI_EXECUTABLE_LINK_FLAGS "")

if(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  find_package(Slicer 4.9 REQUIRED)
  include(${Slicer_USE_FILE})
  mark_as_superbuild(Slicer_DIR)

  # Link flags used in "spharm_add_executable" function
  # On macOS, these are of the form "-Wl,-rpath,@loader_path/..."
  set(${LOCAL_PROJECT_NAME}_CLI_EXECUTABLE_LINK_FLAGS "${Slicer_INSTALL_THIRDPARTY_EXECUTABLE_LINK_FLAGS}")
endif()

#-----------------------------------------------------------------------------
# Build option(s)
#-----------------------------------------------------------------------------
set(_default ON)
if(DEFINED Slicer_DIR)
  set(_default OFF)
endif()
option(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT "Install development support include and libraries for external packages." ${_default})
mark_as_advanced(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
mark_as_superbuild(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)

option(BUILD_TESTING "Build testing" OFF)
mark_as_superbuild(BUILD_TESTING)

option(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
mark_as_superbuild(COMPILE_MetaMeshTools)

option(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
mark_as_superbuild(COMPILE_SegPostProcessCLP)

option(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
mark_as_superbuild(COMPILE_GenParaMeshCLP)

option(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
mark_as_superbuild(COMPILE_ParaToSPHARMMeshCLP)

option(COMPILE_SpharmTool "Compile SpharmTool." ON)
mark_as_superbuild(COMPILE_SpharmTool)

#-----------------------------------------------------------------------------
# Platform check
#-----------------------------------------------------------------------------
if(APPLE)
  # See CMake/Modules/Platform/Darwin.cmake)
  #   6.x == Mac OSX 10.2 (Jaguar)
  #   7.x == Mac OSX 10.3 (Panther)
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  #  11.x == Mac OSX 10.7 (Lion)
  #  12.x == Mac OSX 10.8 (Mountain Lion)
  if (DARWIN_MAJOR_VERSION LESS "9")
    message(FATAL_ERROR "Only Mac OSX >= 10.5 are supported !")
  endif()

  if(NOT DEFINED CMAKE_MACOSX_RPATH)
    set(CMAKE_MACOSX_RPATH ON)
  endif()
  mark_as_superbuild(CMAKE_MACOSX_RPATH:BOOL)
endif()

#-----------------------------------------------------------------------------
# Output directories
foreach(type LIBRARY RUNTIME ARCHIVE)
  # Make sure the directory exists
  if(DEFINED ${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY
     AND NOT EXISTS ${${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY})
    message(FATAL_ERROR "${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY is set to a non-existing directory [${${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY}]")
  endif()

  if(${LOCAL_PROJECT_NAME}_SUPERBUILD)
    set(output_dir ${${LOCAL_PROJECT_NAME}_BINARY_DIR}/bin)
    if(NOT DEFINED ${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY)
      set(${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY ${${LOCAL_PROJECT_NAME}_BINARY_DIR}/${LOCAL_PROJECT_NAME}-build/bin)
    endif()
    mark_as_superbuild(${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY:PATH)
    file(MAKE_DIRECTORY ${${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY})
  else()
    if(NOT DEFINED ${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY)
      set(output_dir ${${LOCAL_PROJECT_NAME}_BINARY_DIR}/bin)
    else()
      set(output_dir ${${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY})
    endif()
  endif()
  set(CMAKE_${type}_OUTPUT_DIRECTORY ${output_dir} CACHE INTERNAL "Single output directory for building all libraries.")
endforeach()

#-----------------------------------------------------------------------------
if(NOT COMMAND SETIFEMPTY)
  macro(SETIFEMPTY)
    set(KEY ${ARGV0})
    set(VALUE ${ARGV1})
    if(NOT ${KEY})
      set(${ARGV})
    endif()
  endmacro()
endif()

#-------------------------------------------------------------------------
# Install directories
if(NOT ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION lib)
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_INSTALL_ARCHIVE_DESTINATION lib)
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION bin)
else()
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_INSTALL_ARCHIVE_DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION ${Slicer_INSTALL_THIRDPARTY_BIN_DIR})
endif()

#-------------------------------------------------------------------------
# Output and install directories for:
# * configuring SlicerExecutionModel external project default output and install directories
#   when building as a standalone project.
# * configuring CLI-like executables that are not using the SlicerExecutionModel.
SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION})
SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_ARCHIVE_DESTINATION})
SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION})

foreach(type IN ITEMS LIBRARY ARCHIVE RUNTIME)
  mark_as_superbuild(VARS
    CMAKE_INSTALL_${type}_DESTINATION:PATH
    ${LOCAL_PROJECT_NAME}_CLI_${type}_OUTPUT_DIRECTORY
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_${type}_DESTINATION
    )
endforeach()

#-------------------------------------------------------------------------
# Augment compiler flags
#-------------------------------------------------------------------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(ITKSetStandardCompilerFlags)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS}" )
else() # Release, or anything else
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS}" )
endif()

