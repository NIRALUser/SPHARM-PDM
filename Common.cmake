
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#-----------------------------------------------------------------------------
# Build option(s)
#-----------------------------------------------------------------------------
set(PRIMARY_PROJECT_NAME ${LOCAL_PROJECT_NAME})

option(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT "Install development support include and libraries for external packages." OFF)
mark_as_advanced(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)

option(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
option(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
option(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
option(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
option(COMPILE_SpharmTool "Compile SpharmTool." ON)

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH
  ${${PROJECT_NAME}_SOURCE_DIR}/CMake
  ${${PROJECT_NAME}_BINARY_DIR}/CMake
  ${CMAKE_MODULE_PATH}
  )
#-----------------------------------------------------------------------------
# Sanity checks
#------------------------------------------------------------------------------
include(PreventInSourceBuilds)
include(PreventInBuildInstalls)

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
find_package(Git)
if(NOT GIT_FOUND)
  message(WARNING "Git may be needed to download external dependencies: Install Git and try to re-configure")
endif()

option(USE_GIT_PROTOCOL "If behind a firewall turn this off to use http instead." ON)
if(NOT USE_GIT_PROTOCOL)
  set(git_protocol "http")
else()
  set(git_protocol "git")
endif()

#-----------------------------------------------------------------------------
# Platform check
#-----------------------------------------------------------------------------
set(PLATFORM_CHECK true)
if(PLATFORM_CHECK)
  # See CMake/Modules/Platform/Darwin.cmake)
  #   6.x == Mac OSX 10.2 (Jaguar)
  #   7.x == Mac OSX 10.3 (Panther)
  #   8.x == Mac OSX 10.4 (Tiger)
  #   9.x == Mac OSX 10.5 (Leopard)
  #  10.x == Mac OSX 10.6 (Snow Leopard)
  if (DARWIN_MAJOR_VERSION LESS "9")
    message(FATAL_ERROR "Only Mac OSX >= 10.5 are supported !")
  endif()
endif()

#-----------------------------------------------------------------------------
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

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
      set(${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY ${${LOCAL_PROJECT_NAME}_BINARY_DIR}/${${LOCAL_PROJECT_NAME}}-build/bin)
    endif()
    # The variable is manually appended to ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
  else()
    if(NOT DEFINED ${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY)
      set(output_dir ${${LOCAL_PROJECT_NAME}_BINARY_DIR}/bin)
    else()
      set(output_dir ${${LOCAL_PROJECT_NAME}_CMAKE_${type}_OUTPUT_DIRECTORY})
    endif()
  endif()
  set(CMAKE_${type}_OUTPUT_DIRECTORY ${output_dir} CACHE INTERNAL "Single output directory for building all libraries.")
endforeach()

#-------------------------------------------------------------------------
# Install directories
if(NOT ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  SETIFEMPTY(CMAKE_INSTALL_LIBRARY_DESTINATION lib)
  SETIFEMPTY(CMAKE_INSTALL_ARCHIVE_DESTINATION lib)
  SETIFEMPTY(CMAKE_INSTALL_RUNTIME_DESTINATION bin)
else()
  SETIFEMPTY(CMAKE_INSTALL_LIBRARY_DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
  SETIFEMPTY(CMAKE_INSTALL_ARCHIVE_DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
  SETIFEMPTY(CMAKE_INSTALL_RUNTIME_DESTINATION ${Slicer_INSTALL_THIRDPARTY_BIN_DIR})
endif()

#-------------------------------------------------------------------------
# SlicerExecutionModel output and install directories
# * for configuring SlicerExecutionModel external project default directories
# * for building and installing regular executables
if(NOT ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION ${CMAKE_INSTALL_LIBRARY_DESTINATION})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION ${CMAKE_INSTALL_ARCHIVE_DESTINATION})
  SETIFEMPTY(${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION ${CMAKE_INSTALL_RUNTIME_DESTINATION})
endif()

#-------------------------------------------------------------------------
if(NOT DEFINED CMAKE_MACOSX_RPATH)
  set(CMAKE_MACOSX_RPATH ON)
endif()

#-------------------------------------------------------------------------
# Augment compiler flags
#-------------------------------------------------------------------------
include(ITKSetStandardCompilerFlags)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS}" )
else() # Release, or anything else
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS}" )
endif()

