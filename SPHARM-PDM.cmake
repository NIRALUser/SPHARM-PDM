
#-----------------------------------------------------------------------------
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(LAPACKE CONFIG REQUIRED)

# Set Fortran_<id>_RUNTIME_LIBRARIES and CMAKE_Fortran_IMPLICIT_LINK_*
set(Fortran_COMPILER_ID "${LAPACKE_Fortran_COMPILER_ID}")
find_package(Fortran REQUIRED)

# --------------------------------------------------------------------------
# Bundle extensions adding source directories.
# --------------------------------------------------------------------------

#
# Support for bundling "SuperBuild-type" extension:
#
# * An extension is considered to be of type "SuperBuild" if a directory
#   "<extension_dir>/SuperBuild" or "<extension_dir>/Superbuild" exists.
#   Corresponding directory is appended to EXTERNAL_PROJECT_ADDITIONAL_DIRS.
#
# * If variable "<extension_name>_EXTERNAL_PROJECT_EXCLUDE_ALL" is set to TRUE, corresponding SuperBuild directory
#   is not appended to EXTERNAL_PROJECT_ADDITIONAL_DIRS.
#
# * Associated external projects are globbed using expression of the form
#   "<extension_dir>/(SuperBuild|Superbuild)/External_*.cmake".
#
# * List of external project names is extracted from the "External_<projectName>.cmake"
#   files and appended to Slicer_DEPENDENCIES. This ensures they are build before Slicer inner build.
#   Setting variable "<extension_name>_EXTERNAL_PROJECT_DEPENDENCIES" to a list of <projectName> allows
#   to override the list of <projectName> appended to Slicer_DEPENDENCIES.
#
# * Variable Slicer_BUNDLED_EXTENSION_NAMES is updated with the names of external project
#   and passed to Slicer inner build. It is then used in SlicerCPack. to package associated
#   external projects if the cache variable <extensionName>_CPACK_INSTALL_CMAKE_PROJECTS
#   was defined in the extension sources.
#
# Corresponding logic is implemented in SuperBuild.cmake
#

set(extensions_build_dir "${CMAKE_BINARY_DIR}/E")

function(_add_extension_source_dir extension_source_dir what)
  get_filename_component(extension_source_dir ${extension_source_dir} REALPATH)
  get_filename_component(extension_source_dirname ${extension_source_dir} NAME_WE)
  message(STATUS "--------------------------------------------------")
  message(STATUS "Configuring ${what}: ${extension_source_dirname}")
  set(ExternalData_SOURCE_ROOT ${extension_source_dir})
  set(${extension_source_dirname}_SOURCE_DIR ${extension_source_dir})
  set(${extension_source_dirname}_BINARY_DIR ${extensions_build_dir}/${extension_source_dirname})
  add_subdirectory(
    ${${extension_source_dirname}_SOURCE_DIR}
    ${${extension_source_dirname}_BINARY_DIR}
    )
endfunction()

if(NOT Slicer_SOURCE_DIR)
  foreach(extension_source_dir ${Slicer_EXTENSION_SOURCE_DIRS})
    _add_extension_source_dir(${extension_source_dir} "extension directory")
  endforeach()
endif()

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
if(SPHARM-PDM_BUILD_TESTING)
  include(CTest)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
# Packaging
#-----------------------------------------------------------------------------
set(EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS)
list(APPEND EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS "${LAPACK_DIR};LAPACK;RuntimeLibraries;/")
if(NOT Slicer_SOURCE_DIR)
  foreach(extension_item IN LISTS Slicer_BUNDLED_EXTENSION_NAMES)
    list(APPEND EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS ${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS})
  endforeach()
endif()
set(${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS "${EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS}" CACHE STRING "List of external projects to install" FORCE)

# Install fortran runtime libraries
if(DEFINED Slicer_DIR)
  if(NOT APPLE)
    Fortran_InstallLibrary(
      FILES ${Fortran_${LAPACKE_Fortran_COMPILER_ID}_RUNTIME_LIBRARIES}
      DESTINATION ${Slicer_INSTALL_THIRDPARTY_LIB_DIR} COMPONENT RuntimeLibraries
      )
  else()
    set(${EXTENSION_NAME}_FIXUP_BUNDLE_LIBRARY_DIRECTORIES ${Fortran_${LAPACKE_Fortran_COMPILER_ID}_RUNTIME_DIRECTORIES} CACHE STRING "List of fixup bundle library directories" FORCE)
  endif()
endif()

#-----------------------------------------------------------------------------
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};Runtime;/")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};RuntimeLibraries;/")
if(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};Development;/")
endif()
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS}")

#-----------------------------------------------------------------------------
if(DEFINED Slicer_DIR)
  include(${Slicer_EXTENSION_GENERATE_CONFIG})
  include(${Slicer_EXTENSION_CPACK})
endif()
