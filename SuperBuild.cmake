#-----------------------------------------------------------------------------
#If it is build as an extension
#-----------------------------------------------------------------------------
if(${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  if( APPLE )
    set( CMAKE_EXE_LINKER_FLAGS -Wl,-rpath,@loader_path/../../../../../ )
  endif()
endif()

#-----------------------------------------------------------------------------
# Superbuild option(s)
#-----------------------------------------------------------------------------
option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
set(Slicer_USE_SYSTEM_ITK ${USE_SYSTEM_ITK})

option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
set(Slicer_USE_SYSTEM_SlicerExecutionModel ${USE_SYSTEM_SlicerExecutionModel})

option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
set(Slicer_USE_SYSTEM_VTK ${USE_SYSTEM_VTK})

option(USE_SYSTEM_LAPACK "Build using an externally defined version of LAPACK" OFF)
set(Slicer_USE_SYSTEM_LAPACK ${USE_SYSTEM_LAPACK})

#-----------------------------------------------------------------------------
# Configure "external" projects
#-----------------------------------------------------------------------------
set(LAPACK_INSTALL_RUNTIME_DIR ${${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION})
set(LAPACK_INSTALL_LIBRARY_DIR ${${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION})

#------------------------------------------------------------------------------
# Project dependencies
#------------------------------------------------------------------------------

set(${LOCAL_PROJECT_NAME}_DEPENDS
  LAPACK
  )
if(NOT SPHARM-PDM_BUILD_SLICER_EXTENSION)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDS
    ITK
    SlicerExecutionModel
    VTK
    )
endif()

#------------------------------------------------------------------------------
# Superbuild-type bundled extensions
#------------------------------------------------------------------------------

# The following logic is documented in the "Bundle extensions adding source directories"
# section found in the top-level CMakeLists.txt

set(_all_extension_depends )

# Build only inner-build for superbuild-type extensions
set(Slicer_BUNDLED_EXTENSION_NAMES)
foreach(extension_dir ${Slicer_EXTENSION_SOURCE_DIRS})
  get_filename_component(extension_dir ${extension_dir} ABSOLUTE)
  get_filename_component(extension_name ${extension_dir} NAME) # The assumption is that source directories are named after the extension project
  if(EXISTS ${extension_dir}/SuperBuild OR EXISTS ${extension_dir}/Superbuild)
    set(${extension_name}_SUPERBUILD 0)
    mark_as_superbuild(${extension_name}_SUPERBUILD:BOOL)

    if(NOT DEFINED ${extension_name}_EXTERNAL_PROJECT_EXCLUDE_ALL)
      set(${extension_name}_EXTERNAL_PROJECT_EXCLUDE_ALL FALSE)
    endif()
    if(NOT ${extension_name}_EXTERNAL_PROJECT_EXCLUDE_ALL)
      list(APPEND EXTERNAL_PROJECT_ADDITIONAL_DIRS "${extension_dir}/SuperBuild")
      list(APPEND EXTERNAL_PROJECT_ADDITIONAL_DIRS "${extension_dir}/Superbuild")
    endif()
    if(NOT DEFINED ${extension_name}_EXTERNAL_PROJECT_DEPENDENCIES)
      set(${extension_name}_EXTERNAL_PROJECT_DEPENDENCIES )
    endif()

    set(_external_project_cmake_files)

    # SuperBuild
    file(GLOB _external_project_cmake_files1 RELATIVE "${extension_dir}/SuperBuild" "${extension_dir}/SuperBuild/External_*.cmake")
    list(APPEND _external_project_cmake_files ${_external_project_cmake_files1})

    # Superbuild
    file(GLOB _external_project_cmake_files2 RELATIVE "${extension_dir}/Superbuild" "${extension_dir}/Superbuild/External_*.cmake")
    list(APPEND _external_project_cmake_files ${_external_project_cmake_files2})

    list(REMOVE_DUPLICATES _external_project_cmake_files)

    set(_extension_depends)
    set(_msg_extension_depends)
    foreach (_external_project_cmake_file ${_external_project_cmake_files})
      string(REGEX MATCH "External_(.+)\.cmake" _match ${_external_project_cmake_file})
      set(_additional_project_name "${CMAKE_MATCH_1}")
      if(${extension_name}_EXTERNAL_PROJECT_EXCLUDE_ALL)
        set(_include FALSE)
      else()
        set(_include TRUE)
        if(NOT "${${extension_name}_EXTERNAL_PROJECT_DEPENDENCIES}" STREQUAL "")
          list(FIND ${extension_name}_EXTERNAL_PROJECT_DEPENDENCIES ${_additional_project_name} _index)
          if(_index EQUAL -1)
            set(_include FALSE)
          endif()
        endif()
      endif()
      if(_include)
          list(APPEND _extension_depends ${_additional_project_name})
          list(APPEND _msg_extension_depends ${_additional_project_name})
      else()
        list(APPEND _msg_extension_depends "exclude(${_additional_project_name})")
      endif()
    endforeach()

    list(APPEND Slicer_BUNDLED_EXTENSION_NAMES ${extension_name})

    message(STATUS "SuperBuild - ${extension_name} extension => ${_msg_extension_depends}")

    list(APPEND _all_extension_depends ${_extension_depends})
  endif()
endforeach()

if(_all_extension_depends)
  list(REMOVE_DUPLICATES _all_extension_depends)
endif()

list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDS ${_all_extension_depends})

mark_as_superbuild(Slicer_BUNDLED_EXTENSION_NAMES:STRING)

#-----------------------------------------------------------------------------
# Top-level "external" project
#-----------------------------------------------------------------------------

# Extension dependencies
foreach(dep ${EXTENSION_DEPENDS})
  mark_as_superbuild(${dep}_DIR)
endforeach()

set(proj ${SUPERBUILD_TOPLEVEL_PROJECT})

ExternalProject_Include_Dependencies(${proj}
  PROJECT_VAR proj
  SUPERBUILD_VAR ${LOCAL_PROJECT_NAME}_SUPERBUILD
  )

ExternalProject_Add(${proj}
  ${${proj}_EP_ARGS}
  DOWNLOAD_COMMAND ""
  INSTALL_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${EXTENSION_BUILD_SUBDIRECTORY}
  CMAKE_CACHE_ARGS
    # Compiler settings
    -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
    -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
    -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
    -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
    -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
    # Packaging
    -DMIDAS_PACKAGE_EMAIL:STRING=${MIDAS_PACKAGE_EMAIL}
    -DMIDAS_PACKAGE_API_KEY:STRING=${MIDAS_PACKAGE_API_KEY}
    # Superbuild
    -D${LOCAL_PROJECT_NAME}_SUPERBUILD:BOOL=OFF
    -DEXTENSION_SUPERBUILD_BINARY_DIR:PATH=${${LOCAL_PROJECT_NAME}_BINARY_DIR}
    -DSlicer_EXTENSION_SOURCE_DIRS:STRING=${Slicer_EXTENSION_SOURCE_DIRS}
  INSTALL_COMMAND ""
  DEPENDS
    ${${LOCAL_PROJECT_NAME}_DEPENDS}
  )

ExternalProject_AlwaysConfigure(${proj})
