find_package(Git REQUIRED)
include(ExternalProject)
enable_language(C)
enable_language(CXX)

if(NOT SETIFEMPTY)
macro(SETIFEMPTY)
  set(KEY ${ARGV0})
  set(VALUE ${ARGV1})
  if(NOT ${KEY})
    set(${ARGV})
  endif(NOT ${KEY})
endmacro(SETIFEMPTY KEY VALUE)
endif(NOT SETIFEMPTY)
SETIFEMPTY(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
SETIFEMPTY(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
SETIFEMPTY(CMAKE_BUNDLE_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_BINARY_DIR}/bin)
link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})

#-----------------------------------------------------------------------------
# Update CMake module path
#------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/SuperBuild
  ${CMAKE_BINARY_DIR}/CMake
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake #  CMake directory
  ${CMAKE_CURRENT_SOURCE_DIR}/src/CMake # CMake directory
  ${CMAKE_MODULE_PATH}
  )

include(SlicerMacroEmptyExternalProject)

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
# Update CMake module path
#------------------------------------------------------------------------------

set(CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/CMake
  ${CMAKE_SOURCE_DIR}/SuperBuild
  ${CMAKE_BINARY_DIR}/CMake
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/CMake #  CMake directory
  ${CMAKE_CURRENT_SOURCE_DIR}/src/CMake # CMake directory
  ${CMAKE_MODULE_PATH}
  )

#-----------------------------------------------------------------------------
# Prerequisites
#------------------------------------------------------------------------------
#
# BRAINS4 Addition: install to the common library
# directory, so that all libs/include etc ends up
# in one common tree
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Where all the prerequisite libraries go" FORCE)
set(${CMAKE_PROJECT_NAME}_BUILD_TESTING ON CACHE BOOL "Turn on Testing for BRAINS")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Statically Link Everything")

# Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()


#-------------------------------------------------------------------------
# augment compiler flags
#-------------------------------------------------------------------------
include(CompilerFlagSettings)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_DEBUG_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_DEBUG_DESIRED_FLAGS}" )
else() # Release, or anything else
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C_RELEASE_DESIRED_FLAGS}" )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_RELEASE_DESIRED_FLAGS}" )
endif()


option(USE_NO_SLICER "Build Slicer and the projects it depends on via SuperBuild.cmake." OFF)
#mark_as_advanced(USE_NO_SLICER))

if(NOT USE_NO_SLICER)
  find_package(Slicer3 REQUIRED NO_DEFAULT_PATH)
  include(${Slicer3_USE_FILE})
  get_filename_component(Slicer3_BUILD_DIR ${Slicer3_USE_FILE} PATH CACHE)
  include(${Slicer3_BUILD_DIR}/Slicer3Config.cmake)
  slicer3_set_default_install_prefix_for_external_projects()
  set(SLICER_ARGS -DSlicer3_DIR:PATH=${Slicer3_DIR})
endif()

#------------------------------------------------------------------------------
# Conditionnaly include ExternalProject Target
#------------------------------------------------------------------------------

set(ep_common_args
  --no-warn-unused-cli
  -DMAKECOMMAND:STRING=${MAKECOMMAND}
  -DCMAKE_SKIP_RPATH:BOOL=ON
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DCMAKE_CXX_FLAGS_RELEASE:STRING=${CMAKE_CXX_FLAGS_RELEASE}
  -DCMAKE_CXX_FLAGS_DEBUG:STRING=${CMAKE_CXX_FLAGS_DEBUG}
  -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
  -DCMAKE_C_FLAGS_RELEASE:STRING=${CMAKE_C_FLAGS_RELEASE}
  -DCMAKE_C_FLAGS_DEBUG:STRING=${CMAKE_C_FLAGS_DEBUG}
  -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=${BUILD_TESTING}
  -DCMAKE_GENERATOR:STRING=${CMAKE_GENERATOR}
  -DCMAKE_EXTRA_GENERATOR:STRING=${CMAKE_EXTRA_GENERATOR}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
  -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
  -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  -DCMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH=${CMAKE_BUNDLE_OUTPUT_DIRECTORY}
  -DCTEST_NEW_FORMAT:BOOL=ON
  -DMEMORYCHECK_COMMAND_OPTIONS:STRING=${MEMORYCHECK_COMMAND_OPTIONS}
  -DMEMORYCHECK_COMMAND:PATH=${MEMORYCHECK_COMMAND}
  -DCMAKE_SHARED_LINKER_FLAGS:STRING=${CMAKE_SHARED_LINKER_FLAGS}
  -DCMAKE_EXE_LINKER_FLAGS:STRING=${CMAKE_EXE_LINKER_FLAGS}
  -DCMAKE_MODULE_LINKER_FLAGS:STRING=${CMAKE_MODULE_LINKER_FLAGS}
  -DSITE:STRING=${SITE}
  -DBUILDNAME:STRING=${BUILDNAME}
)

# Boost
option(USE_SYSTEM_BOOST "Use pre-installed version of the boost libraries" OFF)

set(git_protocol "git")

if(NOT USE_SYSTEM_BOOST)
  set(boost_version 1.41.0)
  set(boost_file
    "http://www.vtk.org/files/support/boost-${boost_version}.cmake-kitware.tar.gz")
  set(boost_md5 "f09997a2dad36627579b3e2215c25a48")
endif(NOT USE_SYSTEM_BOOST)

include(External_Boost)
include(External_CLAPACK)
include(External_ITKv4)
include(External_VTK)
include(External_BatchMake)

set(SlicerExecutionModel_DEPENDENCIES ITKv4)
include(External_SlicerExecutionModel)

set(spharmpdm_DEPENDENCIES CLAPACK Boost ITKv4 VTK SlicerExecutionModel)
#------------------------------------------------------------------------------
# Configure and build Slicer
#------------------------------------------------------------------------------
set(proj spharmpdm)

ExternalProject_Add(${proj}
  DEPENDS ${spharmpdm_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR spharmpdm-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${ep_common_args}
    -Dspharmpdm_SUPERBUILD:BOOL=OFF
    -DADDITIONAL_CXX_FLAGS:STRING=${ADDITIONAL_CXX_FLAGS}
    -DGIT_EXECUTABLE:FILEPATH=${GIT_EXECUTABLE}
    # ITK
    -DITK_DIR:PATH=${ITK_DIR}
    # VTK
    -DVTK_DIR:PATH=${VTK_DIR}
    # GenerateCLP_DIR
    -DGenerateCLP_DIR:PATH=${GenerateCLP_DIR}
    # Boost
    -DUSE_SYSTEM_BOOST:BOOL=OFF
    -DBoost_DIR:PATH=${Boost_DIR}
    # CLAPACK
    -DCLAPACK_DIR:PATH=${CLAPACK_DIR}
    ${trilinos_blas_args}
    ${SLICER_ARGS}
  INSTALL_COMMAND ""
  )
