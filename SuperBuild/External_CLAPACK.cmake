if( NOT EXTERNAL_SOURCE_DIRECTORY )
  set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
endif()
# Make sure this file is included only once
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()

set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

# Include dependent projects if any
set(extProjName CLAPACK)         #The find_package known name
set(proj        ${extProjName} ) #This local name

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory")
endif()


# Set dependency list
set(${proj}_DEPENDENCIES "")

#SlicerMacroCheckExternalProjectDependency(${proj})


if(NOT DEFINED ${extProjName}_DIR AND NOT ${USE_SYSTEM_${extProjName}})
  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

set(CLAPACK_version 3.2.1)
# We needed to patch CLAPACK to avoid building tests when BUILD_TESTING was set to OFF
# Instead of patching the source code each time, we copied the source code from 
# http://svn.slicer.org/Slicer3-lib-mirrors/trunk/clapack-3.2.1-CMAKE.tgz
# to https://github.com/NIRALUser/CLAPACK
set(${proj}_REPOSITORY ${git_protocol}://github.com/NIRALUser/CLAPACK.git)
set(${proj}_GIT_TAG b12609ee7addcad38d6b9265a7d1ae5446b05255)

# Turn off the warnings for CLAPACK on windows
string(REPLACE "/W3" "/W0" CMAKE_CXX_FLAGS_CLAPACK "${ep_common_cxx_flags}")
string(REPLACE "/W4" "/W0" CMAKE_CXX_FLAGS_CLAPACK "${CMAKE_CXX_FLAGS_CLAPACK}")
string(REPLACE "/W3" "/W0" CMAKE_C_FLAGS_CLAPACK "${ep_common_c_flags}")
string(REPLACE "/W4" "/W0" CMAKE_C_FLAGS_CLAPACK "${CMAKE_C_FLAGS_CLAPACK}")

#
# To fix compilation problem: relocation R_X86_64_32 against `a local symbol' can not be
# used when making a shared object; recompile with -fPIC
# See http://www.cmake.org/pipermail/cmake/2007-May/014350.html
#
if(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
  set(CMAKE_C_FLAGS_CLAPACK "-fPIC ${CMAKE_C_FLAGS_CLAPACK}")
endif()
  ### --- End Project specific additions
 ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${EXTERNAL_SOURCE_DIRECTORY}/${proj}
    BINARY_DIR ${proj}-build
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      # Compiler settings
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      # Options
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      -DBUILD_TESTING:BOOL=OFF
    #PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild/CLAPACK-PatchTesting.txt ${CMAKE_CURRENT_BINARY_DIR}/CLAPACK/CMakeLists.txt # No testing if BUILD_TESTING is set to OFF
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
    BUILD_COMMAND ${BUILD_COMMAND_STRING}
    )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
else()

  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} REQUIRED)
    if(NOT ${extProjName}_DIR)
      message(FATAL_ERROR "To use the system ${extProjName}, set ${extProjName}_DIR")
    endif()
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}v4, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)
