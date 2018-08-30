if( NOT EXTERNAL_SOURCE_DIRECTORY )
  set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
endif()

# Make sure this file is included only once by creating globally unique varibles
# based on the name of this included file.
get_filename_component(CMAKE_CURRENT_LIST_FILENAME ${CMAKE_CURRENT_LIST_FILE} NAME_WE)
if(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED)
  return()
endif()
set(${CMAKE_CURRENT_LIST_FILENAME}_FILE_INCLUDED 1)

## External_${extProjName}.cmake files can be recurisvely included,
## and cmake variables are global, so when including sub projects it
## is important make the extProjName and proj variables
## appear to stay constant in one of these files.
## Store global variables before overwriting (then restore at end of this file.)
ProjectDependancyPush(CACHED_extProjName ${extProjName})
ProjectDependancyPush(CACHED_proj ${proj})

# Make sure that the ExtProjName/IntProjName variables are unique globally
# even if other External_${ExtProjName}.cmake files are sourced by
# SlicerMacroCheckExternalProjectDependency
set(extProjName ITK) #The find_package known name
set(proj      ITKv4) #This local name
set(${extProjName}_REQUIRED_VERSION ${${extProjName}_VERSION_MAJOR})  #If a required version is necessary, then set this, else leave blank

#if(${USE_SYSTEM_${extProjName}})
#  unset(${extProjName}_DIR CACHE)
#endif()

# Sanity checks
#if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
#  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
#endif()
set(${proj}_DEPENDENCIES "")

if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  #message(STATUS "${__indent}Adding project ${proj}")
  # Set dependency list

  # Set CMake OSX variable to pass down the external project
  set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
  if(APPLE)
    list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
      -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
      -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
      -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
  endif()

  # Include dependent projects if any
  SlicerMacroCheckExternalProjectDependency(${proj})

  # HACK This code fixes a loony problem with HDF5 -- it doesn't
  #      link properly if -fopenmp is used.
  string(REPLACE "-fopenmp" "" ITK_CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "-fopenmp" "" ITK_CMAKE_CXX_FLAGS "${CMAKE_CX_FLAGS}")

  set(${proj}_CMAKE_OPTIONS
      # Compiler settings
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      # Options
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
      -DITK_LEGACY_REMOVE:BOOL=ON
      -DITKV3_COMPATIBILITY:BOOL=OFF
      -DITK_BUILD_DEFAULT_MODULES:BOOL=ON
      -DModule_ITKReview:BOOL=OFF
      #-DITK_INSTALL_NO_DEVELOPMENT:BOOL=ON
      -DKWSYS_USE_MD5:BOOL=ON # Required by SlicerExecutionModel
      -DITK_WRAPPING:BOOL=OFF #${BUILD_SHARED_LIBS} ## HACK:  QUICK CHANGE
      -DModule_MGHIO:BOOL=ON  # Allow building of the MGHIO classes
    )
  ### --- End Project specific additions
  set(${proj}_REPOSITORY ${git_protocol}://itk.org/ITK.git)
  set(${proj}_GIT_TAG v4.13.1)
  set(ITK_VERSION_ID ITK-4.13)

  ExternalProject_Add(${proj}
    GIT_REPOSITORY ${${proj}_REPOSITORY}
    GIT_TAG ${${proj}_GIT_TAG}
    SOURCE_DIR ${EXTERNAL_SOURCE_DIRECTORY}/${proj}
    BINARY_DIR ${proj}-build
    LOG_CONFIGURE 0  # Wrap configure in script to ignore log output from dashboards
    LOG_BUILD     0  # Wrap build in script to to ignore log output from dashboards
    LOG_TEST      0  # Wrap test in script to to ignore log output from dashboards
    LOG_INSTALL   0  # Wrap install in script to to ignore log output from dashboards
    CMAKE_GENERATOR ${gen}
    CMAKE_CACHE_ARGS
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      ${${proj}_CMAKE_OPTIONS}
## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDENCIES}
  )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/cmake/${ITK_VERSION_ID})
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} ${ITK_VERSION_MAJOR} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)
