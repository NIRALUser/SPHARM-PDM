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
set(extProjName BatchMake) #The find_package known name
set(proj      ${extProjName}) #This local name
set(${extProjName}_REQUIRED_VERSION ${${extProjName}_VERSION_MAJOR})  #If a required version is necessary, then set this, else leave blank


# Sanity checks
#if(DEFINED ${extProjName}_DIR AND NOT EXISTS ${${extProjName}_DIR})
#  message(FATAL_ERROR "${extProjName}_DIR variable is defined but corresponds to non-existing directory (${${extProjName}_DIR})")
#endif()
set(${proj}_DEPENDENCIES ITKv4)


if(NOT ( DEFINED "USE_SYSTEM_${extProjName}" AND "${USE_SYSTEM_${extProjName}}" ) )
  # Include dependent projects if any
  SlicerMacroCheckExternalProjectDependency(${proj})

ExternalProject_Add(${proj}
  GIT_REPOSITORY ${git_protocol}://github.com/NIRALUser/BatchMake.git
  GIT_TAG "5e1c185b564ac2acef35ba68fbde370c19ba33de" 
  SOURCE_DIR BatchMake
  BINARY_DIR BatchMake-build
  ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
    ${COMMON_EXTERNAL_PROJECT_ARGS}
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
    -DUSE_FLTK:BOOL=OFF
    -DGRID_SUPPORT:BOOL=OFF
    -DUSE_SPLASHSCREEN:BOOL=OFF
    -DITK_DIR:PATH=${ITK_DIR}
    -DBUILD_SHARED_LIBS:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
    -DDASHBOARD_SUPPORT:BOOL=OFF
    ${BatchMakeCURLCmakeArg}
  DEPENDS ${${proj}_DEPENDENCIES}
  )
  set(${extProjName}_DIR ${CMAKE_BINARY_DIR}/${proj}-install/lib/BatchMake)
  #message(FATAL_ERROR "${extProjName}_DIR:" ${${extProjName}_DIR})
else()
  if(${USE_SYSTEM_${extProjName}})
    find_package(${extProjName} REQUIRED)
    message("USING the system ${extProjName}, set ${extProjName}_DIR=${${extProjName}_DIR}")
  endif()
  # The project is provided using ${extProjName}_DIR, nevertheless since other
  # project may depend on ${extProjName}, let's add an 'empty' one
  SlicerMacroEmptyExternalProject(${proj} "${${proj}_DEPENDENCIES}")
endif()

list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS ${extProjName}_DIR:PATH)

ProjectDependancyPop(CACHED_extProjName extProjName)
ProjectDependancyPop(CACHED_proj proj)

