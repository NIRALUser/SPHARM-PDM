
#-----------------------------------------------------------------------------
set(verbose FALSE)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
enable_language(C)
enable_language(CXX)

#-----------------------------------------------------------------------------
enable_testing()
include(CTest)

set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

#-----------------------------------------------------------------------------
include(${CMAKE_CURRENT_SOURCE_DIR}/Common.cmake)
#-----------------------------------------------------------------------------
#If it is build as an extension
#-----------------------------------------------------------------------------
if( ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  if( APPLE )
    set( CMAKE_EXE_LINKER_FLAGS -Wl,-rpath,@loader_path/../../../../../ )
  endif()
endif()

#-----------------------------------------------------------------------------
# Enable and setup External project global properties
#-----------------------------------------------------------------------------
include(ExternalProject)
include(SlicerMacroEmptyExternalProject)
include(SlicerMacroCheckExternalProjectDependency)

# Compute -G arg for configuring external projects with the same CMake generator:
if(CMAKE_EXTRA_GENERATOR)
  set(gen "${CMAKE_EXTRA_GENERATOR} - ${CMAKE_GENERATOR}")
else()
  set(gen "${CMAKE_GENERATOR}")
endif()


# With CMake 2.8.9 or later, the UPDATE_COMMAND is required for updates to occur.
# For earlier versions, we nullify the update state to prevent updates and
# undesirable rebuild.
option(FORCE_EXTERNAL_BUILDS "Force rebuilding of external project (if they are updated)" OFF)
if(CMAKE_VERSION VERSION_LESS 2.8.9 OR NOT FORCE_EXTERNAL_BUILDS)
  set(cmakeversion_external_update UPDATE_COMMAND)
  set(cmakeversion_external_update_value "" )
else()
  set(cmakeversion_external_update LOG_UPDATE )
  set(cmakeversion_external_update_value 1)
endif()

#-----------------------------------------------------------------------------
# Superbuild option(s)
#-----------------------------------------------------------------------------
option(BUILD_STYLE_UTILS "Build uncrustify, cppcheck, & KWStyle" OFF)
CMAKE_DEPENDENT_OPTION(
  USE_SYSTEM_Uncrustify "Use system Uncrustify program" OFF
  "BUILD_STYLE_UTILS" OFF
  )
CMAKE_DEPENDENT_OPTION(
  USE_SYSTEM_KWStyle "Use system KWStyle program" OFF
  "BUILD_STYLE_UTILS" OFF
  )
CMAKE_DEPENDENT_OPTION(
  USE_SYSTEM_Cppcheck "Use system Cppcheck program" OFF
  "BUILD_STYLE_UTILS" OFF
  )

set(EXTERNAL_PROJECT_BUILD_TYPE "Release" CACHE STRING "Default build type for support libraries")

option(USE_SYSTEM_ITK "Build using an externally defined version of ITK" OFF)
option(USE_SYSTEM_SlicerExecutionModel "Build using an externally defined version of SlicerExecutionModel"  OFF)
option(USE_SYSTEM_VTK "Build using an externally defined version of VTK" OFF)
option(USE_SYSTEM_zlib "Build using external zlib" ON)
option(USE_SYSTEM_BatchMake "Build using an externally defined version of BatchMake" OFF)
option(USE_SYSTEM_CLAPACK "Build using an externally defined version of CLAPACK" OFF)
option(BUILD_SHARED_LIBS "Use shared libraries" OFF) #to give the user the option to configure their builds as they want

#------------------------------------------------------------------------------
# ${LOCAL_PROJECT_NAME} modules and external tools
#------------------------------------------------------------------------------

option(COMPILE_StatNonParamTestPDM "Compile StatNonParam and ShapeMancova" OFF)
option(COMPILE_shapeworks "Compile shapeworks." OFF)
option(COMPILE_ImageMath "Compile ImageMath." OFF)
option(COMPILE_MetaMeshTools "Compile MetaMeshTools." ON)
option(COMPILE_SegPostProcessCLP "Compile SegPostProcessCLP." ON)
option(COMPILE_GenParaMeshCLP "Compile GenParaMeshCLP." ON)
option(COMPILE_RadiusToMesh "Compile RadiusToMesh." ON)
option(COMPILE_ParaToSPHARMMeshCLP "Compile ParaToSPHARMMeshCLP." ON)
option(COMPILE_SpharmTool "Compile SpharmTool." ON)
option(COMPILE_ShapeAnalysisModule "Compile ShapeAnalysisModule." ON)
option(COMPILE_ParticleModule "Compile ParticleModule." OFF)

#------------------------------------------------------------------------------
# ${LOCAL_PROJECT_NAME} dependency list
#------------------------------------------------------------------------------

set(${LOCAL_PROJECT_NAME}_DEPENDENCIES ITKv4 SlicerExecutionModel VTK BatchMake CLAPACK)

if(BUILD_STYLE_UTILS)
  list(APPEND ${LOCAL_PROJECT_NAME}_DEPENDENCIES Cppcheck KWStyle Uncrustify)
endif()


#-----------------------------------------------------------------------------
# Define Superbuild global variables
#-----------------------------------------------------------------------------

# This variable will contain the list of CMake variable specific to each external project
# that should passed to ${CMAKE_PROJECT_NAME}.
# The item of this list should have the following form: <EP_VAR>:<TYPE>
# where '<EP_VAR>' is an external project variable and TYPE is either BOOL, STRING, PATH or FILEPATH.
# TODO Variable appended to this list will be automatically exported in ${LOCAL_PROJECT_NAME}Config.cmake,
# prefix '${LOCAL_PROJECT_NAME}_' will be prepended if it applies.
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS)

# The macro '_expand_external_project_vars' can be used to expand the list of <EP_VAR>.
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS) # List of CMake args to configure BRAINS
set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES) # List of CMake variable names

# Convenient macro allowing to expand the list of EP_VAR listed in ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
# The expanded arguments will be appended to the list ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS
# Similarly the name of the EP_VARs will be appended to the list ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES.
macro(_expand_external_project_vars)
  set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS "")
  set(${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES "")
  foreach(arg ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS})
    string(REPLACE ":" ";" varname_and_vartype ${arg})
    set(target_info_list ${target_info_list})
    list(GET varname_and_vartype 0 _varname)
    list(GET varname_and_vartype 1 _vartype)
    list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS -D${_varname}:${_vartype}=${${_varname}})
    list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES ${_varname})
  endforeach()
endmacro()

#-----------------------------------------------------------------------------
# Common external projects CMake variables
#-----------------------------------------------------------------------------
list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
  MAKECOMMAND:STRING
  CMAKE_SKIP_RPATH:BOOL
  CMAKE_BUILD_TYPE:STRING
  BUILD_SHARED_LIBS:BOOL
  CMAKE_CXX_COMPILER:PATH
  CMAKE_CXX_FLAGS_RELEASE:STRING
  CMAKE_CXX_FLAGS_DEBUG:STRING
  CMAKE_CXX_FLAGS:STRING
  CMAKE_C_COMPILER:PATH
  CMAKE_C_FLAGS_RELEASE:STRING
  CMAKE_C_FLAGS_DEBUG:STRING
  CMAKE_C_FLAGS:STRING
  CMAKE_SHARED_LINKER_FLAGS:STRING
  CMAKE_EXE_LINKER_FLAGS:STRING
  CMAKE_MODULE_LINKER_FLAGS:STRING
  CMAKE_GENERATOR:STRING
  CMAKE_EXTRA_GENERATOR:STRING
  CMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH
  CMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH
  CMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH
  CMAKE_BUNDLE_OUTPUT_DIRECTORY:PATH
  CTEST_NEW_FORMAT:BOOL
  MEMORYCHECK_COMMAND_OPTIONS:STRING
  MEMORYCHECK_COMMAND:PATH
  SITE:STRING
  BUILDNAME:STRING
  Subversion_SVN_EXECUTABLE:FILEPATH
  GIT_EXECUTABLE:FILEPATH
  )


_expand_external_project_vars()
set(COMMON_EXTERNAL_PROJECT_ARGS ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})
set(extProjName ${LOCAL_PROJECT_NAME})
set(proj        ${LOCAL_PROJECT_NAME})
SlicerMacroCheckExternalProjectDependency(${proj})

#-----------------------------------------------------------------------------
# Set CMake OSX variable to pass down the external project
#-----------------------------------------------------------------------------
set(CMAKE_OSX_EXTERNAL_PROJECT_ARGS)
if(APPLE)
  list(APPEND CMAKE_OSX_EXTERNAL_PROJECT_ARGS
    -DCMAKE_OSX_ARCHITECTURES=${CMAKE_OSX_ARCHITECTURES}
    -DCMAKE_OSX_SYSROOT=${CMAKE_OSX_SYSROOT}
    -DCMAKE_OSX_DEPLOYMENT_TARGET=${CMAKE_OSX_DEPLOYMENT_TARGET})
endif()

set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION  bin)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION  lib)
set(${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION  lib)

#-----------------------------------------------------------------------------
# Add external project CMake args
#-----------------------------------------------------------------------------
list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
  BUILD_EXAMPLES:BOOL
  BUILD_TESTING:BOOL
  LOCAL_PROJECT_NAME:STRING
  ITK_VERSION_MAJOR:STRING
  ITK_DIR:PATH
  VTK_DIR:PATH
  GenerateCLP_DIR:PATH
  SlicerExecutionModel_DIR:PATH
  CLAPACK_DIR:PATH
  BatchMake_DIR:PATH
  BOOST_INCLUDE_DIR:PATH
  BOOST_ROOT:PATH
  COMPILE_StatNonParamTestPDM:BOOL
  COMPILE_shapeworks:BOOL
  COMPILE_ImageMath:BOOL
  COMPILE_MetaMeshTools:BOOL
  COMPILE_SegPostProcessCLP:BOOL
  COMPILE_GenParaMeshCLP:BOOL
  COMPILE_RadiusToMesh:BOOL
  COMPILE_ParaToSPHARMMeshCLP:BOOL
  COMPILE_ShapeAnalysisModule:BOOL
  COMPILE_ParticleModule:BOOL
  COMPILE_SpharmTool:BOOL
  INSTALL_RUNTIME_DESTINATION:STRING
  INSTALL_LIBRARY_DESTINATION:STRING
  INSTALL_ARCHIVE_DESTINATION:STRING
  ${LOCAL_PROJECT_NAME}_CLI_INSTALL_LIBRARY_DESTINATION:PATH
  ${LOCAL_PROJECT_NAME}_CLI_INSTALL_ARCHIVE_DESTINATION:PATH
  ${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION:PATH
  )

if( SPHARM-PDM_BUILD_SLICER_EXTENSION )
  list(APPEND ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARS
    MIDAS_PACKAGE_API_KEY:STRING
    MIDAS_PACKAGE_EMAIL:STRING
    MIDAS_PACKAGE_URL:STRING
    Slicer_DIR:PATH
    SPHARM-PDM_BUILD_SLICER_EXTENSION:BOOL
    )
endif()
_expand_external_project_vars()
set(COMMON_EXTERNAL_PROJECT_ARGS ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})

if(verbose)
  message("Inner external project args:")
  foreach(arg ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_ARGS})
    message("  ${arg}")
  endforeach()
endif()

string(REPLACE ";" "^" ${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES "${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES}")

if(verbose)
  message("Inner external project argnames:")
  foreach(argname ${${CMAKE_PROJECT_NAME}_SUPERBUILD_EP_VARNAMES})
    message("  ${argname}")
  endforeach()
endif()


#------------------------------------------------------------------------------
# Configure and build
#------------------------------------------------------------------------------
set(proj ${LOCAL_PROJECT_NAME})
ExternalProject_Add(${proj}
  DEPENDS ${${LOCAL_PROJECT_NAME}_DEPENDENCIES}
  DOWNLOAD_COMMAND ""
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  BINARY_DIR ${LOCAL_PROJECT_NAME}-build
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
    --no-warn-unused-cli # HACK Only expected variables should be passed down.
    ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
    ${COMMON_EXTERNAL_PROJECT_ARGS}
    -D${LOCAL_PROJECT_NAME}_SUPERBUILD:BOOL=OFF
    -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/${proj}-install
  )

## Force rebuilding of the main subproject every time building from super structure
ExternalProject_Add_Step(${proj} forcebuild
    COMMAND ${CMAKE_COMMAND} -E remove
    ${CMAKE_CURRENT_BUILD_DIR}/${proj}-prefix/src/${proj}-stamp/${proj}-build
    DEPENDEES configure
    DEPENDERS build
    ALWAYS 1
  )

if( ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  unsetForSlicer( NAMES VTK_DIR SlicerExecutionModel_DIR DCMTK_DIR ITK_DIR CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS CMAKE_C_FLAGS)
  # Create fake imported target to avoid importing Slicer target: See SlicerConfig.cmake:line 820
  add_library(SlicerBaseLogic SHARED IMPORTED)
  #find_package(Slicer REQUIRED)
  #include(${Slicer_USE_FILE})
  set( EXTENSION_NO_CLI asc2meta asc2vtk BYU2VTK MeshMath Meta2STL Meta2VTK SpharmTool STL2Meta VTK2Meta )
  set( EXTENSION_CLI GenParaMeshCLP ParaToSPHARMMeshCLP RadiusToMesh SegPostProcessCLP ShapeAnalysisModule  )
  set( EXTENSION_TESTS ${EXTENSION_NO_CLI} ${EXTENSION_CLI} GenParaMeshCLPTest SegPostProcessCLPTest )
  foreach( VAR ${EXTENSION_TESTS})
    add_executable(${VAR} IMPORTED)
    set_property(TARGET ${VAR} PROPERTY IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}/${VAR}${fileextension})
  endforeach()
  if(BUILD_TESTING)
    add_subdirectory(Testing)
  endif()
  resetForSlicer( NAMES VTK_DIR ITK_DIR SlicerExecutionModel_DIR CMAKE_C_COMPILER CMAKE_CXX_COMPILER CMAKE_CXX_FLAGS CMAKE_C_FLAGS)
  set( NOCLI_INSTALL_DIR ${SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION}/../ExternalBin)
  foreach( VAR ${EXTENSION_NO_CLI})
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}/${VAR}${fileextension} DESTINATION ${NOCLI_INSTALL_DIR} )
  endforeach()
  foreach( VAR ${EXTENSION_CLI})
    install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${proj}-install/${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}/${VAR}${fileextension} DESTINATION ${SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION} )
  endforeach()
  install( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bmm DESTINATION ${SlicerExecutionModel_DEFAULT_CLI_INSTALL_RUNTIME_DESTINATION}/.. )
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
  #include(${Slicer_EXTENSION_CPACK})
endif()
