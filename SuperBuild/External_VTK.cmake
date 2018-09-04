
set(proj VTK)

# Set dependency list
set(${proj}_DEPENDS
  ""
  )

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(Slicer_USE_SYSTEM_${proj})
  unset(VTK_DIR CACHE)
  find_package(VTK 8.1.1 REQUIRED)
endif()

# Sanity checks
if(DEFINED VTK_DIR AND NOT EXISTS ${VTK_DIR})
  message(FATAL_ERROR "VTK_DIR [${VTK_DIR}] variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED VTK_DIR AND NOT Slicer_USE_SYSTEM_${proj})

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/kitware/VTK.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_TAG
    "v8.1.1"
    QUIET
    )

  set(EP_VERSION_ID "vtk-8.1")

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(EP_INSTALL_DIR ${CMAKE_BINARY_DIR}/${proj}-install)

  set(log_command_output "0")
  if(NOT "$ENV{DASHBOARD_TEST_FROM_CTEST}" STREQUAL "")
    set(log_command_output "1")
  endif()

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${Slicer_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${Slicer_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    CMAKE_CACHE_ARGS
      # Compiler settings
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      # Output directories
      ## NA
      # Install directories
      -DCMAKE_INSTALL_PREFIX:PATH=${EP_INSTALL_DIR}
      # Options
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DVTK_USE_PARALLEL:BOOL=ON
      -DVTK_DEBUG_LEAKS:BOOL=${${PROJECT_NAME}_USE_VTK_DEBUG_LEAKS}
      -DVTK_LEGACY_REMOVE:BOOL=ON
      -DVTK_WRAP_TCL:BOOL=OFF
      -DVTK_WRAP_PYTHON:BOOL=OFF
      -DVTK_USE_GUISUPPORT:BOOL=OFF
      -DVTK_USE_QT:BOOL=OFF
    # Install
    ## We really do want to install in order to limit # of include paths INSTALL_COMMAND ""
    # Wrap commands to ignore log output from dashboards
    LOG_CONFIGURE ${log_command_output}
    LOG_BUILD     ${log_command_output}
    LOG_INSTALL   1
    DEPENDS
      ${${proj}_DEPENDS}
    )
  set(VTK_DIR ${EP_INSTALL_DIR}/lib/cmake/${EP_VERSION_ID})
else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(VTK_DIR:PATH)
