
set(proj CLAPACK)

# Set dependency list
set(${proj}_DEPENDS
  ""
  )

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(Slicer_USE_SYSTEM_${proj})
  message(FATAL_ERROR "Enabling Slicer_USE_SYSTEM_${proj} is not supported !")
endif()

# Sanity checks
if(DEFINED CLAPACK_DIR AND NOT EXISTS ${CLAPACK_DIR})
  message(FATAL_ERROR "CLAPACK_DIR [${CLAPACK_DIR}] variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED ${proj}_DIR AND NOT Slicer_USE_SYSTEM_${proj})

  # We needed to patch CLAPACK to avoid building tests when BUILD_TESTING was set to OFF
  # Instead of patching the source code each time, we copied the source code from
  # http://svn.slicer.org/Slicer3-lib-mirrors/trunk/clapack-3.2.1-CMAKE.tgz
  # to https://github.com/NIRALUser/CLAPACK

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/NIRALUser/CLAPACK.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_TAG
    "b12609ee7addcad38d6b9265a7d1ae5446b05255"
    QUIET
    )

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  # Turn off the warnings for CLAPACK on windows
  string(REPLACE "/W3" "/W0" CMAKE_C_FLAGS_CLAPACK "${CMAKE_C_FLAGS}")
  string(REPLACE "/W4" "/W0" CMAKE_C_FLAGS_CLAPACK "${CMAKE_C_FLAGS_CLAPACK}")

  ExternalProject_Add(${proj}
    ${${proj}_EP_ARGS}
    GIT_REPOSITORY "${Slicer_${proj}_GIT_REPOSITORY}"
    GIT_TAG "${Slicer_${proj}_GIT_TAG}"
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    CMAKE_CACHE_ARGS
      # Compiler settings
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS_CLAPACK}
      # Options
      -DBUILD_TESTING:BOOL=OFF
    INSTALL_COMMAND ""
    # Wrap commands to ignore log output from dashboards
    LOG_CONFIGURE 1
    LOG_BUILD     1
    DEPENDS
      ${${proj}_DEPENDS}
    )
  set(CLAPACK_DIR ${EP_BINARY_DIR})

else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(CLAPACK_DIR:PATH)
