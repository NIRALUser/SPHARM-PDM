
set(proj LAPACK)

# Set dependency list
set(${proj}_DEPENDS
  ""
  )

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})
  unset(LAPACK_DIR CACHE)
  find_package(LAPACK REQUIRED)
  unset(LAPACKE_DIR CACHE)
  find_package(LAPACKE REQUIRED)
endif()

# Sanity checks
if(DEFINED LAPACK_DIR AND NOT EXISTS ${LAPACK_DIR})
  message(FATAL_ERROR "LAPACK_DIR [${LAPACK_DIR}] variable is defined but corresponds to nonexistent directory")
endif()
if(DEFINED LAPACKE_DIR AND NOT EXISTS ${LAPACKE_DIR})
  message(FATAL_ERROR "LAPACKE_DIR [${LAPACKE_DIR}] variable is defined but corresponds to nonexistent directory")
endif()

if(NOT DEFINED LAPACK_DIR AND NOT Slicer_USE_SYSTEM_LAPACK
   AND NOT DEFINED LAPACKE_DIR AND NOT Slicer_USE_SYSTEM_LAPACKE)

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_REPOSITORY
    "${EP_GIT_PROTOCOL}://github.com/Reference-LAPACK/lapack.git"
    QUIET
    )

  ExternalProject_SetIfNotDefined(
    Slicer_${proj}_GIT_TAG
    "eb6fac41396e383e5c780b1afd01512848a7bcf9"
    QUIET
    )

  # Default values
  if(DEFINED Slicer_DIR OR DEFINED Slicer_SOURCE_DIR)
    set(LAPACK_INSTALL_RUNTIME_DIR ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
    set(LAPACK_INSTALL_LIBRARY_DIR ${Slicer_INSTALL_THIRDPARTY_LIB_DIR})
  else()
    if(NOT DEFINED LAPACK_INSTALL_RUNTIME_DIR)
      set(LAPACK_INSTALL_RUNTIME_DIR "bin")
    endif()
    if(NOT DEFINED LAPACK_INSTALL_LIBRARY_DIR)
      set(LAPACK_INSTALL_LIBRARY_DIR "lib")
    endif()
  endif()

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  if(MSVC)
    # Generator
    find_program(NINJA_EXECUTABLE ninja)
    if(NINJA_EXECUTABLE)
      set(generator "Ninja")
    else()
      set(generator "NMake Makefiles")
    endif()

    # flang
    set(Fortran_COMPILER_ID "Flang")
    find_package(Fortran REQUIRED)
    mark_as_superbuild(Fortran_COMPILER_ID)
    mark_as_superbuild(Fortran_${Fortran_COMPILER_ID}_EXECUTABLE)

    set(build_type "Release")

    set(cmake_cache_args
      # Compiler settings
      -DCMAKE_C_COMPILER:FILEPATH=${Fortran_Flang_CLANG_CL_EXECUTABLE}
      -DCMAKE_CXX_COMPILER:FILEPATH=${Fortran_Flang_CLANG_CL_EXECUTABLE}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
      -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
      -DCMAKE_Fortran_COMPILER:FILEPATH=${Fortran_Flang_EXECUTABLE}
      # Output directories
      ## NA
      # Install directories
      -DLAPACK_INSTALL_RUNTIME_DIR:STRING=${LAPACK_INSTALL_RUNTIME_DIR}
      -DLAPACK_INSTALL_LIBRARY_DIR:STRING=${LAPACK_INSTALL_LIBRARY_DIR}
      # Options
      -DBUILD_TESTING:BOOL=OFF
      -DCBLAS:BOOL=OFF
      -DLAPACKE:BOOL=ON
      -DCMAKE_BUILD_TYPE:STRING=${build_type}
      )
    if(NINJA_EXECUTABLE)
      list(APPEND cmake_cache_args
        -DCMAKE_MAKE_PROGRAM:FILEPATH=${NINJA_EXECUTABLE}
        )
    endif()

    # Configure command
    find_package(Vcvars REQUIRED)
    set(EXTERNAL_PROJECT_CONFIGURE_COMMAND CONFIGURE_COMMAND
      ${CMAKE_COMMAND} -E env
        CC=${Fortran_Flang_CLANG_CL_EXECUTABLE}
        CXX=${Fortran_Flang_CLANG_CL_EXECUTABLE}
        LIB=${Fortran_${Fortran_COMPILER_ID}_IMPLICIT_LINK_DIRECTORIES}
      ${Vcvars_WRAPPER_BATCH_FILE} ${CMAKE_COMMAND} -G ${generator}
        ${cmake_cache_args}
        ${EP_SOURCE_DIR}
      )

    # Build command
    set(EXTERNAL_PROJECT_BUILD_COMMAND BUILD_COMMAND
      ${CMAKE_COMMAND} --build "." --config ${build_type}
      )

    ExternalProject_Add(${proj}
      USES_TERMINAL_DOWNLOAD 1
      USES_TERMINAL_UPDATE 1
      USES_TERMINAL_CONFIGURE 1
      USES_TERMINAL_BUILD 1
      USES_TERMINAL_TEST 1
      USES_TERMINAL_INSTALL 1
      GIT_REPOSITORY "${Slicer_${proj}_GIT_REPOSITORY}"
      GIT_TAG "${Slicer_${proj}_GIT_TAG}"
      SOURCE_DIR ${EP_SOURCE_DIR}
      BINARY_DIR ${EP_BINARY_DIR}
      ${EXTERNAL_PROJECT_CONFIGURE_COMMAND}
      ${EXTERNAL_PROJECT_BUILD_COMMAND}
      INSTALL_COMMAND ""
      DEPENDS
        ${${proj}_DEPENDS}
      )

    #-----------------------------------------------------------------------------
    # Launcher setting specific to build tree
    set(${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
      ${Fortran_${Fortran_COMPILER_ID}_RUNTIME_DIRECTORIES}
      )
    mark_as_superbuild(
      VARS ${proj}_LIBRARY_PATHS_LAUNCHER_BUILD
      LABELS "LIBRARY_PATHS_LAUNCHER_BUILD"
      )

  else()

    find_package(Fortran REQUIRED)
    mark_as_superbuild(Fortran_COMPILER_ID:STRING)
    mark_as_superbuild(Fortran_${Fortran_COMPILER_ID}_EXECUTABLE)

    set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS)
    if(APPLE)
      # Since Detect/Verify tests associated with FortranCInterface CMake module
      # fails when at least the CMAKE_OSX_ARCHITECTURES option is set to x86_64,
      # we explicitly disable the CMAKE_OSX_* options here.
      # Failure was observed when using gfortran_osx-64 from conda.
      list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
        -DCMAKE_OSX_ARCHITECTURES:STRING=
        -DCMAKE_OSX_SYSROOT:PATH=
        -DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=
        )
    endif()

    ExternalProject_Add(${proj}
      ${${proj}_EP_ARGS}
      GIT_REPOSITORY "${Slicer_${proj}_GIT_REPOSITORY}"
      GIT_TAG "${Slicer_${proj}_GIT_TAG}"
      SOURCE_DIR ${EP_SOURCE_DIR}
      BINARY_DIR ${EP_BINARY_DIR}
      CMAKE_CACHE_ARGS
        ${EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS}
        # Compiler settings
        -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
        -DCMAKE_CXX_STANDARD_REQUIRED:BOOL=${CMAKE_CXX_STANDARD_REQUIRED}
        -DCMAKE_CXX_EXTENSIONS:BOOL=${CMAKE_CXX_EXTENSIONS}
        -DCMAKE_Fortran_COMPILER:FILEPATH=${Fortran_${Fortran_COMPILER_ID}_EXECUTABLE}
        -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=1
        # Output directories
        ## NA
        # Install directories
        -DLAPACK_INSTALL_RUNTIME_DIR:STRING=${LAPACK_INSTALL_RUNTIME_DIR}
        -DLAPACK_INSTALL_LIBRARY_DIR:STRING=${LAPACK_INSTALL_LIBRARY_DIR}
        # Options
        -DBUILD_TESTING:BOOL=OFF
        -DCBLAS:BOOL=OFF
        -DLAPACKE:BOOL=ON
      INSTALL_COMMAND ""
      DEPENDS
        ${${proj}_DEPENDS}
      )
    endif()

  set(LAPACK_DIR ${EP_BINARY_DIR})
  set(LAPACKE_DIR ${EP_BINARY_DIR})
else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(LAPACK_DIR:PATH)
mark_as_superbuild(LAPACKE_DIR:PATH)
