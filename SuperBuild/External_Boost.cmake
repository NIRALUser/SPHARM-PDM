# Boost
set(boost_version 1.41.0)
set(boost_file
  "http://www.vtk.org/files/support/boost-${boost_version}.cmake-kitware.tar.gz")
set(boost_md5 "f09997a2dad36627579b3e2215c25a48")

set(boost_source  "${CMAKE_CURRENT_BINARY_DIR}/Boost")
set(boost_binary  "${CMAKE_CURRENT_BINARY_DIR}/Boost-Build")
set(boost_install "${CMAKE_CURRENT_BINARY_DIR}/Boost-install")

if(MSVC)
  set(boost_lib_args
    -DENABLE_SHARED:BOOL=OFF
    -DENABLE_STATIC:BOOL=ON
  )
else()
  if(${BUILD_SHARED_LIBS})
    set(boost_lib_args
      -DENABLE_SHARED:BOOL=ON
      -DENABLE_STATIC:BOOL=OFF
      )
  else()
    set(boost_lib_args
      -DENABLE_SHARED:BOOL=OFF
      -DENABLE_STATIC:BOOL=ON
      )
  endif()
endif()

  ### --- Project specific additions here
  set(${proj}_CMAKE_OPTIONS
      ${boost_lib_args}
      ${titan_compiler_args}
      -DBUILD_VERSIONED:BOOL=OFF
      -DINSTALL_VERSIONED:BOOL=OFF
      -DWITH_MPI:BOOL=OFF
      -DWITH_PYTHON:BOOL=OFF
      ${Boost_EXTRA_ARGS}
    )
  ### --- End Project specific additions
  ExternalProject_Add(Boost
    DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
    URL ${boost_file}
    URL_MD5 ${boost_md5}
    SOURCE_DIR ${boost_source}
    BINARY_DIR ${boost_binary}
    UPDATE_COMMAND ""
    CMAKE_GENERATOR ${gen}
    CMAKE_ARGS
      ${CMAKE_OSX_EXTERNAL_PROJECT_ARGS}
      ${COMMON_EXTERNAL_PROJECT_ARGS}
      ## HACK: NEED TO INSTALL BOOST SO THAT WE CAN FIND IT LATER
      ## INSTALL ""
      -DCMAKE_INSTALL_PREFIX:PATH=${boost_install}
      -DBUILD_EXAMPLES:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      ${${proj}_CMAKE_OPTIONS}
    DEPENDS
      ${${proj}_DEPENDENCIES}
    BUILD_COMMAND ${BUILD_COMMAND_STRING}
    )

# These variables are used to find Boost by other projects
set(Boost_DIR ${boost_install})
