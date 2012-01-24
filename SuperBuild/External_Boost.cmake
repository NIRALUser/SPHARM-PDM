set(boost_source  "${CMAKE_CURRENT_BINARY_DIR}/Boost")
set(boost_binary  "${CMAKE_CURRENT_BINARY_DIR}/Boost-Build")

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

ExternalProject_Add(Boost
  DOWNLOAD_DIR ${CMAKE_CURRENT_BINARY_DIR}
  URL ${boost_file}
  URL_MD5 ${boost_md5}
  UPDATE_COMMAND ""
  SOURCE_DIR ${boost_source}
  BINARY_DIR ${boost_binary}
  CMAKE_ARGS
  ${boost_lib_args}
  ${titan_compiler_args}
  -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_CURRENT_BINARY_DIR}/Boost-install
  -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
  -DBUILD_EXAMPLES:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DBUILD_VERSIONED:BOOL=OFF
  -DINSTALL_VERSIONED:BOOL=OFF
  -DWITH_MPI:BOOL=OFF
  -DWITH_PYTHON:BOOL=OFF
  ${Boost_EXTRA_ARGS}
  )

# These variables are used to find Boost by other projects
set(Boost_DIR ${CMAKE_CURRENT_BINARY_DIR}/Boost-install)
