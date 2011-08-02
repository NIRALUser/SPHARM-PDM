# The Boost external project for Titan

set(boost_source  "${CMAKE_CURRENT_BINARY_DIR}/Boost")
set(boost_binary  "${CMAKE_CURRENT_BINARY_DIR}/Boost-Build")
set(Boost_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/Boost/boost")

if(MSVC)
  set(boost_lib_args
    -DENABLE_SHARED:BOOL=OFF
    -DENABLE_STATIC:BOOL=ON
  )
else()
  set(boost_lib_args
    -DENABLE_SHARED:BOOL=ON
    -DENABLE_STATIC:BOOL=OFF
  )
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
    -DBUILD_SHARED_LIBS:BOOL=ON
    -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG:PATH=${boost_binary}/lib
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_DEBUG:PATH=${boost_binary}/lib
    -DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE:PATH=${boost_binary}/lib
    -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY_RELEASE:PATH=${boost_binary}/lib
    -DBUILD_EXAMPLES:BOOL=OFF
    -DBUILD_TESTING:BOOL=OFF
    -DBUILD_VERSIONED:BOOL=OFF
    -DINSTALL_VERSIONED:BOOL=OFF
    -DWITH_MPI:BOOL=OFF
    -DWITH_PYTHON:BOOL=OFF
    ${Boost_EXTRA_ARGS}
  INSTALL_COMMAND ""
  )

# These variables are used to find Boost by other projects
set(Boost_INCLUDE_DIR "${boost_source}" CACHE PATH "" FORCE)
set(BOOST_LIBRARYDIR "${boost_binary}/lib" CACHE PATH "" FORCE)
