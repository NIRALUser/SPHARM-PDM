# Download FindVcvars.cmake
set(dest_file "${CMAKE_CURRENT_BINARY_DIR}/FindVcvars.cmake")
set(expected_hash "9072eac4ca2d7a06c6a69cefc315338d322954184a7410892e9afdb2486d9fb7")
set(url "https://raw.githubusercontent.com/scikit-build/cmake-FindVcvars/v1.0/FindVcvars.cmake")
if(NOT EXISTS ${dest_file})
  file(DOWNLOAD ${url} ${dest_file} EXPECTED_HASH SHA256=${expected_hash})
else()
  file(SHA256 ${dest_file} current_hash)
  if(NOT ${current_hash} STREQUAL ${expected_hash})
    file(DOWNLOAD ${url} ${dest_file} EXPECTED_HASH SHA256=${expected_hash})
  endif()
endif()
