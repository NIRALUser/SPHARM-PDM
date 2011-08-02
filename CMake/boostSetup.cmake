# sets the following variables
# BOOST_FOUND         - Set to true when boost is found.
# BOOST_INCLUDE_DIR   - Include directory for boost headers.
# BOOST_LIBRARY_DIR   - Link directory for boost libraries.
# BOOST_LIBRARIES     - list of boost libraries

find_path(BOOST_INCLUDE_DIR boost/config.hpp
  /opt/local/Boost/boost_1_39_0
  /usr/include
  /usr/local/include
)
include_directories(${BOOST_INCLUDE_DIR})


find_path(BOOST_LIBRARY_DIR libboost_filesystem.a
  /opt/local/Boost/boost_1_39_0_linux64/lib
  /usr/lib64 /usr/lib /opt/local/lib /usr/local/lib /tools/lib
  DOC " directory for boost, where libboost can be found"
)

if(BOOST_LIBRARY_DIR)
  set(BOOST_FOUND 1)
else(BOOST_LIBRARY_DIR)
  message(FATAL_ERROR "BOOST library not found!\n")
  set(BOOST_INCLUDE_DIR NOT_FOUND)
  set(BOOST_FOUND 0)
endif(BOOST_LIBRARY_DIR)

# Win32 has automatic linking of boost libraries so only add boost
# libaries on unix
if(NOT WIN32 OR CYGWIN)
  set(BOOST_LIBRARIES boost_program_options  boost_system boost_filesystem)
endif(NOT WIN32 OR CYGWIN)

include_directories(${BOOST_INCLUDE_DIR})
link_directories(${BOOST_LIBRARY_DIR})
link_libraries(${BOOST_LIBRARIES})
