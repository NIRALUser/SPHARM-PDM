# sets the following variables
# BOOST_FOUND         - Set to true when boost is found.
# BOOST_INCLUDE_DIR   - Include directory for boost headers.
# BOOST_LIBRARY_DIR   - Link directory for boost libraries.
# BOOST_LIBRARIES     - list of boost libraries

FIND_PATH(BOOST_INCLUDE_DIR boost/config.hpp
  /opt/local/Boost/boost_1_39_0
  /usr/include
  /usr/local/include
)
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})


FIND_PATH(BOOST_LIBRARY_DIR libboost_filesystem.a
  /opt/local/Boost/boost_1_39_0_linux64/lib
  /usr/lib64 /usr/lib /opt/local/lib /usr/local/lib /tools/lib  
  DOC " directory for boost, where libboost can be found"
)

IF(BOOST_LIBRARY_DIR)
  SET(BOOST_FOUND 1)
ELSE(BOOST_LIBRARY_DIR)
  MESSAGE(FATAL_ERROR "BOOST library not found!\n")
  SET(BOOST_INCLUDE_DIR NOT_FOUND)
  SET(BOOST_FOUND 0)
ENDIF(BOOST_LIBRARY_DIR)

# Win32 has automatic linking of boost libraries so only add boost
# libaries on unix
IF(NOT WIN32 OR CYGWIN)
  SET(BOOST_LIBRARIES boost_program_options  boost_system boost_filesystem)
ENDIF(NOT WIN32 OR CYGWIN)

INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
LINK_DIRECTORIES(${BOOST_LIBRARY_DIR})
LINK_LIBRARIES(${BOOST_LIBRARIES})
