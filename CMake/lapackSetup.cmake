# sets the following variables
# LAPACK_FOUND         - Set to true when lapack is found.
# LAPACK_INCLUDE_DIR   - Include directory for lapack headers.
# LAPACK_LIBRARY_DIR   - Link directory for lapack libraries.
# LAPACK_LIBRARIES     - list of lapack libraries

FIND_PATH(LAPACK_LIBRARY_DIR liblapack.a
  ${NeuroLib_LIBS_PATH}
  /usr/lib /opt/local/lib /usr/local/lib /tools/lib  /usr/lib64
  DOC " directory for lapack, where liblapack can be found"
  )
 
IF(LAPACK_LIBRARY_DIR)
  SET(LAPACK_FOUND 1)
  SET(LAPACK_INCLUDE_DIR ${LAPACK_LIBRARY_DIR}/../include)
   
  IF(WIN32)
    SET(LAPACK_LIBRARIES
      lapack
      blas
      )
  ENDIF(WIN32)
  IF(CMAKE_SYSTEM MATCHES "Darwin*")
    SET(LAPACK_LIBRARIES
      lapack
      blas
      )
  ENDIF(CMAKE_SYSTEM MATCHES "Darwin*")
  IF(CMAKE_SYSTEM MATCHES "Sun*")
    SET(LAPACK_LIBRARIES
      lapack
      blas
      I77
      F77
      )
  ENDIF(CMAKE_SYSTEM MATCHES "Sun*")
  IF(CMAKE_SYSTEM MATCHES "Linux*")
	EXEC_PROGRAM(gcc ARGS --version OUTPUT_VARIABLE CMAKE_C_COMPILER_VERSION)
	IF(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
		SET(GCC_FORTRAN_LIBRARY gfortran)
	ELSE(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
		SET(GCC_FORTRAN_LIBRARY g2c)
	ENDIF(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
    SET(LAPACK_LIBRARIES
      lapack
      blas
      ${GCC_FORTRAN_LIBRARY}
      )
  ENDIF(CMAKE_SYSTEM MATCHES "Linux*")
   
ELSE(LAPACK_LIBRARY_DIR)
  MESSAGE(FATAL_ERROR "LAPACK library not found!\n" "Please go to http://www.ia.unc.edu/dev/tutorials/InstallLib")
  SET(LAPACK_INCLUDE_DIR NOT_FOUND)
  SET(LAPACK_FOUND 0)
ENDIF(LAPACK_LIBRARY_DIR)
 
INCLUDE_DIRECTORIES(${LAPACK_INCLUDE_DIR})
LINK_DIRECTORIES(${LAPACK_LIBRARY_DIR})
LINK_LIBRARIES(${LAPACK_LIBRARIES})
