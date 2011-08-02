# sets the following variables
# LAPACK_FOUND         - Set to true when lapack is found.
# LAPACK_INCLUDE_DIR   - Include directory for lapack headers.
# LAPACK_LIBRARY_DIR   - Link directory for lapack libraries.
# LAPACK_LIBRARIES     - list of lapack libraries

find_path(LAPACK_LIBRARY_DIR liblapack.a
  ${NeuroLib_LIBS_PATH}
  /usr/lib /opt/local/lib /usr/local/lib /tools/lib  /usr/lib64
  DOC " directory for lapack, where liblapack can be found"
  )

if(LAPACK_LIBRARY_DIR)
  set(LAPACK_FOUND 1)
  set(LAPACK_INCLUDE_DIR ${LAPACK_LIBRARY_DIR}/../include)

  if(WIN32)
    set(LAPACK_LIBRARIES
      lapack
      blas
      )
  endif(WIN32)
  if(CMAKE_SYSTEM MATCHES "Darwin*")
    set(LAPACK_LIBRARIES
      lapack
      blas
      )
  endif(CMAKE_SYSTEM MATCHES "Darwin*")
  if(CMAKE_SYSTEM MATCHES "Sun*")
    set(LAPACK_LIBRARIES
      lapack
      blas
      I77
      F77
      )
  endif(CMAKE_SYSTEM MATCHES "Sun*")
  if(CMAKE_SYSTEM MATCHES "Linux*")
	exec_program(gcc ARGS --version OUTPUT_VARIABLE CMAKE_C_COMPILER_VERSION)
	if(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
		set(GCC_FORTRAN_LIBRARY gfortran)
	else(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
		set(GCC_FORTRAN_LIBRARY g2c)
	endif(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
    set(LAPACK_LIBRARIES
      lapack
      blas
      ${GCC_FORTRAN_LIBRARY}
      )
  endif(CMAKE_SYSTEM MATCHES "Linux*")

else(LAPACK_LIBRARY_DIR)
  message(FATAL_ERROR "LAPACK library not found!\n" "Please go to http://www.ia.unc.edu/dev/tutorials/InstallLib")
  set(LAPACK_INCLUDE_DIR NOT_FOUND)
  set(LAPACK_FOUND 0)
endif(LAPACK_LIBRARY_DIR)

include_directories(${LAPACK_INCLUDE_DIR})
link_directories(${LAPACK_LIBRARY_DIR})
link_libraries(${LAPACK_LIBRARIES})
