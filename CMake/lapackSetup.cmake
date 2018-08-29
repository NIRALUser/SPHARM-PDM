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
  endif()
  if(CMAKE_SYSTEM MATCHES "Darwin*")
    set(LAPACK_LIBRARIES
      lapack
      blas
      )
  endif()
  if(CMAKE_SYSTEM MATCHES "Sun*")
    set(LAPACK_LIBRARIES
      lapack
      blas
      I77
      F77
      )
  endif()
  if(CMAKE_SYSTEM MATCHES "Linux*")
	exec_program(gcc ARGS --version OUTPUT_VARIABLE CMAKE_C_COMPILER_VERSION)
	if(CMAKE_C_COMPILER_VERSION MATCHES ".*4\\.[0-9]\\.[0-9].*")
		set(GCC_FORTRAN_LIBRARY gfortran)
	else()
		set(GCC_FORTRAN_LIBRARY g2c)
	endif()
    set(LAPACK_LIBRARIES
      lapack
      blas
      ${GCC_FORTRAN_LIBRARY}
      )
  endif()

else()
  message(FATAL_ERROR "LAPACK library not found!\n" "Please go to http://www.ia.unc.edu/dev/tutorials/InstallLib")
  set(LAPACK_INCLUDE_DIR NOT_FOUND)
  set(LAPACK_FOUND 0)
endif()

include_directories(${LAPACK_INCLUDE_DIR})
link_directories(${LAPACK_LIBRARY_DIR})
link_libraries(${LAPACK_LIBRARIES})
