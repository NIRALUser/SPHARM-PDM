
IF(WIN32)
  ADD_DEFINITIONS("-DCOMPLEX=std::complex<double>")
ELSE(WIN32)
  ADD_DEFINITIONS('-DCOMPLEX=std::complex<double>')
ENDIF(WIN32)

SET(SparseLib_SRCS compcol_double.cxx
  diagpre_double.cxx mvblasf.cxx mvmd.cxx mvvc.cxx
  mvvf.cxx qsort_int.cxx comprow_double.cxx
  iotext_double.cxx mvblasi.cxx mvmf.cxx
  mvvcio.cxx mvvi.cxx coord_double.cxx icpre_double.cxx
  ilupre_double.cxx mvblasc.cxx  mvmi.cxx mvvd.cxx 
  spmm.cxx iohb_double.cxx mvblasd.cxx
  mvmc.cxx mvvdio.cxx qsort_double.cxx spsm.cxx
)

ADD_LIBRARY(SparseMatrixLib STATIC  ${SparseLib_SRCS})

INSTALL_TARGETS(/lib SparseMatrixLib)

INSTALL_FILES(/include/Shape/Numerics "(\\.h|\\.txx)$")

