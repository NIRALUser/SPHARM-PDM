
set(SparseLib_SRCS
  compcol_double.cxx
  comprow_double.cxx
  coord_double.cxx
  diagpre_double.cxx
  icpre_double.cxx
  ilupre_double.cxx
  iohb_double.cxx
  iotext_double.cxx
  mvblasc.cxx
  mvblasd.cxx
  mvblasf.cxx
  mvblasi.cxx
  mvmc.cxx
  mvmd.cxx
  mvmf.cxx
  mvmi.cxx
  mvvc.cxx
  mvvcio.cxx
  mvvd.cxx
  mvvdio.cxx
  mvvf.cxx
  mvvi.cxx
  qsort_double.cxx
  qsort_int.cxx
  spmm.cxx
  spsm.cxx
  )

add_library(SparseMatrixLib STATIC ${SparseLib_SRCS})
