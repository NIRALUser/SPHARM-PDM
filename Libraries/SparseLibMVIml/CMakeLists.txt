
if(WIN32)
  add_definitions("-DCOMPLEX=std::complex<double>")
else()
  add_definitions('-DCOMPLEX=std::complex<double>')
endif()

set(SparseLib_SRCS compcol_double.cxx
  diagpre_double.cxx mvblasf.cxx mvmd.cxx mvvc.cxx
  mvvf.cxx qsort_int.cxx comprow_double.cxx
  iotext_double.cxx mvblasi.cxx mvmf.cxx
  mvvcio.cxx mvvi.cxx coord_double.cxx icpre_double.cxx
  ilupre_double.cxx mvblasc.cxx  mvmi.cxx mvvd.cxx
  spmm.cxx iohb_double.cxx mvblasd.cxx
  mvmc.cxx mvvdio.cxx qsort_double.cxx spsm.cxx
)

add_library(SparseMatrixLib STATIC  ${SparseLib_SRCS})

if( NOT SPHARM-PDM_BUILD_SLICER_EXTENSION )
  install_targets(/lib SparseMatrixLib)
  install_files(/include/Shape/Numerics "(\\.h|\\.txx)$")
endif()
