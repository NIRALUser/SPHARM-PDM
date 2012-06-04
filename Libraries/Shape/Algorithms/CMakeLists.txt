
#Link to CLAPACK
find_package(CLAPACK NO_MODULE REQUIRED)
include_directories(${CLAPACK_DIR}/../CLAPACK/INCLUDE)


if(WIN32)
  add_definitions("-DCOMPLEX=std::complex<double>")
else(WIN32)
  add_definitions('-DCOMPLEX=std::complex<double>')
endif(WIN32)

set(ShapeAlgorithms_SRCS
 SphericalHarmonicMeshSource.cxx TSurfaceNet.cxx TVoxelVolume.cxx
 EqualAreaParametricMeshNewtonIterator.cxx
 ParametricMeshToSPHARMSpatialObjectFilter.cxx
 SphericalHarmonicMedialAxisMeshSource.cxx
)

add_library(ShapeAlgorithms  STATIC ${ShapeAlgorithms_SRCS})
