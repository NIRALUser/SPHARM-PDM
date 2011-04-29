
IF(WIN32)
  ADD_DEFINITIONS("-DCOMPLEX=std::complex<double>")
ELSE(WIN32)
  ADD_DEFINITIONS('-DCOMPLEX=std::complex<double>')
ENDIF(WIN32)

SET(ShapeAlgorithms_SRCS
 SphericalHarmonicMeshSource.cxx TSurfaceNet.cxx TVoxelVolume.cxx
 EqualAreaParametricMeshNewtonIterator.cxx
 ParametricMeshToSPHARMSpatialObjectFilter.cxx       
)

ADD_LIBRARY(ShapeAlgorithms  STATIC ${ShapeAlgorithms_SRCS})
