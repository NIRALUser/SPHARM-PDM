
set(SpatialObject_SRCS
SphericalHarmonicSpatialObject.cxx
)

add_library(SpatialObject ${SpatialObject_SRCS})
target_link_libraries(SpatialObject ShapeAlgorithms )
