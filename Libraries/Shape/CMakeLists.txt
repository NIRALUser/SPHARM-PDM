CMAKE_MINIMUM_REQUIRED(VERSION 2.6)





LINK_LIBRARIES( ${ITK_LIBRARIES})

SUBDIRS(
     	SpatialObject
     	IO
     	Numerics
     	Algorithms
)

