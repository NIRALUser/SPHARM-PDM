cmake_minimum_required(VERSION 2.6)

link_libraries( ${ITK_LIBRARIES})

subdirs(
     	Algorithms
     	SpatialObject
     	IO
     	Numerics
)

