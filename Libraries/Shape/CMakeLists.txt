cmake_minimum_required(VERSION 2.6)

link_libraries( ${ITK_LIBRARIES})

add_subdirectory(SpatialObject)
add_subdirectory(IO)
add_subdirectory(Numerics)
add_subdirectory(Algorithms)

