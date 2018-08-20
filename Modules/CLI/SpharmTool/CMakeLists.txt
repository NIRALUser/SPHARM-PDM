cmake_minimum_required(VERSION 2.6)

project(SpharmTool)

set (SPHARMTOOL SpharmTool)

set(${SPHARMTOOL}_SOURCE
SpharmTool.cxx
)

add_executable(${SPHARMTOOL} ${${SPHARMTOOL}_SOURCE})
target_link_libraries(${SPHARMTOOL} ${ITK_LIBRARIES} ShapeIO SpatialObject ShapeAlgorithms SpatialObject)

if( NOT SPHARM-PDM_BUILD_SLICER_EXTENSION )
 install(TARGETS
   SpharmTool
   BUNDLE DESTINATION  bin
   RUNTIME DESTINATION bin
   LIBRARY DESTINATION lib
   ARCHIVE DESTINATION lib/static  
 )
else()
  install(TARGETS SpharmTool DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} )
endif()
