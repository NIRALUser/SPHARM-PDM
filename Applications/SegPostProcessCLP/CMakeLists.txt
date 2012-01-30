#project(SegPostProcessCLP)

set(SegPostProcessCLP_SOURCE SegPostProcessCLP.cxx )
GENERATECLP(SegPostProcessCLP_SOURCE SegPostProcessCLP.xml)
add_executable(SegPostProcessCLP ${SegPostProcessCLP_SOURCE})

target_link_libraries(SegPostProcessCLP ${ITK_LIBRARIES})

install(TARGETS
  SegPostProcessCLP
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )
