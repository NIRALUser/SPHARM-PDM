add_executable(SpharmTool SpharmTool.cxx)
target_link_libraries(SpharmTool
  PRIVATE
    ${ITK_LIBRARIES}
    ShapeIO SpatialObject
    ShapeAlgorithms
    SpatialObject
  )

install(TARGETS SpharmTool
  BUNDLE DESTINATION  ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} COMPONENT Runtime
  RUNTIME DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} COMPONENT Runtime
  )
