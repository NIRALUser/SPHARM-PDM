add_library(ShapeNumerics INTERFACE)

if(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/Numerics COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
