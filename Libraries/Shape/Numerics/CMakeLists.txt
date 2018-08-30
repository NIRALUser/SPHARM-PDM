add_library(ShapeNumerics INTERFACE)

if(NOT SPHARM-PDM_BUILD_SLICER_EXTENSION)
  install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/"
    DESTINATION include/Shape/Numerics COMPONENT Development
    FILES_MATCHING PATTERN "*.h" PATTERN "*.txx"
    )
endif()
