project(ImageMath)

find_package(Slicer3 QUIET NO_DEFAULT_PATH)

if(Slicer3_FOUND)
  include(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

else(Slicer3_FOUND)
  find_package(ITK REQUIRED)
  if(ITK_FOUND)
    include(${USE_ITK_FILE})
  else(ITK_FOUND)
    message(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  endif(ITK_FOUND)

endif(Slicer3_FOUND)


add_executable( ImageMath ImageMath.cxx)

add_executable( CropImage CropImage.cxx)


target_link_libraries( ImageMath ${ITK_LIBRARIES} )

target_link_libraries( CropImage ${ITK_LIBRARIES} )

if(Slicer3_FOUND)
 slicer3_set_plugins_output_path(ImageMath)
 slicer3_install_plugins(ImageMath)
 slicer3_set_plugins_output_path(CropImage)
 slicer3_install_plugins(CropImage)
endif(Slicer3_FOUND)

install(TARGETS
  ImageMath
  CropImage
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

