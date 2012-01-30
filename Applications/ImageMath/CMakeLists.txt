#project(ImageMath)

add_executable( ImageMath ImageMath.cxx)
add_executable( CropImage CropImage.cxx)

target_link_libraries( ImageMath ${ITK_LIBRARIES} )
target_link_libraries( CropImage ${ITK_LIBRARIES} )

install(TARGETS
  ImageMath
  CropImage
  BUNDLE DESTINATION  bin
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib/static  )

