add_subdirectory(CLI)
if( ${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION )
  add_subdirectory(Scripted)
endif()
