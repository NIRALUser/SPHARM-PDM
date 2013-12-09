if( NOT EXTERNAL_SOURCE_DIRECTORY )
  set( EXTERNAL_SOURCE_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/ExternalSources )
endif()

if(OPT_USE_SYSTEM_BatchMake)
  find_package(BatchMake REQUIRED)
  include(${BatchMake_USE_FILE})
  set(BatchMake_DEPEND "")
  return()
endif()

set(proj BatchMake)
if(NOT DEFINED git_protocol)
  set(git_protocol "git")
endif()

ExternalProject_Add(${proj}
  GIT_REPOSITORY ${git_protocol}://batchmake.org/BatchMake.git
  GIT_TAG "8addbdb62f0135ba01ffe12ddfc32121b6d66ef5" 
  SOURCE_DIR BatchMake
  BINARY_DIR BatchMake-build
  ${cmakeversion_external_update} "${cmakeversion_external_update_value}"
  CMAKE_GENERATOR ${gen}
  CMAKE_ARGS
  ${COMMON_EXTERNAL_PROJECT_ARGS}
  -DUSE_FLTK:BOOL=OFF
  -DGRID_SUPPORT:BOOL=OFF
  -DUSE_SPLASHSCREEN:BOOL=OFF
  -DITK_DIR:PATH=${ITK_DIR}
  -DBUILD_SHARED_LIBS:BOOL=OFF
  -DBUILD_TESTING:BOOL=OFF
  -DDASHBOARD_SUPPORT:BOOL=OFF
  ${BatchMakeCURLCmakeArg}
  INSTALL_COMMAND ""
  PATCH_COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild/BatchMakePatchedZip.c ${CMAKE_CURRENT_BINARY_DIR}/BatchMake/Utilities/Zip/zip.c # No "" # Patch for windows compilation error (declaration of variable after beginning of block - "uLong year")
    DEPENDS ${ITK_DEPEND}
  )

set(BatchMake_DEPEND BatchMake )
set(BatchMake_DIR  ${CMAKE_CURRENT_BINARY_DIR}/BatchMake-build)
