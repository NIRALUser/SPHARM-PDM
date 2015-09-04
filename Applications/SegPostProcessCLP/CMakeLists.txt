#project(SegPostProcessCLP)

set( CLP SegPostProcessCLP )

set(${CLP}_SOURCE ${CLP}.cxx )
GENERATECLP(${CLP}_SOURCE ${CLP}.xml)
add_executable(${CLP} ${${CLP}_SOURCE})
target_link_libraries(${CLP} ${ITK_LIBRARIES})

add_library(${CLP}Lib STATIC ${${CLP}_SOURCE})
set_target_properties (${CLP}Lib PROPERTIES COMPILE_FLAGS "-Dmain=ModuleEntryPoint")
target_link_libraries(${CLP}Lib ${ITK_LIBRARIES})

if( NOT SPHARM-PDM_BUILD_SLICER_EXTENSION )
  install(TARGETS
    SegPostProcessCLP
    BUNDLE DESTINATION  bin
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib/static  
  )
else()
  install(TARGETS ${CLP}
          RUNTIME DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION}
         )
endif()
