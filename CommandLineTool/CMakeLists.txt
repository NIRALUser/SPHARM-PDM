# Export of the files SPHARM-PDM-parameters.ini and SPHARM-PDM.py in the SPHARM package share folder in order to run SPHARM-PDM by command-lines via SlicerSALT or Slicer
set(SPHARM-PDM_python-script ${CMAKE_SOURCE_DIR}/CommandLineTool/SPHARM-PDM.py)
set(SPHARM-PDM_parameters-file ${CMAKE_SOURCE_DIR}/CommandLineTool/SPHARM-PDM-parameters.ini)

set(SPHARM_commandline_install_dir ${Slicer_INSTALL_ROOT}${Slicer_SHARE_DIR}/CommandLineTool)

if(APPLE)
  set(SPHARM_commandline_install_dir ${Slicer_INSTALL_ROOT}${Slicer_BUNDLE_EXTENSIONS_LOCATION}${Slicer_SHARE_DIR}/CommandLineTool)
endif()

install(FILES ${SPHARM-PDM_python-script} DESTINATION ${SPHARM_commandline_install_dir} COMPONENT RuntimeLibraries)
install(FILES ${SPHARM-PDM_parameters-file} DESTINATION ${SPHARM_commandline_install_dir} COMPONENT RuntimeLibraries)
