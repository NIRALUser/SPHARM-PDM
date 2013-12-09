###Declares a list that contains the name of all the variables set in SlicerConfig.cmake
set( SlicerConfigVarList
Slicer_USE_BatchMake
Slicer_USE_MIDAS
Slicer_USE_NUMPY
Slicer_USE_OpenIGTLink
Slicer_USE_PYTHONQT
Slicer_USE_PYTHONQT_WITH_TCL
Slicer_USE_QtTesting 
Slicer_BUILD_CLI_SUPPORT
Slicer_BUILD_DICOM_SUPPORT
Slicer_BUILD_DIFFUSION_SUPPORT
Slicer_BUILD_EXTENSIONMANAGER_SUPPORT
Slicer_BUILD_TESTING
Slicer_REQUIRED_QT_VERSION
Slicer_REQUIRED_QT_MODULES
Slicer_CMAKE_DIR
Slicer_EXTENSIONS_CMAKE_DIR
Slicer_LAUNCHER_EXECUTABLE
Slicer_LAUNCH_COMMAND
Slicer_WC_REVISION
Slicer_WC_URL
Slicer_WC_ROOT
Slicer_OS
Slicer_ARCHITECTURE
Slicer_MAIN_PROJECT
Slicer_MAIN_PROJECT_APPLICATION_NAME
Slicer_LICENSE_FILE
Slicer_README_FILE
Slicer_CXX_MODULE_TEST_TEMPLATES_DIR
Slicer_PYTHON_MODULE_TEST_TEMPLATES_DIR
Slicer_EXTENSION_CPACK
Slicer_EXTENSION_CPACK_BUNDLE_FIXUP
Slicer_EXTENSIONS_DIRNAME
Slicer_BUILD_CLI
Slicer_BUILD_QTLOADABLEMODULES
Slicer_BUILD_QTSCRIPTEDMODULES
Slicer_BUILD_SHARED
Slicer_LIBRARY_PROPERTIES
Slicer_EXPORT_HEADER_TEMPLATE
Slicer_LOGOS_RESOURCE
Slicer_HOME
Slicer_CORE_LIBRARY
Slicer_GUI_LIBRARY
MRML_LIBRARIES
Slicer_Libs_VTK_WRAPPED_LIBRARIES
Slicer_Libs_INCLUDE_DIRS
Slicer_Base_INCLUDE_DIRS
Slicer_ModuleLogic_INCLUDE_DIRS
Slicer_ModuleMRML_INCLUDE_DIRS
Slicer_ModuleWidgets_INCLUDE_DIRS
Slicer_USE_FILE
Slicer_BINARY_INNER_SUBDIR
Slicer_BIN_DIR
Slicer_LIB_DIR
Slicer_INCLUDE_DIR
Slicer_SHARE_DIR
Slicer_ITKFACTORIES_DIR
Slicer_CLIMODULES_SUBDIR
Slicer_CLIMODULES_BIN_DIR
Slicer_CLIMODULES_LIB_DIR
Slicer_CLIMODULES_SHARE_DIR
Slicer_QTLOADABLEMODULES_SUBDIR
Slicer_QTLOADABLEMODULES_BIN_DIR
Slicer_QTLOADABLEMODULES_LIB_DIR
Slicer_QTLOADABLEMODULES_INCLUDE_DIR
Slicer_QTLOADABLEMODULES_SHARE_DIR
Slicer_QTLOADABLEMODULES_PYTHON_LIB_DIR
Slicer_USE_PYTHONQT
Slicer_QTSCRIPTEDMODULES_SUBDIR
Slicer_QTSCRIPTEDMODULES_BIN_DIR
Slicer_QTSCRIPTEDMODULES_LIB_DIR
Slicer_QTSCRIPTEDMODULES_INCLUDE_DIR
Slicer_QTSCRIPTEDMODULES_SHARE_DIR
Slicer_INSTALL_ROOT
Slicer_INSTALL_BIN_DIR
Slicer_INSTALL_LIB_DIR
Slicer_INSTALL_INCLUDE_DIR
Slicer_INSTALL_SHARE_DIR
Slicer_INSTALL_ITKFACTORIES_DIR
Slicer_INSTALL_CLIMODULES_BIN_DIR
Slicer_INSTALL_CLIMODULES_LIB_DIR
Slicer_INSTALL_CLIMODULES_SHARE_DIR
Slicer_INSTALL_QTLOADABLEMODULES_BIN_DIR
Slicer_INSTALL_QTLOADABLEMODULES_LIB_DIR
Slicer_INSTALL_QTLOADABLEMODULES_PYTHON_LIB_DIR
Slicer_INSTALL_QTLOADABLEMODULES_INCLUDE_DIR
Slicer_INSTALL_QTLOADABLEMODULES_SHARE_DIR
Slicer_USE_PYTHONQT
Slicer_INSTALL_QTSCRIPTEDMODULES_BIN_DIR
Slicer_INSTALL_QTSCRIPTEDMODULES_LIB_DIR
Slicer_INSTALL_QTSCRIPTEDMODULES_INCLUDE_DIR
Slicer_INSTALL_QTSCRIPTEDMODULES_SHARE_DIR
Slicer_INSTALL_PREFIX
Slicer_ExternalData_OBJECT_STORES
Slicer_ExternalData_URL_TEMPLATES
QT_QMAKE_EXECUTABLE
BatchMake_DIR
DCMTK_DIR
CTK_DIR
QtTesting_DIR
ITK_DIR
OpenIGTLink_DIR
PYTHON_EXECUTABLE
PYTHON_INCLUDE_DIR
PYTHON_LIBRARY
PYTHON_DEBUG_LIBRARY
qMidasAPI_DIR
SLICERLIBCURL_DIR
SlicerExecutionModel_EXTRA_INCLUDE_DIRECTORIES
SlicerExecutionModel_EXTRA_EXECUTABLE_TARGET_LIBRARIES
SlicerExecutionModel_DIR
Teem_DIR
VTK_DIR
Slicer_CMAKE_CXX_COMPILER
Slicer_CMAKE_C_COMPILER
Slicer_CMAKE_CXX_FLAGS
Slicer_CMAKE_C_FLAGS
CMAKE_C_COMPILER
CMAKE_CXX_COMPILER
CMAKE_CXX_FLAGS
CMAKE_C_FLAGS
Slicer_EXTERNAL_PROJECTS
Slicer_EXTERNAL_PROJECTS_NO_USEFILE
CMAKE_MODULE_PATH
)

###Unset variables and saves them as ${var}_TMP to be able to run without any error:
#find_package(slicer)
#include(${Slicer_USE_FILE})
#Typically used in SuperBuild.cmake before find_package(slicer)
macro( unsetForSlicer )
  set(options VERBOSE )
  set(oneValueArgs )
  set(multiValueArgs NAMES )
  CMAKE_PARSE_ARGUMENTS(TMP
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  foreach( var ${TMP_NAMES} )
    if( TMP_VERBOSE )
      message( "Unsetting ${var} - old value: ${${var}}" )
    endif()
    set( ${var}_TMP ${${var}} )
    unset( ${var} CACHE )
    unset( ${var} )
  endforeach()
endmacro()


###Reset variables that had been saved by "unsetForSlicer". Typically used in SuperBuild.cmake after
#unsetForSlicer(ITK_DIR VTK_DIR)
#find_package(slicer)
#include(${Slicer_USE_FILE})
macro( resetForSlicer )
  set(options VERBOSE )
  set(oneValueArgs )
  set(multiValueArgs NAMES )
  CMAKE_PARSE_ARGUMENTS(TMP
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  foreach( var ${TMP_NAMES} )
    if( "${${var}_TMP}" STREQUAL "" )
      message( "Does not reset ${var}. ${var}_TMP was empty" )
    else()
      set( ${var} ${${var}_TMP} CACHE PATH "${var} PATH" FORCE )
      unset( ${var}_TMP )
      if( TMP_VERBOSE )
        message( "resetting ${var} - new value: ${${var}}" )
      endif()
    endif()
  endforeach()
endmacro()

###unset all variables set by SlicerConfig.cmake except the ones explicitly specified
macro( unsetAllForSlicerBut )
  set(options VERBOSE )
  set(oneValueArgs )
  set(multiValueArgs NAMES )
  CMAKE_PARSE_ARGUMENTS(TMP
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )
  foreach( var ${TMP_NAMES} )
    if( TMP_VERBOSE )
      message( "Saving ${var} - value: ${${var}}" )
    endif()
    set( ${var}_TMP ${${var}} )
  endforeach()
  foreach( var ${SlicerConfigVarList} )
    unset( ${var} CACHE )
    unset( ${var} )
  endforeach()
  foreach( var ${TMP_NAMES} )
    if( TMP_VERBOSE )
      message( "Writing back ${var} - value: ${${var}_TMP}")
    endif()
    set( ${var} ${${var}_TMP} )
  endforeach()
endmacro()

