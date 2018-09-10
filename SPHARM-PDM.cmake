
#-----------------------------------------------------------------------------
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(LAPACKE CONFIG REQUIRED)

# Set Fortran_<id>_RUNTIME_LIBRARIES and CMAKE_Fortran_IMPLICIT_LINK_*
set(Fortran_COMPILER_ID "${LAPACKE_Fortran_COMPILER_ID}")
find_package(Fortran REQUIRED)

#-----------------------------------------------------------------------------
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE ON)

add_subdirectory(Libraries)
add_subdirectory(Modules)
if (${LOCAL_PROJECT_NAME}_BUILD_SLICER_EXTENSION)
 add_subdirectory(CommandLineTool)
endif()

#-----------------------------------------------------------------------------
# Testing
#-----------------------------------------------------------------------------
if(BUILD_TESTING)
  include(CTest)
  add_subdirectory(Testing)
endif()

#-----------------------------------------------------------------------------
# Packaging
#-----------------------------------------------------------------------------
set(EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS)
list(APPEND EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS "${LAPACK_DIR};LAPACK;RuntimeLibraries;/")
set(${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS "${EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS}" CACHE STRING "List of external projects to install" FORCE)

# Install fortran runtime libraries
if(DEFINED Slicer_DIR AND NOT APPLE)
  Fortran_InstallLibrary(
    FILES ${Fortran_${LAPACKE_Fortran_COMPILER_ID}_RUNTIME_LIBRARIES}
    DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
    )
endif()

#-----------------------------------------------------------------------------
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};Runtime;/")
set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};RuntimeLibraries;/")
if(${LOCAL_PROJECT_NAME}_INSTALL_DEVELOPMENT)
  set(CPACK_INSTALL_CMAKE_PROJECTS "${CPACK_INSTALL_CMAKE_PROJECTS};${CMAKE_BINARY_DIR};${EXTENSION_NAME};Development;/")
endif()
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS}")

#-----------------------------------------------------------------------------
if(DEFINED Slicer_DIR)
  include(${Slicer_EXTENSION_GENERATE_CONFIG})
  include(${Slicer_EXTENSION_CPACK})
endif()
