
#-----------------------------------------------------------------------------
find_package(ITK 4 REQUIRED)
include(${ITK_USE_FILE})

find_package(SlicerExecutionModel REQUIRED)
include(${SlicerExecutionModel_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(LAPACKE NO_MODULE REQUIRED)

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

#-----------------------------------------------------------------------------
if(UNIX AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # Install gfortran shared libraries
  set(SAVED_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so")
  set(fortran_LIBRARIES ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
  list(REMOVE_DUPLICATES fortran_LIBRARIES)
  list(REMOVE_ITEM fortran_LIBRARIES "c" "m")
  foreach(lib IN LISTS fortran_LIBRARIES)
    find_library(${lib}_LIBRARY ${lib} HINTS ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES} NO_DEFAULT_PATH)
    if(NOT ${lib}_LIBRARY)
      continue()
    endif()
    get_filename_component(${lib}_LIBRARY_REALPATH ${${lib}_LIBRARY} REALPATH)
    get_filename_component(lib_name ${${lib}_LIBRARY} NAME)
    install(FILES ${${lib}_LIBRARY_REALPATH}
      DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_LIBRARY_DESTINATION} COMPONENT RuntimeLibraries
      RENAME ${lib_name}
      )
  endforeach()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${SAVED_CMAKE_FIND_LIBRARY_SUFFIXES})
elseif(WIN32)
  # Install flang runtime libraries
  foreach(lib IN LISTS Fortran_Flang_RUNTIME_LIBS)
    install(FILES ${${lib}_LIBRARY}
      DESTINATION ${${LOCAL_PROJECT_NAME}_INSTALL_RUNTIME_DESTINATION} COMPONENT RuntimeLibraries
      )
  endforeach()
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
