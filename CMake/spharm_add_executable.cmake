#
# This function builds and adds install rule for a standalone SPHARM-PDM executable
#
# It must be used for standalone SPHARM-PDM executables that are not built using
# the SEMMacroBuildCLI. This applies for executables that do not provide an XML
# file describing their interface.
#
# Usage:
#
#   spharm_add_executable(
#     NAME name
#     [SRC aname.cpp]
#     [ADDITIONAL_SRCS other1.cpp [other2.cpp [...]]]
#     [TARGET_LIBRARIES lib1 [lib2 [...]]]
#     [NO_INSTALL]
#     )
#
# If only the NAME is specified, the function will lookup for source files like name.cpp, name.cxx,
# name.cc, ... (complete list of extensions considered are the one associated with add_executable CMake command)
#
# If the main executable source file has different name, the SRC allows to explicitly specify it.
#
# ADDITIONAL_SRCS allows to list any other source files that should be compiled in.
#
# TARGET_LIBRARIES allows to list all libraries to link the executable with.
#
# NO_INSTALL allows to skip the install rule
#

function(spharm_add_executable)
  set(options NO_INSTALL)
  set(oneValueArgs NAME SRC)
  set(multiValueArgs ADDITIONAL_SRCS TARGET_LIBRARIES)
  cmake_parse_arguments(MY "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(MY_UNPARSED_ARGUMENTS)
    message(AUTHOR_WARNING "Unparsed arguments given [${MY_UNPARSED_ARGUMENTS}]")
  endif()

  # Sanity checks
  set(expected_defined_vars
    LOCAL_PROJECT_NAME
    ${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION
    ${LOCAL_PROJECT_NAME}_CLI_EXECUTABLE_LINK_FLAGS
    )
  foreach(var ${expected_defined_vars})
    if(NOT DEFINED ${var})
      message(FATAL_ERROR "Variable ${var} is not defined !")
    endif()
  endforeach()

  message(STATUS "Configuring executable: ${MY_NAME}")
  set(_src ${MY_NAME}.cxx)
  if(DEFINED MY_SRC)
    set(_src ${MY_SRC})
  endif()
  add_executable(${MY_NAME} ${_src} ${MY_ADDITIONAL_SRCS})
  if(MY_TARGET_LIBRARIES)
    target_link_libraries(${MY_NAME} ${MY_TARGET_LIBRARIES})
  endif()

  if(NOT "${${LOCAL_PROJECT_NAME}_CLI_EXECUTABLE_LINK_FLAGS}" STREQUAL "")
    set_target_properties(${MY_NAME} PROPERTIES LINK_FLAGS ${${LOCAL_PROJECT_NAME}_CLI_EXECUTABLE_LINK_FLAGS})
  endif()

  if(NOT MY_NO_INSTALL)
    install(TARGETS ${MY_NAME}
      RUNTIME DESTINATION ${${LOCAL_PROJECT_NAME}_CLI_INSTALL_RUNTIME_DESTINATION} COMPONENT Runtime
    )
  endif()
endfunction()
