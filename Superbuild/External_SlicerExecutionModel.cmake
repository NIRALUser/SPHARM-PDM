if(DEFINED SlicerExecutionModel_DIR AND NOT EXISTS ${SlicerExecutionModel_DIR})
  message(FATAL_ERROR "SlicerExecutionModel_DIR variable is defined but corresponds to non-existing directory")
endif()

if(NOT DEFINED SlicerExecutionModel_DIR)

set(proj SlicerExecutionModel)
    ##HACK: Need to replace with Slicer github repository.
    ExternalProject_Add(${proj}
      SVN_REPOSITORY "http://svn.slicer.org/Slicer3/trunk/Libs/SlicerExecutionModel"
      SOURCE_DIR ${proj}
      BINARY_DIR ${proj}-build
      DEPENDS ${SlicerExecutionModel_DEPENDENCIES}
      CMAKE_GENERATOR ${gen}
      CMAKE_ARGS
      --no-warn-unused-cli
        ${LOCAL_CMAKE_BUILD_OPTIONS}
        -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX}
        -DCMAKE_CXX_COMPILER:STRING=${CMAKE_CXX_COMPILER}
        -DCMAKE_CXX_COMPILER_ARG1:STRING=${CMAKE_CXX_COMPILER_ARG1}
        -DCMAKE_C_COMPILER:STRING=${CMAKE_C_COMPILER}
        -DCMAKE_C_COMPILER_ARG1:STRING=${CMAKE_C_COMPILER_ARG1}
        -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
        -DCMAKE_C_FLAGS:STRING=${CMAKE_C_FLAGS}
        -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
        -DBUILD_EXAMPLES:BOOL=OFF
        -DBUILD_TESTING:BOOL=OFF
        -DBUILD_SHARED_LIBS:BOOL=${BUILD_SHARED_LIBS}
        -DITK_DIR:PATH=${ITK_DIR}
    )
    set(GenerateCLP_DIR ${CMAKE_INSTALL_PREFIX}/lib/GenerateCLP)
    set(GenerateCLP_DEPEND "${proj}")
    set(ModuleDescriptionParser_DIR
      ${CMAKE_INSTALL_PREFIX}/lib/ModuleDescriptionParser)
    message(STATUS "GenerateCLP_DIR=${GenerateCLP_DIR}")
    set(SlicerExecutionModel_DEPEND "${proj}")
endif(NOT DEFINED SlicerExecutionModel_DIR)
