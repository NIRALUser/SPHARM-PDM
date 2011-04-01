CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(Slicer3 QUIET NO_DEFAULT_PATH)

IF (Slicer3_FOUND)
  INCLUDE(${Slicer3_USE_FILE})
  slicer3_set_default_install_prefix_for_external_projects()

ELSE (Slicer3_FOUND)
  FIND_PACKAGE(ITK REQUIRED)
  IF (ITK_FOUND)
    INCLUDE(${USE_ITK_FILE})
  ELSE(ITK_FOUND)
    MESSAGE(FATAL_ERROR, "ITK not found. Please set ITK_DIR.")
  ENDIF (ITK_FOUND)

ENDIF (Slicer3_FOUND)



LINK_LIBRARIES( ${ITK_LIBRARIES})

SUBDIRS(
     	IO
     	SpatialObject
     	Numerics
     	Algorithms
)

