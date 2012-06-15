CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(RadiusToMesh)

FIND_PACKAGE(VTK REQUIRED)
IF (VTK_FOUND)
  SET(VTK_USE_QVTK TRUE)
  SET(VTK_USE_GUISUPPORT TRUE)
  INCLUDE(${VTK_USE_FILE})
ELSE(VTK_FOUND)
   MESSAGE(FATAL_ERROR, "VTK not found. Please set VTK_DIR.")
ENDIF (VTK_FOUND)

FIND_PACKAGE(GenerateCLP REQUIRED)
IF (GenerateCLP_FOUND)
  INCLUDE(${GenerateCLP_USE_FILE})
ELSE(GenerateCLP_FOUND)
   MESSAGE(FATAL_ERROR, "GenerateCLP not found. Please set GenerateCLP_DIR.")
ENDIF (GenerateCLP_FOUND)

GENERATECLP(RadiusToMesh.cxx RadiusToMesh.xml)
ADD_EXECUTABLE(RadiusToMesh RadiusToMesh.cxx)

TARGET_LINK_LIBRARIES(RadiusToMesh   vtkWidgets
  vtkRendering
  vtkGraphics
  vtkImaging
  vtkIO
  vtkFiltering
  vtkCommon
  vtkHybrid
  vtksys
  vtkViews
  vtkInfovis
  vtklibxml2
  vtkDICOMParser
  vtkpng
  vtkzlib
  vtkjpeg
  vtkalglib
  vtkexpat
  vtkverdict
  vtkmetaio
  vtkNetCDF
  vtksqlite
  vtkexoIIc
  vtkftgl
  vtkfreetype)