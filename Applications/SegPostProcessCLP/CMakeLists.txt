PROJECT(SegPostProcess)

#NeuroLib_USE_ITK()
## NeuroLib_USE_SHAPE()

FIND_PACKAGE(GenerateCLP QUIET)
IF (GenerateCLP_FOUND)

   INCLUDE(${GenerateCLP_USE_FILE})
   SET ( SegPostProcessCLP_SOURCE SegPostProcessCLP.cxx )
   GENERATECLP(SegPostProcessCLP_SOURCE SegPostProcessCLP.xml)
   ADD_EXECUTABLE ( SegPostProcessCLP ${SegPostProcessCLP_SOURCE})
   TARGET_LINK_LIBRARIES ( SegPostProcessCLP ITKIO ITKBasicFilters ITKCommon)

ELSE(GenerateCLP_FOUND)
   MESSAGE(WARNING, " GenerateCLP not found. SegPostProcessCLP will not be build. If you also want to compile the CLP version (for a Slicer use), please set GenerateCLP_DIR.")
ENDIF(GenerateCLP_FOUND)

