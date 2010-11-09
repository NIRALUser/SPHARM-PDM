
NeuroLib_USE_ITK() 
NeuroLib_USE_SHAPE()
NeuroLib_USE_VTK() 


ADD_EXECUTABLE( GenParaMesh  GenParaMesh.cxx)
TARGET_LINK_LIBRARIES( GenParaMesh ${ITK_LIBRARIES})

ADD_EXECUTABLE( ParaToSPHARMMesh  ParaToSPHARMMesh.cxx itkMeshTovtkPolyData.cxx)
TARGET_LINK_LIBRARIES( ParaToSPHARMMesh 
  ${ITK_LIBRARIES}
  vtkIO
  vtkHybrid
  vtkCommon
  vtkFiltering
)



FIND_PACKAGE(GenerateCLP QUIET)
IF (GenerateCLP_FOUND)
   INCLUDE(${GenerateCLP_USE_FILE})

   SET (GenParaMeshCLP_SOURCE GenParaMeshCLP.cxx itkMeshTovtkPolyData.cxx)
    GENERATECLP(GenParaMeshCLP_SOURCE GenParaMeshCLP.xml)
    ADD_EXECUTABLE (GenParaMeshCLP ${GenParaMeshCLP_SOURCE})
    TARGET_LINK_LIBRARIES ( GenParaMeshCLP
		   ${ITK_LIBRARIES}
		   vtkIO)

    SET (ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.cxx itkMeshTovtkPolyData.cxx vtkPolyDataToitkMesh.cxx)
    GENERATECLP(ParaToSPHARMMeshCLP_SOURCE ParaToSPHARMMeshCLP.xml)
    ADD_EXECUTABLE (ParaToSPHARMMeshCLP ${ParaToSPHARMMeshCLP_SOURCE})
    TARGET_LINK_LIBRARIES ( ParaToSPHARMMeshCLP
       ${ITK_LIBRARIES}
       vtkIO
       vtkHybrid
       vtkCommon
       vtkFiltering
       )
ENDIF(GenerateCLP_FOUND)

