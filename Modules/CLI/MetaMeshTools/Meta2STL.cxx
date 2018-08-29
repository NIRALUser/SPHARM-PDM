/*=========================================================================

  Author: Martin Styner

  This is a converter from the Meta ITK spatial object file into STL format
  via VTK

  The converter requires the name of the input STL and output META files
  as command line arguments.

=========================================================================*/

#include "itkMeshTovtkPolyData.h"
#include "argio.hh"
#include <iostream>
#include <string>

#include <itkMesh.h>
#include <itkMetaMeshConverter.h>
#include <itkLineCell.h>
#include <itkTriangleCell.h>
#include <itkMeshSpatialObject.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTriangleCell.h>
#include <itkDefaultDynamicMeshTraits.h>

#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTriangleFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkCellArray.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>


using namespace std;

int main(int argc, const char * *argv)
{
  if( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " Input.meta Output.stl" << endl;
    std::cout << " Writes STL triangulation  Output from a meta file" << endl;
    return 0;
    }

  string inputFilename = argv[1];
  string outputFilename = argv[2];
  // 0. read Meta mesh

  typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
  typedef itk::Mesh<double, 3, MeshTraitsType>                        itkMeshType;
  typedef itk::MeshSpatialObject<itkMeshType>                         itkMeshSOType;
  typedef itk::MetaMeshConverter<3, double, MeshTraitsType>           MeshConverterType;

  // read the data in meta format
  MeshConverterType::Pointer itkConverter = MeshConverterType::New();
  itkMeshSOType::Pointer meshSO =
    dynamic_cast<itkMeshSOType *>(itkConverter->ReadMeta(inputFilename.c_str() ).GetPointer() );
  itkMeshType::Pointer   mesh = meshSO->GetMesh();

  // convert to vtk format
  itkMeshTovtkPolyData ITKVTKConverter;
  ITKVTKConverter.SetInput( mesh );

  // 2. write vtk mesh as STL
  vtkSmartPointer<vtkSTLWriter> STLFileWriter = vtkSmartPointer<vtkSTLWriter>::New();
  // Gets and writes each Label-mesh in a different .stl file
  STLFileWriter->SetFileName(outputFilename.c_str() );
  STLFileWriter->SetInputData( ITKVTKConverter.GetOutput() );
  STLFileWriter->SetFileName( outputFilename.c_str() );
  STLFileWriter->Update();

  //  STLFileWriter->Write();
}
