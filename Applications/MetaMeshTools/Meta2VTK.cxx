/*
 * convert Meta To VTK
 *
 * author:  Ipek Oguz
 *
 */

#include <itkDefaultDynamicMeshTraits.h>
#include <itkMetaMeshConverter.h>
#include "itkMeshTovtkPolyData.h"

#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPolyDataReader.h"

#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "argio.hh"

int main(int argc, const char * *argv)
{
  if( argc <= 2 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help") )
    {
    std::cout << "Meta2VTK " << endl;
    std::cout << "usage: Meta2VTK infile outfile [-v]" << endl;
    std::cout << endl;
    std::cout << "infile               input meta file" << endl;
    std::cout << "outfile              output vtk  file" << endl;
    std::cout << endl << endl;
    exit(0);
    }

  char * infile = strdup(argv[1]);
  char * outfile = strdup(argv[2]);
  bool   debug      = ipExistsArgument(argv, "-v");

  typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
  typedef itk::Mesh<double, 3, MeshTraitsType>                        itkMeshType;
  typedef itk::MeshSpatialObject<itkMeshType>                         itkMeshSOType;
  typedef itk::MetaMeshConverter<3, double, MeshTraitsType>           MeshConverterType;

  // read the data in meta format
  MeshConverterType *    itkConverter = new MeshConverterType();
  itkMeshSOType::Pointer meshSO = itkConverter->ReadMeta(infile);
  itkMeshType::Pointer   mesh = meshSO->GetMesh();
  delete (itkConverter);

  // convert to vtk format
  itkMeshTovtkPolyData * ITKVTKConverter = new itkMeshTovtkPolyData;
  ITKVTKConverter->SetInput( mesh );

  // write out the vtk mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput( ITKVTKConverter->GetOutput() );
  writer->SetFileName( outfile );
  writer->Update();

  writer->Delete();
  delete (ITKVTKConverter);

  return 0;
}
