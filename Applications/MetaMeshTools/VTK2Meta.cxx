/*
 * convert VTK Polydata To Meta
 *
 * author:  Ipek Oguz
 *
 */

#include <itkDefaultDynamicMeshTraits.h>
#include <itkMetaMeshConverter.h>
#include "vtkPolyDataToitkMesh.h"

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
    std::cout << "VTK2Meta " << endl;
    std::cout << "usage: VTK2Meta infile outfile [-v]" << endl;
    std::cout << endl;
    std::cout << "infile               input vtk file" << endl;
    std::cout << "outfile              output meta file" << endl;
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

  // read in the vtk polydata file
  vtkPolyDataReader *reader = vtkPolyDataReader::New();
  reader->SetFileName( infile );
  reader->Update();
  vtkPolyData *polydata = reader->GetOutput();

  // convert to itk mesh data structure
  vtkPolyDataToitkMesh *vtkItkConverter = new vtkPolyDataToitkMesh();
  vtkItkConverter->SetInput( polydata );

  if (debug) {
    std::cout << "converting mesh " << infile <<std::endl
  }

  // write out the itk meta mesh file
  itkMeshSOType::Pointer meshSO = itkMeshSOType::New();
  meshSO->SetMesh( vtkItkConverter->GetOutput() );
  MeshConverterType::Pointer itkConverter = MeshConverterType::New();
  itkConverter->WriteMeta( meshSO, outfile );

  // cleanup memory
  reader->Delete();

  return 0;
}
