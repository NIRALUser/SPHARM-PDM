/*=========================================================================

  Author: Ipek Oguz

  This is a converter from STL file format into ITK spatial object file
  format (.meta).

  The converter requires the name of the input STL and output META files
  as command line arguments.

  This file includes a chunk of code borrowed from Insight Applications,
  to convert a VTK polydata object into an ITK mesh object.

=========================================================================*/

#include "argio.hh"
#include <iostream>
#include <string>

#include <vtkSTLReader.h>
#include "itkMesh.h"
#include <itkMetaMeshConverter.h>
#include "itkLineCell.h"
#include "itkTriangleCell.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPoints.h"
#include <vtkSmartPointer.h>
#include "vtkCellArray.h"
#include <vtkVersion.h>

using namespace std;

const unsigned int PointDimension   = 3;
const unsigned int MaxCellDimension = 2;

typedef itk::DefaultStaticMeshTraits<double, PointDimension,
                                     MaxCellDimension, double,
                                     double> MeshTraits;
typedef itk::Mesh<double, PointDimension, MeshTraits> MeshType;

int main(int argc, const char * *argv)
{
  if( argc < 2 )
    {
    std::cout << "Usage: " << argv[0] << " Input Output" << endl;
    return 0;
    }

  string inputFilename = argv[1];
  string outputFilename = argv[2];

  // load the stl file
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(inputFilename.c_str() );
  std::cout << "Reading input" << endl;
  reader->Update();

  std::cout << "Converting mesh" << endl;
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();

  MeshType::Pointer mesh = MeshType::New();

  // Transfer the points from the vtkPolyData into the itk::Mesh
  //
  const unsigned int numberOfPoints = polyData->GetNumberOfPoints();

  vtkSmartPointer<vtkPoints> vtkpoints = polyData->GetPoints();

  mesh->GetPoints()->Reserve( numberOfPoints );
  for( unsigned int p = 0; p < numberOfPoints; p++ )
    {
    double * apoint = vtkpoints->GetPoint( p );
    mesh->SetPoint( p, MeshType::PointType( apoint ) );
    }

  //
  // Transfer the cells from the vtkPolyData into the itk::Mesh
  //
  vtkSmartPointer<vtkCellArray> triangleStrips = polyData->GetStrips();

#if VTK_MAJOR_VERSION >= 9 || (VTK_MAJOR_VERSION >= 8 && VTK_MINOR_VERSION >= 90)
  const vtkIdType*  cellPoints;
#else
  vtkIdType*  cellPoints;
#endif
  vtkIdType   numberOfCellPoints;

  //
  // First count the total number of triangles from all the triangle strips.
  //
  unsigned int numberOfTriangles = 0;

  triangleStrips->InitTraversal();

  while( triangleStrips->GetNextCell( numberOfCellPoints, cellPoints ) )
    {
    numberOfTriangles += numberOfCellPoints - 2;
    }

  vtkSmartPointer<vtkCellArray> polygons = polyData->GetPolys();

  polygons->InitTraversal();

  while( polygons->GetNextCell( numberOfCellPoints, cellPoints ) )
    {
    if( numberOfCellPoints == 3 )
      {
      numberOfTriangles++;
      }
    }

  //
  // Reserve memory in the itk::Mesh for all those triangles
  //
  mesh->GetCells()->Reserve( numberOfTriangles );

  //
  // Copy the triangles from vtkPolyData into the itk::Mesh
  //
  //

  typedef MeshType::CellType CellType;

  typedef itk::TriangleCell<CellType> TriangleCellType;

  int cellId = 0;

  // first copy the triangle strips
  triangleStrips->InitTraversal();
  while( triangleStrips->GetNextCell( numberOfCellPoints, cellPoints ) )
    {

    unsigned int numberOfTrianglesInStrip = numberOfCellPoints - 2;

    uint64_t pointIds[3];
    pointIds[0] = cellPoints[0];
    pointIds[1] = cellPoints[1];
    pointIds[2] = cellPoints[2];
    for( unsigned int t = 0; t < numberOfTrianglesInStrip; t++ )
      {
      MeshType::CellAutoPointer c;
      TriangleCellType * tcell = new TriangleCellType;
      TriangleCellType::PointIdentifier itkPts[3];
      for (int ii = 0; ii < 3; ++ii)
        {
        itkPts[ii] = static_cast<TriangleCellType::PointIdentifier>(pointIds[ii]);
        }
      tcell->SetPointIds( itkPts );
      c.TakeOwnership( tcell );
      mesh->SetCell( cellId, c );
      cellId++;
      pointIds[0] = pointIds[1];
      pointIds[1] = pointIds[2];
      pointIds[2] = cellPoints[t + 3];
      }
    }

  // then copy the normal triangles
  polygons->InitTraversal();
  while( polygons->GetNextCell( numberOfCellPoints, cellPoints ) )
    {
    if( numberOfCellPoints != 3 ) // skip any non-triangle.
      {
      continue;
      }
    MeshType::CellAutoPointer c;
    TriangleCellType *  t = new TriangleCellType;
    TriangleCellType::PointIdentifier itkPts[3];
    for (int ii = 0; ii < numberOfCellPoints; ++ii)
      {
      itkPts[ii] = static_cast<TriangleCellType::PointIdentifier>(cellPoints[ii]);
      }
    t->SetPointIds( itkPts );
    c.TakeOwnership( t );
    mesh->SetCell( cellId, c );
    cellId++;
    }

  // write meta file
  typedef itk::MeshSpatialObject<MeshType>              itkMeshSOType;
  typedef itk::MetaMeshConverter<3, double, MeshTraits> MeshConverterType;

  std::cout << "Writing output" << endl;
  itkMeshSOType::Pointer meshSO = itkMeshSOType::New();
  meshSO->SetMesh(mesh);
  MeshConverterType::Pointer itkConverter = MeshConverterType::New();
  itkConverter->WriteMeta(meshSO, outputFilename.c_str() );

  std::cout << "Conversion Completed" << endl;

  return 0;
}
