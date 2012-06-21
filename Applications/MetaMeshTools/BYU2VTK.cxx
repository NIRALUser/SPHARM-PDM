/* 
 * convert BYU To VTK
 *
 * author:  Ipek Oguz 
 *
 */

#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkPoints.h"
#include "vtkIdList.h"
#include "vtkCellArray.h"

#include <fstream>
#include <iostream>
#include <vector>

#include "argio.hh"

int main(int argc, const char **argv)
{
  if ( argc <= 2 || ipExistsArgument ( argv, "-usage" ) || ipExistsArgument ( argv, "-help" ) )
    {
      std::cout << "Usage: " << argv[0] << "infile outfile [-v]" << endl << endl ;
      std::cout << "infile               input BYU file" << endl  ;
      std::cout << "outfile              output VTK file" << endl << endl << endl ;
      exit(0) ;
    }

  char * infile = strdup ( argv[1] ) ;
  char * outfile = strdup ( argv[2] ) ;
  bool debug = ipExistsArgument ( argv, "-v" ) ;

  // read the data in byu format, building the polydata as we go

  // read the header
  std::ifstream in ;
  in.open ( infile ) ;
  int dummy, nPts, nTris, nTris_3, nTris_repeat ;
  in >> dummy >> nPts >> nTris >> nTris_3 ;
  if ( nTris_3 != 3 * nTris ) 
    {
      std::cout << "I don't understand the header. " << nTris_3 << " is not equal to 3 * " << nTris << ". " << std::endl ;
      exit ( 0 ) ;
    }
  in >> dummy >> nTris_repeat ;
  if ( nTris_repeat != nTris )
    {
      std::cout << "I don't understand the header. " << nTris_repeat << " is not equal to " << nTris << ". " << std::endl ;
      exit ( 0 ) ;
    }

  // create the polydata
  vtkPolyData *polydata = vtkPolyData::New () ;
  vtkPoints *points = vtkPoints::New () ;
  points->SetNumberOfPoints ( nPts ) ;

  // read the vertices
  if ( debug ) 
    std::cout << "Reading vertices. " << std::endl ;
  double x, y, z ;
  for ( int i = 0 ; i < nPts ; i++ )
    {
      in >> x >> y >> z ;
      points->SetPoint ( i, x, y, z ) ;
    }
  polydata->SetPoints ( points ) ;
  if ( debug ) 
    std::cout << "Finished reading vertices." << std::endl ;

  // read the tris
  if ( debug ) 
    std::cout << "Reading triangles. " << std::endl ;
  int a, b, c ;
  vtkIdList *ids = vtkIdList::New () ;
  ids->SetNumberOfIds ( 3 ) ;
  vtkCellArray *cells = vtkCellArray::New () ;
  for ( int i = 0 ; i < nTris ; i++ )
    {
      in >> a >> b >> c ;
      ids->SetId ( 0, a - 1 ) ;
      ids->SetId ( 1, b - 1 ) ;
      ids->SetId ( 2, - c - 1 ) ;
      cells->InsertNextCell ( ids ) ;
    }
  in.close () ;
  polydata->SetPolys ( cells ) ;
  if ( debug ) 
    std:: cout << "Finished reading triangles. Writing polydata to file. " << std::endl ;

  // write out the vtk mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New () ;
  writer->SetInput ( polydata ) ;
  writer->SetFileName ( outfile ) ;
  writer->Update () ;
  writer->Delete () ;
  ids->Delete () ;
  cells->Delete () ;
  points->Delete () ;
  polydata->Delete () ;
  if ( debug ) 
    std::cout << "Success." << std::endl ;

  return 0 ;
}
