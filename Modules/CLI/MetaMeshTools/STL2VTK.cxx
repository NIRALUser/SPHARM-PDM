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

#include <vtkVersion.h>
#include <vtkSTLReader.h>
#include "vtkPolyDataWriter.h"
#include "vtkPolyData.h"
#include <vtkSmartPointer.h>

#ifndef vtkFloatingPointType
#define vtkFloatingPointType float
#endif

int main(int argc, const char **argv)
{
  if ( argc < 2 ) 
  {
    std::cout << "Usage: " << argv[0] << " Input Output" << endl ;
    return 0 ;
  }

  std::string inputFilename = argv[1] ;
  std::string outputFilename = argv[2] ;

  // load the stl file
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName (inputFilename.c_str()) ;
  std::cout << "Reading input" << endl ;
  reader->Update() ;
  
  std::cout << "Converting mesh" << endl ;
  vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput() ;

  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New () ;
  writer->SetFileName ( outputFilename.c_str () ) ;
  #if VTK_MAJOR_VERSION > 5
  writer->SetInputData ( polyData ) ;
  #else
  writer->SetInput ( polyData ) ;
  #endif
  writer->Update () ;
  
  std::cout << "Conversion Completed" << endl ;

  return 0 ;
}
