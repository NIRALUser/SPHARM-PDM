#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdlib>
#include <vtkXMLPolyDataReader.h>
#include <vtkVersion.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSphereSource.h>
#include "vtkPolyDataMapper.h"
#include "vtkPoints.h"
#include <vtkSTLReader.h>
#include <vtkDecimatePro.h>
#include <vtkSmartPointer.h> 
#include <vtkGenericDataObjectReader.h>
using namespace std;

int main(int argc, const char* argv[])
{

if (argc != 3) 
{
	cout << "2 arguments : 2 image to compare " <<endl;
	return EXIT_FAILURE;
}
 std::string inputFilename1 = argv[1];
 std::string inputFilename2 = argv[2];

vtkSmartPointer<vtkPolyDataReader> polyReader1 =
vtkSmartPointer<vtkPolyDataReader>::New();
std::cout<<inputFilename1<<std::endl;
polyReader1->SetFileName(inputFilename1.c_str());
polyReader1->Update();
vtkSmartPointer<vtkPolyData> vtkpolydata1 =
vtkSmartPointer<vtkPolyData>::New(); 
vtkpolydata1= polyReader1->GetOutput();
 

vtkSmartPointer<vtkPolyDataReader> polyReader2 =
vtkSmartPointer<vtkPolyDataReader>::New();
std::cout<<inputFilename2<<std::endl;
polyReader2->SetFileName(inputFilename2.c_str());
polyReader2->Update();
vtkSmartPointer<vtkPolyData> vtkpolydata2 =
vtkSmartPointer<vtkPolyData>::New(); 
vtkpolydata2= polyReader2->GetOutput();

if (vtkpolydata2->GetNumberOfCells() != vtkpolydata1->GetNumberOfCells()) return EXIT_FAILURE;
if (vtkpolydata2->GetNumberOfStrips() != vtkpolydata1->GetNumberOfStrips()) return EXIT_FAILURE;
if (vtkpolydata2->GetNumberOfPolys() != vtkpolydata1->GetNumberOfPolys()) return EXIT_FAILURE;
if (vtkpolydata2->GetNumberOfLines() != vtkpolydata1->GetNumberOfLines()) return EXIT_FAILURE;
if (vtkpolydata2->GetNumberOfVerts() != vtkpolydata1->GetNumberOfVerts()) return EXIT_FAILURE;
if (vtkpolydata2->GetNumberOfPieces() != vtkpolydata1->GetNumberOfPieces()) return EXIT_FAILURE;
if (vtkpolydata2->GetMaxCellSize() != vtkpolydata1->GetMaxCellSize()) return EXIT_FAILURE;

return EXIT_SUCCESS;
}
