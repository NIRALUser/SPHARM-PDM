 // convert FreeSurfer .asc Files to .meta Files
 // author:  Cassian MARC cassian.marc@cpe.fr
/*
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkPolyDataWriter.h"
#include "vtkFSSurfaceReader.h"
#include "vtkFSSurfaceScalarReader.h"
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>

using namespace std;


typedef struct
{
	float x;
	float y;
	float z;
} t_point;

typedef struct
{
	int x;
	int y;
	int z;
} t_cell;

int main(int argc, const char **argv)
{	
	t_point point;
	t_cell cell;
	char asc_title[1000];
        string meta_header="ObjectType = Mesh\nNDims = 3\nID = 0\nTransformMatrix = 1 0 0 0 1 0 0 0 1\nOffset = 0 0 0\nCenterOfRotation = 0 0 0\nElementSpacing = 1 1 1\nPointType = MET_FLOAT\nPointDataType = MET_FLOAT\nCellDataType = MET_FLOAT\nNCellTypes = 1\nPointDim = ID x y ...\nNPoints = ";
        string meta_points="Points = \n";
        string meta_cell_header="CellType = TRI\nNCells = ";
        string meta_cell="Cells = \n";
	int n_points, n_cells;
	float temp;
	int i,j=0;
	
	// Check arguments (input and output file needed)
	if (argc !=3) 
	{
		cout << endl << endl << "####  Wrong number of arguments  ####" << endl << endl;
		cout << "example :  ./asc2vtk  inputfilepath/name.asc  outputfilepath/name.meta" << endl << endl;
		exit(0);
	}
	
	// Open input and output files
	ifstream ifs(argv[1]);
	ofstream ofs(argv[2]);
	
	// Get 1st line : title
	ifs.getline(asc_title,1000);
	
	// Get 2nd line : number of points and cells
	ifs >> n_points >> n_cells;
	
	// Set title of vtk file with number of points
	ofs << meta_header << n_points <<'\n'<< meta_points;
	
	// Copy all vertices in output file
	// in .vtk format each line has the coords of 3 points
	for(i=0;i<n_points;i++)
	{
		ifs >> point.x >> point.y >> point.z >> temp; 
		ofs << i <<" "<< point.x <<" "<< point.y <<" "<< point.z <<" \n";
        }
	/* temp is used only to avoid problems with ifstream operator >>
	As there are 4 "datas" on each .asc line, we should get the 4.
	The patch given by the 4th data isn't used by .vtk so it get lost (int temp) */	
	
	
	// .vtk format separate vertices from celles with the numbers of cells and cells*4
        ofs << meta_cell_header << n_cells <<'\n'<< meta_cell;
	
	// Copy all cells in output file
	for(i=0;i<n_cells;i++)
	{
		ifs >> cell.x >> cell.y >> cell.z >> temp;
		ofs << i <<" "<< cell.x <<" "<< cell.y <<" "<< cell.z <<" \n";
	}	
	ifs.close();
	ofs.close();
	return 0;
}
