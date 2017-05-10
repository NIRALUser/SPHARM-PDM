/*=========================================================================

  Program:   NeuroLib (asc2vtk command line tool: convert FreeSurfer .asc files to .vtk files)
  Language:  C++
  Date:      $Date: 2010/12/06 17:30:59 $
  Version:   $Revision: 1.2 $
  Author:    Cassian Marc (cassian.marc@cpe.fr)

  Copyright (c)  Cassian Marc. All rights reserved.
  See NeuroLibCopyright.txt or http://www.ia.unc.edu/dev/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

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
} t_polygon;

int main(int argc, const char **argv)
{	
  t_point point;
  t_polygon polygon;
  char asc_title[1000];
  string vtk_title="# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET POLYDATA\nPOINTS ";
  string vtk_float=" float\n";
  string vtk_polygon="POLYGONS ";
  int n_points, n_polygons;
  float temp;
  int i,j=0;
  
  // Check arguments (input and output file needed)
  if (argc !=3) 
    {
      cout << endl << endl << "####  Wrong number of arguments  ####" << endl << endl;
      cout << "example :  ./asc2vtk  inputfilepath/name.asc  outputfilepath/name.vtk" << endl << endl;
      exit(0);
    }
  
  // Open input and output files
  ifstream ifs(argv[1]);
  ofstream ofs(argv[2]);
  
  // Get 1st line : title
  ifs.getline(asc_title,1000);
  
  // Get 2nd line : number of points and polygons
  ifs >> n_points >> n_polygons;
  
  // Set title of vtk file with number of points
  ofs << vtk_title << n_points << vtk_float;
  
  // Copy all vertices in output file
  // in .vtk format each line has the coords of 3 points
  for(i=0;i<n_points;i++)
    {
      if(j<2)
	{
	  ifs >> point.x >> point.y >> point.z >> temp; 
	  ofs << point.x <<" "<< point.y <<" "<< point.z <<" ";
	  j++;
	}
      else
	{
	  ifs >> point.x >> point.y >> point.z >> temp;
	  ofs << point.x <<" "<< point.y <<" "<< point.z <<'\n';
	  j=0;
	}
      /* temp is used only to avoid problems with ifstream operator >>
	 As there are 4 "datas" on each .asc line, we should get the 4.
	 The patch given by the 4th data isn't used by .vtk so it get lost (int temp) */
      
    }
  
  // .vtk format separate vertices from polygones with the numbers of polygons and polygons*4
  if(j!=0)
    ofs << '\n';
  ofs << vtk_polygon <<" "<< n_polygons <<" "<< n_polygons*4 <<'\n';
  
  // Copy all polygons in output file
  for(i=0;i<n_polygons;i++)
    {
      ifs >> polygon.x >> polygon.y >> polygon.z >> temp;
      ofs << 3 <<" "<< polygon.x <<" "<< polygon.y <<" "<< polygon.z <<'\n';
    }
  
  ofs << '\n';	
  ifs.close();
  ofs.close();
  return 0;
}
