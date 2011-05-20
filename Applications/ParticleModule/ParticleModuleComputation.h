#ifndef PARTICLEMODULECOMPUTATION_H
#define PARTICLEMODULECOMPUTATION_H

#include "ParticleModuleParameters.h"


#include <strstream> 
#include <iostream>
#include <cstring>
#include <fstream>
#include <stdio.h>
#include <itksys/Process.h>
#include <vector>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <itksys/Glob.hxx>
#include <dirent.h>

class ParticleModuleComputation: public ParticleModuleParameters
{
 public:
	
	ParticleModuleComputation();
	~ParticleModuleComputation();

void Computation();



private:



};

#endif

