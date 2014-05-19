#include <iostream>
#include <fstream>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkFieldData.h>
#include "RadiusToMeshCLP.h"

std::vector<std::string> ReadNameFromCSV(std::string Filename)
{
	std::vector<std::string> Names;
	std::ifstream File(Filename.c_str(), std::ios::in);
	std::string Buffer;
	getline(File,Buffer);
	while(Buffer.size()>0)
	{
		Names.push_back(Buffer.substr(0,Buffer.find_first_of(",")));
		if(Buffer.find_first_of(",")!=std::string::npos)
			Buffer=Buffer.substr(Buffer.find_first_of(",")+1,Buffer.size()-Buffer.find_first_of(",")-1);
		else
			Buffer="";
	}
	return Names;	
}

std::vector<std::vector<double> > ReadDataFromCSV(std::string Filename, std::vector<std::string> FieldsName)
{
	std::vector<std::vector<double> > Data;
	for(size_t i=0; i<FieldsName.size(); i++)
          {
          Data.push_back(std::vector<double>());
          }
	std::ifstream File(Filename.c_str(), std::ios::in);
	std::string Buffer;
	getline(File,Buffer);
	getline(File,Buffer);
	while(!File.eof())
	{
		std::istringstream iss(Buffer);
		for(size_t i=0; i<FieldsName.size(); i++)
		{
			double Value;
			Value=atof((Buffer.substr(0,Buffer.find_first_of(","))).c_str());
			Buffer=Buffer.substr(Buffer.find_first_of(",")+1,Buffer.size()-Buffer.find_first_of(",")-1);
			Data[i].push_back(Value);
		}
		getline(File,Buffer);
	}
	return Data;
}

void GetBounds(std::vector<double> Vector, double Bounds[])
{
	double min=100000,max=-1000000;
	for(size_t i=0; i<Vector.size(); i++)
	{
		if(Vector[i]<min)
			min=Vector[i];
		if(Vector[i]>max)
			max=Vector[i];
	}
	Bounds[0]=min;
	Bounds[1]=max;
}

int main(int argc, char* argv[])
{
	PARSE_ARGS;
	
	std::cout<<"Reading CSV data..."<<std::endl;
	std::vector<std::vector<double> > FieldsData;
	std::vector<std::string> FieldsName;
	FieldsName=ReadNameFromCSV(CSV);
	FieldsData=ReadDataFromCSV(CSV, FieldsName);
	std::cout<<"CSV File read successfuly."<<std::endl;
	
	std::cout<<"Reading VTK data..."<<std::endl;
	vtkSmartPointer<vtkPolyData> PolyData;
	vtkSmartPointer<vtkPolyDataReader> reader(vtkPolyDataReader::New());
	reader->SetFileName(VTK.c_str());
	if(reader->IsFilePolyData())
	{
		reader->Update();
		PolyData=reader->GetOutput();
		std::cout<<"VTK File read successfuly."<<std::endl;
	}
	else
		std::cout<<"Error reading VTK File. Check VTK Format."<<std::endl;
	
	for(size_t i=0; i<FieldsName.size(); i++)
	{
        if(PolyData->GetPoints()->GetNumberOfPoints() != (vtkIdType)FieldsData[i].size()*phiIteration)
		{
			std::cout<<"Data samples and points size doesn't match!"<<std::endl;
			std::cout<<"Abort."<<std::endl;
			return 0;
		}
	}
	
	for(size_t i=0; i<FieldsName.size(); i++)
	{
		vtkSmartPointer<vtkFloatArray> Field(vtkFloatArray::New());
		vtkSmartPointer<vtkFloatArray> FieldSlicer(vtkFloatArray::New());
		Field->SetNumberOfComponents(FieldsData[i].size()*phiIteration);
		FieldSlicer->SetNumberOfComponents(FieldsData[i].size()*phiIteration);
		for(size_t j = 0; j<FieldsData[i].size(); j++)
		{
			double* Bounds=new double[2];
			GetBounds(FieldsData[i],Bounds);
			for(int k=0; k<phiIteration; k++)
			{
				Field->InsertComponent(0,j*phiIteration+k,FieldsData[i][j]);
				double SlicerValue=(int)((FieldsData[i][j]-Bounds[0])/(Bounds[1]-Bounds[0])*100);
				FieldSlicer->InsertComponent(0,j*phiIteration+k,SlicerValue);
			}
		}
		Field->SetName(FieldsName[i].c_str());
		FieldSlicer->SetName((FieldsName[i]+"_Slicer").c_str());
		PolyData->GetPointData()->AddArray(Field);
		PolyData->GetPointData()->AddArray(FieldSlicer);
	}
	
	vtkSmartPointer<vtkPolyDataWriter> writer(vtkPolyDataWriter::New());
	writer->SetFileName(output.c_str());
  #if VTK_MAJOR_VERSION > 5
	writer->SetInputData(PolyData);
  #else
	writer->SetInput(PolyData);
  #endif
	writer->Update();
	
	return 0;
}
