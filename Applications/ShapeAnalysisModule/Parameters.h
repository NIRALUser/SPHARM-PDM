#ifndef PARAMETERS_H
#define PARAMETERS_H


#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <itksys/Process.h>
#include <vector>
#include <dirent.h> 

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkVTKImageIO.h"
#include "itkImageIOBase.h"

#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPointSet.h"


#include <string>
#include <stdlib.h>
#include <itksys/Glob.hxx>

using namespace std;

class Parameters
{
  public:
  
	Parameters();
	 ~Parameters();
 
	void ReadFile(const char *_FileName);
	void SetModulePath(char*);
	void SetParameterFile(const char *_ParamFile);
	void SetDataList(vector < vector<string> > _List);

	void SetOutputDirectory(const char *_OutputDirectory);
	void SetChooseInputColumnState(bool inputColumn);
	void SetColumnVolumeFile(int column);
	void SetDataNumber(int _DataNumber);
	void SetEnforcedSpace(float _sx, float _sy, float _sz);
	void SetLabel(double _Label);
	void SetLabelState(bool _LabelState);
	void SetNumberOfIterations(int _NumIter);
	void SetSubdivLevel(int _SubdivLevel);
	void SetSPHARMDegree(int _SPHARMDegree);
	void SetGaussianFilteringState(bool _GaussianFilteringState);
	void SetVarianceBox(int _vx, int _vy, int _vz);
	void SetRegTemplateState(bool _RegTemplateState); 
	void SetFlipTemplateState(bool _FlipTemplateState); 
	void SetTemplateMState(bool _TemplateMState);
	void SetFlipTemplate(const char *_FlipTemplate);
	void SetRegTemplate(const char *_FlipTemplate);
	void SetParaOut1State(bool _ParaOut1State);
	void SetParaOut2State(bool _ParaOut2State);
	void SetFinalFlip(int _None_Flip, int _X_Flip, int _Y_Flip, int _Z_Flip, int _XY_Flip, int _YZ_Flip, int _XZ_Flip, int _XYZ_Flip);

	char* GetModulePath();
	vector<string> GetOutputFileHeaders();
	string GetVolumeFileExtension();
	int GetColumnVolumeFile();
	bool GetChooseInputColumnState();
	const char* GetParameterFile();
	int MeanTemplateExist();
	char* GetOutputDirectory();
	int GetDataNumber();
	float GetEnforcedSpaceX();
	float GetEnforcedSpaceY();
	float GetEnforcedSpaceZ();
  	double GetLabel();
	bool GetLabelState();
	int GetNumberOfIterations();
	double GetSubdivLevel();
	int GetSPHARMDegree();

	bool GetGaussianFilteringState();
   
	int GetVarianceBoxX();
	int GetVarianceBoxY();
	int GetVarianceBoxZ();
 
	bool GetFlipTemplateState();
	bool GetRegTemplateState();
	bool GetTemplateMState();
	char* GetFlipTemplate();
	char* GetRegTemplate(); 
	
	bool GetParaOut1State();
	bool GetParaOut2State();
	
	int GetFinalFlipN();
	int GetFinalFlipX();
	int GetFinalFlipY();
	int GetFinalFlipZ();
	int GetFinalFlipXY();
	int GetFinalFlipYZ();
	int GetFinalFlipXZ();
	int GetFinalFlipXYZ();
 
	string GetNthDataListValue(int line,int column);

	void  SetAllFilesName();
	char* GetAllFilesName(int);
	char* GetAllSurfSPHARMFiles(int);
	char* GetAllSurfSPHARMellalignFiles(int);
	char* GetAllSurfSPHARMprocalignFiles(int);
	char* GetAllSurfmeanSPHARMFiles(int);
	char* GetAllSurfmeanSPHARMellalignFiles(int);
	char* GetAllSurfmeanSPHARMprocalignFiles(int);
	char* GetAllPhiFiles(int);
	char* GetAllThetaFiles(int);
	void SetImageDimensions(char*);
	vector <double> GetImageDimensions();
	string GetDirectionToDisplay();

	void SetOverwriteSegPostProcess(bool);
	void SetOverwriteGenParaMesh(bool);
	void SetOverwriteParaToSPHARMMesh(bool);

	bool GetOverwriteSegPostProcess();
	bool GetOverwriteGenParaMesh();
	bool GetOverwriteParaToSPHARMMesh();

	bool DirectoryIsEmpty(const char*);
	
private:
 
	char m_ModulePath[512];
	vector < vector<string> > m_List;
	vector<string> m_OutputFileHeaders;
	char m_FilesName[512];
	char m_MRLMFile[512];
	string m_VolumeFileExtension;
	int m_ColumnVolume;
	bool m_ChooseColumnState;

	int ListSize;
  
	char m_ParamFile[512];
	char m_OutputDirectory[512];

	int m_DataNumber;
	
	float m_sx;
	float m_sy;
	float m_sz;
	
	double m_Label;
	bool m_LabelState;
	int m_NumIter;
	int m_SubdivLevel;
	int m_SPHARMDegree;

	bool m_GaussianFilteringState;
	
	int m_vx;
	int m_vy;
	int m_vz;
	
	bool m_FlipTemplateState;
	bool m_RegTemplateState;
	bool m_TemplateMState;
	char m_FlipTemplate[512];
	char m_RegTemplate[512];

	bool m_ParaOut1State;
	bool m_ParaOut2State;
	
	int m_None_Flip;
	int m_X_Flip;
	int m_Y_Flip;
	int m_Z_Flip;
	int m_XY_Flip;
	int m_YZ_Flip;
	int m_XZ_Flip;
	int m_XYZ_Flip;

	bool m_OverwriteSegPostProcess;
	bool m_OverwriteGenParaMesh;
	bool m_OverwriteParaToSPHARMMesh;

	string m_directionToDisplay;
	vector<double>m_Dims;

	char** m_AllFilesName; 
	char** surfSPHARM_Files; 
	char** surfSPHARM_ellalign_Files; 
	char** surfSPHARM_procalign_Files; 
	char* Phi_Files; 
	char* Theta_Files;
};

#endif
