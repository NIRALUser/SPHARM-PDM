#ifndef DEF_PARTICLEMODULEPARAMETERS
#define DEF_PARTICLEMODULEPARAMETERS

#include <stdio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string.h>
#include <itksys/Glob.hxx>
#include <itksys/Process.h>
#include <string>
#include <stdlib.h>
#include <dirent.h>
#include <iostream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <itksys/SystemTools.hxx>
#include <strstream>
#include <itksys/SystemTools.hxx>
#include <time.h>

#include <algorithm>

#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPointSet.h"
#include "vtkDataSet.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkPolyDataWriter.h"
#include <vtkSmartPointer.h>

#include <vtkMath.h>
#include <vtkMergePoints.h>
#include <vtkPointSource.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkDelaunay2D.h>

using namespace std;

class ParticleModuleParameters
{
public:

  ParticleModuleParameters();
  ~ParticleModuleParameters();

  void SetParameterFile(const char *);

  void SetColumnMeshFile(int);

  int GetColumnMeshFile();

  void SetSmooting(double );

  double GetSmoomting();

  void SetEnforcedSpace(float, float, float );

  float GetEnforcedSpaceX();

  float GetEnforcedSpaceY();

  float GetEnforcedSpaceZ();

  void  SetEndRegularization(double );

  double GetEndRegularization();

  void SetStartRegularization(double);

  double GetStartRegularization();

  void SetOptimizationIteration(int );

  int GetOptimizationIteration();

  void SetCheckpointingInterval(int );

  int GetCheckpointingInterval();

  const char * GetParameterFile();

  void SetProcrustesScaling(int );

  void SetProcrustesInterval(int);

  void SetAdaptivityStrength(float );

  void SetRelativeWeighting(float );

  int GetProcrustesScaling();

  int GetProcrustesInterval();

  float GetAdaptivityStrength();

  float GetRelativeWeighting();

  void SetAdaptivityStrengthState(bool );

  void SetProcrustesScalingState(bool );

  void SetRelativeWeightingState(bool );

  void SetProcrustesIntervalState(bool);

  bool GetProcrustesScalingState();

  bool GetRelativeWeightingState();

  bool GetAdaptivityStrengthState();

  bool GetProcrustesIntervalState();

  void SetHorizontalGridPara(int);

  void SetVerticalGridPara(int );

  int GetHorizontalGridPara();

  int GetVerticalGridPara();

  void SetTemplate(int );

  int GetTemplate();

  void WriteCommandLine();

  int ReadFile(const char *);

  void ScaleMeta(std::string, std::string, float, float, float);

  void SetOutputDirectory(const char *);

  char * GetOutputDirectory();

  int GetDataNumber();

  void CreateLptsFiles();

  void SetLptsFilesName(std::string );

  void CreatePreprocessFiles();

  void RunParticleCorrespondencePreprocessing();

  void CommandlineParticleCorrespondencePreprocessing();

  void FindMhaFiles();

  void CreateCorrespondenceFiles();

  void CommandlineShapeWorksRun();

  void FindLptsFiles();

  void CopyVTKFiles(std::string originVTK, int );

  bool VtkInOutputDirectory(std::string );

  void OrganizeOutputDirectory();

  std::string SetVTKFilesName(std::string, std::string );

  void LptsToVTK(std::string, bool );

  void FindPostLptsFiles(std::string, std::vector<std::string> );

  void CreateVTKFiles();

  void MakePPDirectory();

  void SetVTKFilesName(std::string );

  void RunShapeWorksRun();

  void ChangeNames();

  std::string OutputDirectory;
  void VTKToLpts(std::string, std::string );

  void Run(std::vector<const char *>, bool );

  std::string SetVTKName(std::string, bool );

  void ColorMap();

  void CreateMrml();

  void SetImageDimensions(char *);

  vector<double> GetImageDimensions();

  std::string GetDirectionToDisplay();

  void CreatePsotProcessFiles();

private:

  int nbPoints;

  char   m_ParamFile[512];
  int    m_ColumnMesh;
  double m_smoothing;
  double m_endRegularization;
  double m_startRegularization;
  int    m_optimizationIteration;
  int    m_checkpointingInterval;
  int    m_ProcrustesScaling;
  int    m_ProcrustesInterval;
  float  m_AdaptivityStrength;
  float  m_RelativeWeight;
  bool   m_RelativeWeightingState;
  bool   m_AdaptivityStrengthState;
  bool   m_ProcrustesIntervalState;
  bool   m_ProcrustesScalingState;
  int    m_temp;

  vector<vector<string> >   m_List;
  int                       ListSize;
  char                      m_OutputDirectory[512];
  vector<string>            Data;
  vector<string>            Data2;
  int                       m_DataNumber;
  std::vector<const char *> m_args;

  vector<double> m_Dims;
  string         m_directionToDisplay;
  double         m_const_orientation;

  float m_sx;
  float m_sy;
  float m_sz;

  std::string              preprocess_file;
  std::vector<std::string> LptsFile;
  std::vector<std::string> mhaFile;
  std::string              Correspondence_file;
  std::vector<std::string> PostLptsFile;
  std::vector<std::string> PostLptsWptsFile;
  std::vector<std::string> PostModeFile;
  std::vector<std::string> VTKCorrespFile;
  std::vector<std::string> WPTSFile;
  std::string              PP_file;
  std::vector<std::string> VTKPrePross;
  std::vector<std::string> LPTS;
  std::vector<std::string> VTKPostPross;

  char * Convert_Double_To_CharArray(double doubleVariable);

  double string_to_double( const std::string& s );

  std::string PPOutputDirectory;
  int         m_VerticalGridParaGrid;
  int         m_HorizontalGridParaGrid;

};

#endif
