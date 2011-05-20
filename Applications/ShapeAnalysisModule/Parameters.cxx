 #include "Parameters.h"
#include <stdio.h>
#include <string.h>
#include <limits>
#include <sstream>
#include <string>
using namespace std;

 Parameters::Parameters()
 {
 }

Parameters::~Parameters()
 {}



// Set the slicer module path
void Parameters::SetModulePath(char *path)
{
	strcpy(m_ModulePath,path);
	char* p=strrchr(m_ModulePath,'/');
	p[1]='\0';
	strcat(m_ModulePath,"BatchMake_Applications/");
}

// Get the slicer module path
char* Parameters::GetModulePath()
{
	return m_ModulePath;
}

// Read the csv input file
void Parameters::ReadFile(const char *_FileName)
{
	int nbLine=0;
	int column=0;
	int nbColumn=0;	
	bool volumeFileFound=false;
	vector <vector<string> > List;
	
	ifstream file(_FileName);
	string Line;
	
	if(file)
	{
		//read line by line
		while(getline(file,Line))
		{
			istringstream s(Line);
			string value;
			
			vector<string> Data;

			if(!Line.empty())
			{
				while(getline(s,value,','))
				{	
					Data.push_back(value);
					if(nbLine==0)
					{
						m_OutputFileHeaders.push_back(value);// Headers of the input file useed for the ouput file
					}
					if(nbLine==1 && volumeFileFound==false) 
					{	
						if(!m_ChooseColumnState)// the user doesn't choose the input file's column
						{
						char *volume=const_cast<char*>(strchr(value.c_str(),'.'));
						if(volume ) //the first volume file is selected automatically
						{
							if(strcmp(volume,".gipl.gz")==0 || strcmp(volume,".gipl")==0 || strcmp(volume,".mka")==0 || strcmp(volume,".mhd")==0 || strcmp(volume,".nrrd") || strcmp(volume,".nhdr")==0 || strcmp(volume,".nii")==0 || strcmp(volume,".nii")==0 || strcmp(volume,".nii.gz")==0 || strcmp(volume,".hdr")==0 || strcmp(volume,".mnc")==0 ) 
							{
								m_ColumnVolume=column;
								m_VolumeFileExtension=volume;
								volumeFileFound=true;
							}
							else 
							{
								cerr<<"Volume file not found in the input file."<<endl;
								exit(0);
							}	
						}	
						}
						else // the user chooses the input file's column to select a volume file
						{
							if(m_ColumnVolume==column) 
							{
								char *volume=const_cast<char*>(strchr(value.c_str(),'.'));
								if(volume )
								{
									m_VolumeFileExtension=volume;
									volumeFileFound=true;
								}
								else 
								{
									cerr<<"The column you chose doesn't containt a volume file. Please enter a new column number."<<endl;
									exit(0);
								}
							}
						}
					}
					column++;	
						
				}
				nbLine++;
				nbColumn=column;
				column=0;
				List.push_back(Data);
			}
		}
		SetDataNumber(nbLine-1);	
	}
	else
	{ 
		cerr<<"Error read input File"<<endl;
		exit(0);
	}

	SetDataList(List);
}

// Get the headers of the output file
vector<string> Parameters::GetOutputFileHeaders()
{
	return m_OutputFileHeaders;

}

//create a file with all the names of the generated files
void Parameters::SetParameterFile(const char *_ParamFile)
{
   	std::strcpy(m_ParamFile, _ParamFile);
}

const char* Parameters::GetParameterFile()
{
	return m_ParamFile;
}

string Parameters::GetVolumeFileExtension()
{
	return m_VolumeFileExtension;
}

void Parameters::SetChooseInputColumnState(bool inputColumn)
{
	m_ChooseColumnState=inputColumn;
}

bool Parameters::GetChooseInputColumnState()
{
	return m_ChooseColumnState;
}

void Parameters::SetColumnVolumeFile(int column)
{	
	 m_ColumnVolume=column;
}

int Parameters::GetColumnVolumeFile()
{
	return m_ColumnVolume;
}

void Parameters::SetOutputDirectory(const char *_OutputDirectory)
{
 	 std::strcpy(m_OutputDirectory, _OutputDirectory);
}

char* Parameters::GetOutputDirectory()
{
  	return m_OutputDirectory;
}

void Parameters::SetDataNumber(int _DataNumber)
{
  	m_DataNumber = _DataNumber;
}

int Parameters::GetDataNumber()
{
  	return m_DataNumber;
}

void Parameters::SetEnforcedSpace(float _sx, float _sy, float _sz)
{
	if (_sx >= 0.0)
	{m_sx = _sx;}
	
	if (_sy >= 0.0)
	{m_sy = _sy;}
	
	if (_sz >= 0.0)
	{m_sz = _sz;}
}

float Parameters::GetEnforcedSpaceX()
{
  	return m_sx;
}

float Parameters::GetEnforcedSpaceY()
{
  	return m_sy;
}

float Parameters::GetEnforcedSpaceZ()
{
  	return m_sz;
}

void Parameters::SetLabel(double _Label)
{
  	m_Label = _Label;
}

double Parameters::GetLabel()
{
  	return m_Label;
}

void Parameters::SetLabelState(bool _LabelState)
{
  	m_LabelState = _LabelState;
}

bool Parameters::GetLabelState()
{
  	return m_LabelState;
}

void Parameters::SetNumberOfIterations(int _NumIter)
{
  	m_NumIter = _NumIter;
}

int Parameters::GetNumberOfIterations()
{
  	return m_NumIter;
}

void Parameters::SetSubdivLevel(int _SubdivLevel)
{
  	m_SubdivLevel = _SubdivLevel;
}

double Parameters::GetSubdivLevel()
{
  	return m_SubdivLevel;
}


void Parameters::SetSPHARMDegree(int _SPHARMDegree)
{
  	m_SPHARMDegree = _SPHARMDegree;
}

int Parameters::GetSPHARMDegree()
{
  	return m_SPHARMDegree;
}


void Parameters::SetGaussianFilteringState(bool _GaussianFilteringState)
{
 	 m_GaussianFilteringState = _GaussianFilteringState;
}

bool Parameters::GetGaussianFilteringState()
{
  	return m_GaussianFilteringState;
}

void Parameters::SetVarianceBox(int _vx,int _vy,int _vz)
{
	if (_vx >= 0)
	m_vx = _vx;
	
	if (_vy >= 0)
	m_vy = _vy;
	
	if (_vz >= 0)
	m_vz = _vz;
}

int Parameters::GetVarianceBoxX()
{
  	return m_vx;
}

int Parameters::GetVarianceBoxY()
{
  	return m_vy;
}

int Parameters::GetVarianceBoxZ()
{
  	return m_vz;
}

void Parameters::SetFlipTemplateState(bool _FlipTemplateState)
{
  	m_FlipTemplateState = _FlipTemplateState;
}

void Parameters::SetRegTemplateState(bool _RegTemplateState)
{
  	m_RegTemplateState = _RegTemplateState;
}


bool Parameters::GetRegTemplateState()
{
  	return m_RegTemplateState;
}
bool Parameters::GetFlipTemplateState()
{
  	return m_FlipTemplateState;
}

void Parameters::SetTemplateMState(bool _TemplateMState)
{
  	m_TemplateMState = _TemplateMState;
}

bool Parameters::GetTemplateMState()
{
  	return m_TemplateMState;
}

void Parameters::SetFlipTemplate(const char *_FlipTemplate)
{
  	std::strcpy(m_FlipTemplate, _FlipTemplate);
}

char* Parameters::GetFlipTemplate()
{
  	return m_FlipTemplate;
}

void Parameters::SetRegTemplate(const char *_RegTemplate)
{
  	std::strcpy(m_RegTemplate, _RegTemplate);
}

char* Parameters::GetRegTemplate()
{
  	return m_RegTemplate;
}


void Parameters::SetParaOut1State(bool _ParaOut1State)
{
  	m_ParaOut1State = _ParaOut1State;
}
bool Parameters::GetParaOut1State()
{
  	return m_ParaOut1State;
}

void Parameters::SetParaOut2State(bool _ParaOut2State)
{
  	m_ParaOut2State = _ParaOut2State;
}

bool Parameters::GetParaOut2State()
{
  	return m_ParaOut2State;
}

void Parameters::SetFinalFlip(int _None_Flip, int _X_Flip, int _Y_Flip, int _Z_Flip, int _XY_Flip, int _YZ_Flip, int _XZ_Flip, int _XYZ_Flip)
{
	m_None_Flip = _None_Flip;
	m_X_Flip = _X_Flip;
	m_Y_Flip = _Y_Flip;
	m_Z_Flip = _Z_Flip;
	m_XY_Flip = _XY_Flip;
	m_YZ_Flip = _YZ_Flip;
	m_XZ_Flip = _XZ_Flip;
	m_XYZ_Flip = _XYZ_Flip;
}

int Parameters::GetFinalFlipN()
{
  	return m_None_Flip;
}

int Parameters::GetFinalFlipX()
{
  	return m_X_Flip;
}

int Parameters::GetFinalFlipY()
{
  	return m_Y_Flip;
}

int Parameters::GetFinalFlipZ()
{
  	return m_Z_Flip;
}

int Parameters::GetFinalFlipXY()
{
  	return m_XY_Flip;
}

int Parameters::GetFinalFlipYZ()
{
  	return m_YZ_Flip;
}

int Parameters::GetFinalFlipXZ()
{
  	return m_XZ_Flip;
}

int Parameters::GetFinalFlipXYZ()
{
  	return m_XYZ_Flip;
}

int Parameters::SetHorizontalGridPara(int _HorizontalGridPara)
{
	m_HorizontalGridParaGrid = _HorizontalGridPara;
//std::cout<<"m_HorizontalGridParaGrid "<<m_HorizontalGridParaGrid<<std::endl;
}
int Parameters::SetVerticalGridPara(int _VerticalGridPara)
{
	 m_VerticalGridParaGrid = _VerticalGridPara;
//std::cout<<"m_VerticalGridParaGrid "<<m_VerticalGridParaGrid<<std::endl;
}

int Parameters::GetHorizontalGridPara()
{
	return	m_HorizontalGridParaGrid;
}
int Parameters::GetVerticalGridPara()
{
	return m_VerticalGridParaGrid;
}

void Parameters::SetParticlesState(bool _DoParticlesCorrespondence)
{
	m_DoParticlesCorrespondence=_DoParticlesCorrespondence;
}
bool Parameters::GetParticlesState()
{
	return m_DoParticlesCorrespondence;
}

void Parameters::SetUseProcalign(bool _UseProcalign)
{
	m_UseProcalign=_UseProcalign;
}
bool Parameters::GetUseProcalign()
{
	return m_UseProcalign;
}
//Set Data list containing all the data from csv file
void Parameters::SetDataList(vector < vector<string> > _List)
{	
		m_List=_List;
		ListSize = m_List.size();
}

int Parameters::GetListSize()
{	
		return ListSize;
}

int Parameters::MeanTemplateExist()
{
	itksys::Glob glob;
	
	std::vector<std::string> Mean;
	std::string Expression = GetOutputDirectory();
	Expression.insert(Expression.size(),"/Template/*meanAll.vtk");
	glob.RecurseOn();
	glob.FindFiles(Expression);
	Mean = glob.GetFiles();
	
	return Mean.size();
}

//Get data value
string Parameters::GetNthDataListValue(int line,int column)
{
		string file;
		file=m_List[line][column];
		return file;
}

// Set the dimension of the first object to display it in a MRML scene
void Parameters::SetImageDimensions(char *filename)
{
	//read vtk file
	vtkPolyDataReader *meshin = vtkPolyDataReader::New();
	meshin->SetFileName(filename);
	
	try{
		meshin->Update();		
	}
	catch(...)
	{
		std::cout << "Cannot open file: " << filename << " to set dimensions" << std::endl;
		return ;
	}

	vtkPolyData *poly= vtkPolyData::New();
	poly=meshin->GetOutput();
	vtkIdType idNumPointsInFile=poly->GetNumberOfPoints();

	vtkPoints * pts;
	double minCoord[3];
	double maxCoord[3];

	double *firstCoord;
	

	//find the max and min coordinates
	pts=poly->GetPoints();
	firstCoord=pts->GetPoint(0);

	minCoord[0]=firstCoord[0];
	minCoord[1]=firstCoord[1];
	minCoord[2]=firstCoord[2];

	maxCoord[0]=firstCoord[0];
	maxCoord[1]=firstCoord[1];
	maxCoord[2]=firstCoord[2];

  	for(unsigned int i = 1; i < idNumPointsInFile; i++)
   	{	
		double *p;
		p=pts->GetPoint(i);
	
		if(p[0]<=minCoord[0])
		{
			minCoord[0]=p[0];
		}
		
		if(p[1]<=minCoord[1])
		{
			minCoord[1]=p[1];
		}
	
		if(p[2]<=minCoord[2])
		{
			minCoord[2]=p[2];
		}
	
		if(p[0]>=maxCoord[0])
		{
			maxCoord[0]=p[0];
		}
	
		if(p[1]>=maxCoord[1])
		{
			maxCoord[1]=p[1];
		}
	
		if(p[2]>=maxCoord[2])
		{
			maxCoord[2]=p[2];
		}
    	}
	
	m_Dims.clear();

	//find the biggest dimension
	m_Dims.push_back(maxCoord[0]-minCoord[0]);
	m_Dims.push_back(maxCoord[1]-minCoord[1]);
	m_Dims.push_back(maxCoord[2]-minCoord[2]);

	if(m_Dims[0]<m_Dims[1])
	{
		if(m_Dims[1]<m_Dims[2])
			{m_directionToDisplay="ZYX";
			m_const_orientation=m_Dims[0];}
		else
		{
			if(m_Dims[0]<m_Dims[2])
				{m_directionToDisplay="YZX";

				m_const_orientation=m_Dims[0];}
			else 	{m_directionToDisplay="YXZ";
				m_const_orientation=m_Dims[2];}
		}
	}

	else
	{
		if(m_Dims[0]<m_Dims[2])
			{m_directionToDisplay="ZXY";
			m_const_orientation=m_Dims[1];}
		else
		{
			if(m_Dims[1]<m_Dims[2])
				{m_directionToDisplay="XZY";
				m_const_orientation=m_Dims[1];}
			else	 {m_directionToDisplay="XYZ";
				m_const_orientation=m_Dims[2];}
		}			
	}
/*m_Dims.push_back(minCoord[0]);
m_Dims.push_back(minCoord[1]);
m_Dims.push_back(maxCoord[2]);*/





m_Dims.push_back(minCoord[0]);//3
m_Dims.push_back(maxCoord[0]);
m_Dims.push_back(minCoord[1]);//5
m_Dims.push_back(maxCoord[1]);
m_Dims.push_back(minCoord[2]);//7
m_Dims.push_back(maxCoord[2]);

}

vector <double> Parameters::GetImageDimensions()
{
	return m_Dims;
}

double Parameters::GetConstantOrientation()
{
	return m_const_orientation;
}

string Parameters::GetDirectionToDisplay()
{
	return m_directionToDisplay;
}

//Set the names of all the files containninb by the input files 
void Parameters::SetAllFilesName()
{	

	int DataNumber=GetDataNumber();
	m_AllFilesName = new char *[DataNumber];
	//m_ListFiles= new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		m_AllFilesName[i] = new char[512];
		//m_ListFiles[i] = new char[512];
		char c[512];
		char *file;
	
		std::strcpy(c,GetNthDataListValue(i+1,GetColumnVolumeFile()).c_str());
		file=std::strrchr(c,'/');

		for(unsigned int j=0;j<strlen(file);j++)
		{
			file[j]=file[j+1];
			if(file[j]=='.')
				file[j]='\0';
		}
		std::strcpy(m_AllFilesName[i],file);
		//std::strcpy(m_ListFiles[i],file);

	}
}

void Parameters::SetFilesNameMRML(int nummrml)
{
int m_nbHorizontal=GetHorizontalGridPara();
int m_nbVertical=GetVerticalGridPara();
int nbShapesPerMRML= m_nbHorizontal * m_nbVertical;

	int DataNumber=GetDataNumber();
	int DataNumber_csv=GetDataNumber();
	
	SetAllFilesName();

	int begining, end,lastmrml;
	/*begining=nummrml*25;
	end=((nummrml+1)*25)-1;*/
begining=nummrml*nbShapesPerMRML;
	end=((nummrml+1)*nbShapesPerMRML)-1;

	//to know if it's the last mrml
	lastmrml=begining+(end-begining)+1;
	if(DataNumber<=lastmrml)
	{
		end=DataNumber-1;
		DataNumber=end;
	}

	m_ListFiles= new char *[DataNumber];
	
	int listfile_index=0;
/*std::cout<<"end"<<end<<std::endl;
std::cout<<"begining"<<begining<<std::endl;
std::cout<<"DataNumber_csv"<<DataNumber_csv<<std::endl;*/
	for(int i=0;i<DataNumber_csv;i++)
	{
		if(i<=end && i>=begining){

m_ListFiles[listfile_index] = new char[512];

			std::strcpy(m_ListFiles[listfile_index],m_AllFilesName[i]);
			
//std::cout<<"m_ListFiles[listfile_index]"<<m_ListFiles[listfile_index]<<std::endl;
listfile_index++;
		}
	}
//std::cout<<"end"<<std::endl;
}
/*
void Parameters::SetAllFilesName(int nummrml)
{
	int DataNumber_csv=GetDataNumber();
	int DataNumber=GetDataNumber();


	int begining, end,lastmrml;

	if(nummrml>=0){
		
		begining=nummrml*25;
		end=((nummrml+1)*25)-1;
	
		//to know if it's the last mrml
		lastmrml=begining+(end-begining)+1;
		if(DataNumber<=lastmrml)
		{
			end=DataNumber-1;
			DataNumber=end;
		}
	}

	else{
		begining=0;
		end=DataNumber;
	}

	if(nummrml==-1){m_AllFilesName = new char *[DataNumber];}
	else{m_ListFiles= new char *[DataNumber];}
	

	for(int i=0;i<DataNumber_csv;i++)
	{	
		if(nummrml==-1){m_AllFilesName[i] = new char[512];}
		else{m_ListFiles[i] = new char[512];}
		char c[512];
		char *file;
	
		std::strcpy(c,GetNthDataListValue(i+1,GetColumnVolumeFile()).c_str());
		file=std::strrchr(c,'/');

		for(unsigned int j=0;j<strlen(file);j++)
		{
			file[j]=file[j+1];
			if(file[j]=='.')
				file[j]='\0';
		}
	if(nummrml==-1){	std::strcpy(m_AllFilesName[i],file);std::strcpy(m_ListFiles[i],file);}
		else{ 

//std::cout<<"end "<<end<<std::endl;
//std::cout<<"begining "<<begining<<std::endl;

if(i<=end && i>=begining) {
std::cout<<"file "<<file<<std::endl;
std::cout<<"m_ListFiles[i] "<<m_ListFiles[i]<<std::endl;
std::strcpy(m_ListFiles[i],file);}}
	}


}
*/

char* Parameters::GetAllFilesName(int i)
{
	return m_AllFilesName[i];
}

char* Parameters::GetListFiles(int i)
{
	return m_ListFiles[i];
}
string Parameters::GetListFiles_ellalign()
{
	return m_ListFiles_ellalign;
}
string Parameters::GetListFiles_procalign()
{
	return m_ListFiles_procalign;
}


void Parameters::SetOverwriteSegPostProcess(bool overwrite)
{
	m_OverwriteSegPostProcess=overwrite;
}

bool Parameters::GetOverwriteSegPostProcess()
{
	return m_OverwriteSegPostProcess;
}


void Parameters::SetRandomizeInputs(bool random)
{
	m_RandomizeInputs=random;
}

bool Parameters::GetRandomizeInputs()
{
	return m_RandomizeInputs;
}

void Parameters::SetOverwriteGenParaMesh(bool overwrite)
{
	m_OverwriteGenParaMesh=overwrite;
}

bool Parameters::GetOverwriteGenParaMesh()
{
	return m_OverwriteGenParaMesh;
}

void Parameters::SetOverwriteParaToSPHARMMesh(bool overwrite)
{
	m_OverwriteParaToSPHARMMesh=overwrite;
}

bool Parameters::GetOverwriteParaToSPHARMMesh()
{
	return m_OverwriteParaToSPHARMMesh;
}

//To know if the directory is empty or not
bool Parameters::DirectoryIsEmpty(const char * path)
{
	DIR *dir=opendir(path);
	struct dirent *mydir;
	int nbFiles=-2;

	if(dir!=NULL)
	{
		while((mydir=readdir(dir))!=NULL)
		{
			nbFiles++;

		}
	}
	else{ cerr<<" Cannot read the directory "<<path<<endl;}
	if(nbFiles==0)
		return true;
	else
		return false;
}

//Get all the names of the SPHARM files
char* Parameters::GetAllSurfSPHARMFiles(int ind)
{
	int DataNumber=GetDataNumber();
	surfSPHARM_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		surfSPHARM_Files[i] = new char[512];
		std::strcpy(surfSPHARM_Files[i],GetOutputDirectory());
		std::strcat(surfSPHARM_Files[i],"/Mesh/SPHARM/");
		std::strcat(surfSPHARM_Files[i],GetAllFilesName(i));
		std::strcat(surfSPHARM_Files[i],"_pp_surfSPHARM.vtk");
	}
	return surfSPHARM_Files[ind];
}

//Get all the names of the SPHARM ellalign files
char* Parameters::GetAllSurfSPHARMellalignFiles(int ind)
{
	int DataNumber=GetDataNumber();
	surfSPHARM_ellalign_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		surfSPHARM_ellalign_Files[i] = new char[512];
		std::strcpy(surfSPHARM_ellalign_Files[i],GetOutputDirectory());
		std::strcat(surfSPHARM_ellalign_Files[i],"/Mesh/SPHARM/");
		std::strcat(surfSPHARM_ellalign_Files[i],GetAllFilesName(i));
		std::strcat(surfSPHARM_ellalign_Files[i],"_pp_surfSPHARM_ellalign.vtk");

	}
	return surfSPHARM_ellalign_Files[ind];
}

//Get all the names of the SPHARM procalign files
char* Parameters::GetAllSurfSPHARMprocalignFiles(int ind)
{
	int DataNumber=GetDataNumber();
	surfSPHARM_procalign_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		surfSPHARM_procalign_Files[i] = new char[512];
		std::strcpy(surfSPHARM_procalign_Files[i],GetOutputDirectory());
		std::strcat(surfSPHARM_procalign_Files[i],"/Mesh/SPHARM/");
		std::strcat(surfSPHARM_procalign_Files[i],GetAllFilesName(i));
		std::strcat(surfSPHARM_procalign_Files[i],"_pp_surfSPHARM_procalign.vtk");

	}
	return surfSPHARM_procalign_Files[ind];
}



//Get all the names of the SPHARM files
char* Parameters::GetAllSurfmeanSPHARMFiles(int ind)
{
	int DataNumber=GetDataNumber();
	surfSPHARM_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		surfSPHARM_Files[i] = new char[512];
		std::strcpy(surfSPHARM_Files[i],GetOutputDirectory());
		std::strcat(surfSPHARM_Files[i],"/Mesh/SPHARM/");
		std::strcat(surfSPHARM_Files[i],GetAllFilesName(i));
		std::strcat(surfSPHARM_Files[i],"_pp_surf_tMeanSPHARM.vtk");
	}


	return surfSPHARM_Files[ind];
}

//Get all the names of the SPHARM ellalign files
char* Parameters::GetAllSurfmeanSPHARMellalignFiles(int ind)
{
	int DataNumber=GetDataNumber();
	surfSPHARM_ellalign_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		surfSPHARM_ellalign_Files[i] = new char[512];
		std::strcpy(surfSPHARM_ellalign_Files[i],GetOutputDirectory());
		std::strcat(surfSPHARM_ellalign_Files[i],"/Mesh/SPHARM/");
		std::strcat(surfSPHARM_ellalign_Files[i],GetAllFilesName(i));
		std::strcat(surfSPHARM_ellalign_Files[i],"_pp_surf_tMeanSPHARM_ellalign.vtk");
	}
	return surfSPHARM_ellalign_Files[ind];
}

//Get all the names of the SPHARM procalign files
char* Parameters::GetAllSurfmeanSPHARMprocalignFiles(int ind)
{
	int DataNumber=GetDataNumber();
	surfSPHARM_procalign_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		surfSPHARM_procalign_Files[i] = new char[512];
		std::strcpy(surfSPHARM_procalign_Files[i],GetOutputDirectory());
		std::strcat(surfSPHARM_procalign_Files[i],"/Mesh/SPHARM/");
		std::strcat(surfSPHARM_procalign_Files[i],GetAllFilesName(i));
		std::strcat(surfSPHARM_procalign_Files[i],"_pp_surf_tMeanSPHARM_procalign.vtk");
	}
	return surfSPHARM_procalign_Files[ind];
}



//Get all the names of the phi files
char* Parameters::GetAllPhiFiles(int ind)
{
	Phi_Files = new char[512];
	char c[512];
	char *file;
	
	std::strcpy(c,GetNthDataListValue(1,GetColumnVolumeFile()).c_str());

	file=std::strrchr(c,'/');

	for(unsigned int j=0;j<strlen(file);j++)
	{
		file[j]=file[j+1];
		if(file[j]=='.')
			file[j]='\0';
	}
	std::strcpy(Phi_Files,GetOutputDirectory());
	std::strcat(Phi_Files,"/Mesh/SPHARM/");
	std::strcat(Phi_Files,file);

	std::strcat(Phi_Files,"_pp_surf_paraPhi.txt");
	return Phi_Files;
}

//Get all the names of the theta files
char* Parameters::GetAllThetaFiles(int ind)
{
	Theta_Files = new char[512];
	char c[512];
	char *file;
	
	std::strcpy(c,GetNthDataListValue(1,GetColumnVolumeFile()).c_str());
	file=std::strrchr(c,'/');

	for(unsigned int j=0;j<strlen(file);j++)
	{
		file[j]=file[j+1];
		if(file[j]=='.')
			file[j]='\0';
	}
	std::strcpy(Theta_Files,GetOutputDirectory());
	std::strcat(Theta_Files,"/Mesh/SPHARM/");
	std::strcat(Theta_Files,file);
	std::strcat(Theta_Files,"_pp_surf_paraTheta.txt");
	return Theta_Files;
}

void Parameters::SetEulerNumbers()
{
	for(unsigned int  i=0;i<EulerFile.size();i++)
	{
		ifstream file(EulerFile[i].c_str(), ios::in); 
		std::string eulerNumber_string;  
		getline(file, eulerNumber_string);
		m_EulerNumber.push_back (eulerNumber_string);
		if (eulerNumber_string.at(0)=='2') {m_SphericalTopo.push_back ("1");}
				else {m_SphericalTopo.push_back ("0");}	
		file.close();
	}

}

void Parameters::ModifyCSV(int Particles)
{
	FindFiles();
	SetEulerNumbers();
	
	//GetPostCorrespondenceFiles();

	char csv_read[512];
	std::strcpy(csv_read,GetOutputDirectory());
	std::strcat(csv_read,"/OutputGroupFile/");
	std::strcat(csv_read,"ShapeAnalysisModule_OutputFileVersion1.csv");

	char csv_write[512];
	std::strcpy(csv_write,GetOutputDirectory());
	std::strcat(csv_write,"/OutputGroupFile/");
	std::strcat(csv_write,"ShapeAnalysisModule_OutputFile.csv");

	ifstream read(csv_read);
	ofstream write(csv_write, ios::out);

	int nbline=0;
	int in =0;

	if(read){
		if(write)
		{
			std::string line;
			while(getline(read, line))
			{
				if(nbline==0)
				{
					line.append( ", Euler Number, Spherical topology ,Correspondence ");
					write<<line<<std::endl;
				}
				else{
					if( nbline< m_EulerNumber.size())
					{					
						line.append(",");
						line.append(m_EulerNumber.at(nbline-1));
						line.append(",");
						line.append(m_SphericalTopo.at(nbline-1));
					//	in=1;
					//}
					if(Particles ==1 )
					{					
						line.append(",");
						line.append(GetPostCorrespondenceFiles(nbline-1));
std::cout<<line<<std::endl;

						//in=1;
					}
					/*if( in == 1) {*/write<<line<<std::endl; /*in =0;*/}
				}
			  nbline++;
			}
			write.close();
			read.close();
		}

	}
	remove( csv_read );
}

char* Parameters::GetPostCorrespondenceFiles(int ind)
{
	int DataNumber=GetDataNumber();
	Corres_Files = new char *[DataNumber];

	for(int i=0;i<DataNumber;i++)
	{	
		Corres_Files[i] = new char[512];
		std::strcpy(Corres_Files[i],GetOutputDirectory());
		std::strcat(Corres_Files[i],"/ParticleCorrespondence/Corresponding_Meshes/");
		std::strcat(Corres_Files[i],GetAllFilesName(i));
		if(GetUseProcalign()){std::strcat(Corres_Files[i],"_pp_surfSPHARM_procalign_corr.vtk");}
		else{std::strcat(Corres_Files[i],"_pp_surfSPHARM_corr.vtk");}
	}
	
	return Corres_Files[ind];
}



std::string Parameters::readMRML(std::string nameMRML, bool colorMap)
{
	std::string textMRML;


	ifstream read(nameMRML.c_str());
	if(read){
			std::string line;
			while(getline(read, line))
			{
				textMRML.append(line);
				textMRML.append("\n");
			}
			read.close();
		}


	if(colorMap ==0){ //delete the last line </ MRML >
		textMRML.erase (textMRML.end()-9, textMRML.end());}

	else{ //delete the last and first line <MRML >
		textMRML.erase (textMRML.begin(), textMRML.begin()+7);
		textMRML.erase (textMRML.end()-9, textMRML.end());}

	return textMRML;

}

void Parameters::ModifyMRML(std::string nameMRML,std::string nameMRMLPhi, std::string nameMRMLTheta)
{
	std::string help;
		
	std::string textMRML;
	std::string textMRMLPhi;
	std::string textMRMLTheta;
		
	help.append(readMRML(nameMRML,0));
	
	std::string  mrml_write;

	mrml_write.append(nameMRML);

	ofstream filemrml(mrml_write.c_str(), ios::out ); 

		if(filemrml)
		{
			
			filemrml << help <<"\n";
			help.clear();
		
			//add SceneSnapshot phi
			filemrml <<"<SceneSnapshot\n id=\"vtkMRMLSceneSnapshotNode1\" name=\"Color Map Phi\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\">\n <Selection\n id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"NULL\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>\n <Interaction\n id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>\n <Layout\n id=\"vtkMRMLLayoutNode1\" name=\"vtkMRMLLayoutNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentViewArrangement=\"2\" guiPanelVisibility=\"1\" bottomPanelVisibility =\"1\" guiPanelLR=\"0\" numberOfCompareViewRows=\"0\" numberOfCompareViewColumns=\"0\" numberOfLightboxRows=\"1\" numberOfLightboxColumns=\"1\" mainPanelSize=\"400\" secondaryPanelSize=\"400\"></Layout>\n <View\n id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"true\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"true\" backgroundColor=\"0.70196 0.70196 0.90588\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>\n <Camera\n id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" position=\"-726.497 88.984 13.4559\" focalPoint=\"0 0 0\" viewUp=\"0 0 1\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>\n <TGParameters\n id=\"vtkMRMLChangeTrackerNode1\" name=\"vtkMRMLChangeTrackerNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ROIMin=\"-1 -1 -1\" ROIMax=\"-1 -1 -1\" SegmentThresholdMin=\"-1\" SegmentThresholdMax=\"-1\" Analysis_Intensity_Flag=\"0\" Analysis_Deformable_Flag=\"0\" UseITK=\"1\"></TGParameters> <VolumeRenderingSelection\n id=\"vtkMRMLVolumeRenderingSelectionNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" activeVolumeRenderingID=\"NULL\"></VolumeRenderingSelection>\n <Slice\n id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"387.5 250 1\" dimensions=\"496 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n <SliceComposite\n id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\"\n doPropagateVolumeSelection=\"1\"></SliceComposite>\n <Slice\n id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n <SliceComposite\n id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n <Slice\n id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"0 0 1 0 -1 0 0 0 0 1 0 0 0 0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n <SliceComposite\n id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n <Crosshair\n id=\"vtkMRMLCrosshairNode1\" name=\"vtkMRMLCrosshairNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" crosshairMode=\"NoCrosshair\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" crosshairRAS=\"0 0 0\"></Crosshair>\n <ClipModels\n id=\"vtkMRMLClipModelsNode1\" name=\"vtkMRMLClipModelsNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" clipType=\"0\" redSliceClipState=\"0\" yellowSliceClipState=\"0\" greenSliceClipState=\"0\"></ClipModels>\n <ScriptedModule\n id=\"vtkMRMLScriptedModuleNode1\" name=\"vtkMRMLScriptedModuleNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\"></ScriptedModule>\n ";

			help.append(readMRML(nameMRMLPhi,1));
			filemrml << help <<"\n";
			help.clear();
			filemrml <<"</SceneSnapshot>"<<"\n";

			//add SceneSnapshot theta
			filemrml <<"<SceneSnapshot\n id=\"vtkMRMLSceneSnapshotNode1\" name=\"Color Map Theta\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\">\n <Selection\n id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"NULL\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>\n <Interaction\n id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>\n <Layout\n id=\"vtkMRMLLayoutNode1\" name=\"vtkMRMLLayoutNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentViewArrangement=\"2\" guiPanelVisibility=\"1\" bottomPanelVisibility =\"1\" guiPanelLR=\"0\" numberOfCompareViewRows=\"0\" numberOfCompareViewColumns=\"0\" numberOfLightboxRows=\"1\" numberOfLightboxColumns=\"1\" mainPanelSize=\"400\" secondaryPanelSize=\"400\"></Layout>\n <View\n id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"true\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"true\" backgroundColor=\"0.70196 0.70196 0.90588\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>\n <Camera\n id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" position=\"-726.497 88.984 13.4559\" focalPoint=\"0 0 0\" viewUp=\"0 0 1\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>\n <TGParameters\n id=\"vtkMRMLChangeTrackerNode1\" name=\"vtkMRMLChangeTrackerNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ROIMin=\"-1 -1 -1\" ROIMax=\"-1 -1 -1\" SegmentThresholdMin=\"-1\" SegmentThresholdMax=\"-1\" Analysis_Intensity_Flag=\"0\" Analysis_Deformable_Flag=\"0\" UseITK=\"1\"></TGParameters> <VolumeRenderingSelection\n id=\"vtkMRMLVolumeRenderingSelectionNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" activeVolumeRenderingID=\"NULL\"></VolumeRenderingSelection>\n <Slice\n id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"387.5 250 1\" dimensions=\"496 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n <SliceComposite\n id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\"\n doPropagateVolumeSelection=\"1\"></SliceComposite>\n <Slice\n id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n <SliceComposite\n id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n <Slice\n id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"0 0 1 0 -1 0 0 0 0 1 0 0 0 0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n <SliceComposite\n id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\n <Crosshair\n id=\"vtkMRMLCrosshairNode1\" name=\"vtkMRMLCrosshairNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" crosshairMode=\"NoCrosshair\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" crosshairRAS=\"0 0 0\"></Crosshair>\n <ClipModels\n id=\"vtkMRMLClipModelsNode1\" name=\"vtkMRMLClipModelsNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" clipType=\"0\" redSliceClipState=\"0\" yellowSliceClipState=\"0\" greenSliceClipState=\"0\"></ClipModels>\n <ScriptedModule\n id=\"vtkMRMLScriptedModuleNode1\" name=\"vtkMRMLScriptedModuleNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\"></ScriptedModule>\n ";

			help.append(readMRML(nameMRMLTheta,1));
			filemrml << help <<"\n";
			help.clear();
			filemrml <<"</SceneSnapshot>"<<"\n"<<"</MRML>";
			
			filemrml.close();
	
	}

	//remove( nameMRMLPhi.c_str());
	//remove( nameMRMLTheta.c_str() );
}



void Parameters::FindFiles()
{
	itksys::Glob globEulerFile;
	std::string pathFile =GetOutputDirectory();
	std::string path="/EulerFiles/*.txt";
	pathFile=pathFile+path;
	globEulerFile.FindFiles(pathFile);
	EulerFile=globEulerFile.GetFiles();
}

void Parameters::FindTemplateFiles(int type)
{
	itksys::Glob globTemplateFile;
	std::string pathFile =GetOutputDirectory();
	std::string path;

	if(GetTemplateMState()){
		path="/Template/meanAll.vtk";
		
	}
	else{
		if(type==0  ||type==2){path="/Template/*_pp_surfSPHARM.vtk";}
		if(type==1){ path="/Template/*_pp_surfSPHARM_ellalign.vtk";}
		//if(type==2){path="/Template/*_pp_surfSPHARM_procalign.vtk";}
	}
	pathFile=pathFile+path;
	globTemplateFile.FindFiles(pathFile);
	name_template=globTemplateFile.GetFiles();

}

std::string Parameters::GetTemplate( int type)
{

	FindTemplateFiles(type);
size_t found;
found=name_template[0].find("/Template");
std::string mytemplate;
std::string help;
help.append(name_template[0]);
//mytemplate.append("..")
//name_template[0]).erase(0,found);
if(type==1){mytemplate.append("../");}
if(type==2){mytemplate.append("../");}
mytemplate.append("..");
mytemplate.append(help.begin()+found,help.end());
	return mytemplate ;
}

int Parameters::SetNbSnapShot()
{
	int DataNumber=GetDataNumber();
	if(DataNumber>25)
	{int SnapShotNumber= (DataNumber-1)/24+1;
	return SnapShotNumber;}
	else
	{return 0;}
}

void Parameters::DeleteTransformsFolders(int type)
{
	char *dirTransform=NULL;
	dirTransform=new  char[512] ;
	std::strcpy(dirTransform, GetOutputDirectory());
	if(type==0){std::strcat(dirTransform, "/MRML/TransformFiles/");}
	if(type==1){std::strcat(dirTransform, "/MRML/Ellalign/TransformFiles/");}
	if(type==2){std::strcat(dirTransform, "/MRML/Procalign/TransformFiles/");}
	
	bool transformDirectoryEmpty=DirectoryIsEmpty(dirTransform);
	if(!transformDirectoryEmpty)
	{
		int length=strlen(GetOutputDirectory());
		length=length+10;

		DIR *pdir = NULL;
		struct dirent *pent;
		pdir = opendir (dirTransform);
		while ((pent=readdir(pdir)))
		{
			char *file=NULL;
			file= new char[512];
			strcpy(file,dirTransform);
			strcat(file, pent->d_name);
			if(file[length]!='.')
			{	
				remove(file);
				//cerr<<"Error deleting file "<<file<<endl;
			}
		}
	}


}

char * Parameters::Convert_Double_To_CharArray(double doubleVariable) 
{
	char *CharArrayOut;
	CharArrayOut = new char [512];
	std::string stringVariable;
	std::stringstream strStream;
	strStream << doubleVariable;
	stringVariable = strStream.str();
	strcpy(CharArrayOut,stringVariable.c_str());
	return CharArrayOut;
}
