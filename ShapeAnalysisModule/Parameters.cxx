 #include "Parameters.h"

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
						char *volume=strchr(value.c_str(),'.');
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
								char *volume=strchr(value.c_str(),'.');
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

void Parameters::SetTemplateState(bool _TemplateState)
{
  	m_TemplateState = _TemplateState;
}

bool Parameters::GetTemplateState()
{
  	return m_TemplateState;
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

//Set Data list containing all the data from csv file
void Parameters::SetDataList(vector < vector<string> > _List)
{	
		m_List=_List;
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
			m_directionToDisplay="ZYX";
		else
		{
			if(m_Dims[0]<m_Dims[2])
				m_directionToDisplay="YZX";
			else 	m_directionToDisplay="YXZ";
		}
	}

	else
	{
		if(m_Dims[0]<m_Dims[2])
			m_directionToDisplay="ZXY";
		else
		{
			if(m_Dims[1]<m_Dims[2])
				m_directionToDisplay="XZY";
			else	 m_directionToDisplay="XYZ";
		}			
	}
}

vector <double> Parameters::GetImageDimensions()
{
	return m_Dims;
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

	for(int i=0;i<DataNumber;i++)
	{	
		m_AllFilesName[i] = new char[512];
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
	}
}

char* Parameters::GetAllFilesName(int i)
{
	return m_AllFilesName[i];
}

void Parameters::SetOverwriteSegPostProcess(bool overwrite)
{
	m_OverwriteSegPostProcess=overwrite;
}

bool Parameters::GetOverwriteSegPostProcess()
{
	return m_OverwriteSegPostProcess;
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
	std::strcat(Phi_Files,"/Template/");
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
	std::strcat(Theta_Files,"/Template/");
	std::strcat(Theta_Files,file);
	std::strcat(Theta_Files,"_pp_surf_paraTheta.txt");
	return Theta_Files;
}
