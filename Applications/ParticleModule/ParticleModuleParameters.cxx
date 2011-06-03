#include "ParticleModuleParameters.h"

using namespace std;

 ParticleModuleParameters::ParticleModuleParameters()
 {
 }

ParticleModuleParameters::~ParticleModuleParameters()
 {}



//.....................................................................................
//                           Parameters from the command line
//....................................................................................

void ParticleModuleParameters::SetParameterFile(const char *_ParamFile)
{
   	std::strcpy(m_ParamFile, _ParamFile);
}

const char* ParticleModuleParameters::GetParameterFile()
{
	return m_ParamFile;
}
void ParticleModuleParameters::SetOutputDirectory(const char *_OutputDirectory)
{
 	 std::strcpy(m_OutputDirectory, _OutputDirectory);
}

char* ParticleModuleParameters::GetOutputDirectory()
{
  	return m_OutputDirectory;
}


void ParticleModuleParameters::SetColumnMeshFile(int column)
{	
	 m_ColumnMesh=column;
}

int ParticleModuleParameters::GetColumnMeshFile()
{
	return m_ColumnMesh;
}


void ParticleModuleParameters::SetSmooting(double smooting)
{	
	 m_smoothing=smooting;
}

double ParticleModuleParameters::GetSmoomting()
{
	return m_smoothing;
}

void ParticleModuleParameters::SetEnforcedSpace(float _sx, float _sy, float _sz)
{
	if (_sx >= 0.0)
	{m_sx = _sx;}
	
	if (_sy >= 0.0)
	{m_sy = _sy;}
	
	if (_sz >= 0.0)
	{m_sz = _sz;}
}

float ParticleModuleParameters::GetEnforcedSpaceX()
{  	return m_sx;
}

float ParticleModuleParameters::GetEnforcedSpaceY()
{
  	return m_sy;
}

float ParticleModuleParameters::GetEnforcedSpaceZ()
{
  	return m_sz;
}

void ParticleModuleParameters::SetEndRegularization(double endRegularization)
{	
	 m_endRegularization=endRegularization;
}

double ParticleModuleParameters::GetEndRegularization()
{
	return m_endRegularization;
}

void ParticleModuleParameters::SetStartRegularization(double startRegularization)
{	
	 m_startRegularization=startRegularization;
}

double ParticleModuleParameters::GetStartRegularization()
{
	return m_startRegularization;
}

void ParticleModuleParameters::SetOptimizationIteration(int optimizationIteration)
{	
	 m_optimizationIteration=optimizationIteration;
}

int ParticleModuleParameters::GetOptimizationIteration()
{
	return m_optimizationIteration;
}

void ParticleModuleParameters::SetCheckpointingInterval(int checkpointingInterval)
{	
	 m_checkpointingInterval=checkpointingInterval;
}

int ParticleModuleParameters::GetCheckpointingInterval()
{
	return m_checkpointingInterval;
}


void ParticleModuleParameters::SetRelativeWeighting(float RelativeWeighting)
{	
	 m_RelativeWeight=RelativeWeighting;
}

float ParticleModuleParameters::GetRelativeWeighting()
{
  	return m_RelativeWeight;
}

void ParticleModuleParameters::SetAdaptivityStrength(float AdaptivityStrength)
{	
	 m_AdaptivityStrength=AdaptivityStrength;
}

float ParticleModuleParameters::GetAdaptivityStrength()
{
  	return m_AdaptivityStrength;
}

void ParticleModuleParameters::SetProcrustesInterval(int ProcrustesInterval)
{	
	 m_ProcrustesInterval=ProcrustesInterval;
}

int ParticleModuleParameters::GetProcrustesInterval()
{
	return m_ProcrustesInterval;
}

void ParticleModuleParameters::SetProcrustesScaling(int ProcrustesScaling)
{	
	 m_ProcrustesScaling=ProcrustesScaling;
}

int ParticleModuleParameters::GetProcrustesScaling()
{
	return m_ProcrustesScaling;
}

void ParticleModuleParameters::SetRelativeWeightingState(bool RelativeWeightingState)
{
  	m_RelativeWeightingState = RelativeWeightingState;
}

bool ParticleModuleParameters::GetRelativeWeightingState()
{
  	return m_RelativeWeightingState;
}
void ParticleModuleParameters::SetAdaptivityStrengthState(bool AdaptivityStrengthState)
{
  	m_AdaptivityStrengthState = AdaptivityStrengthState;
}

bool ParticleModuleParameters::GetAdaptivityStrengthState()
{
  	return m_AdaptivityStrengthState;
}
void ParticleModuleParameters::SetProcrustesIntervalState(bool ProcrustesIntervalState)
{
  	m_ProcrustesIntervalState = ProcrustesIntervalState;
}

bool ParticleModuleParameters::GetProcrustesIntervalState()
{
  	return m_ProcrustesIntervalState;
}
void ParticleModuleParameters::SetProcrustesScalingState(bool ProcrustesScalingState)
{
  	m_ProcrustesScalingState = ProcrustesScalingState;
}

bool ParticleModuleParameters::GetProcrustesScalingState()
{
  	return m_ProcrustesScalingState;
}

void ParticleModuleParameters::SetHorizontalGridPara(int _HorizontalGridPara)
{
	m_HorizontalGridParaGrid = _HorizontalGridPara;
}
void ParticleModuleParameters::SetVerticalGridPara(int _VerticalGridPara)
{
	 m_VerticalGridParaGrid = _VerticalGridPara;
}

int ParticleModuleParameters::GetHorizontalGridPara()
{
	return	m_HorizontalGridParaGrid;
}
int ParticleModuleParameters::GetVerticalGridPara()
{
	return m_VerticalGridParaGrid;
}

void ParticleModuleParameters::WriteCommandLine()
{
	char fileCommanLine[512];
	std::strcpy (fileCommanLine,GetOutputDirectory());
	std::strcat(fileCommanLine,"/commandline.txt");
	
	std::ofstream file(fileCommanLine);

	if(file)
	{
		file<<"ParticleModule ";
		file<<"--columMeshFile "<<m_ColumnMesh;
		file<<" --smoothing "<<GetSmoomting();
		file<<" --startRegularization "<<GetStartRegularization();
		file<<" --endRegularization "<<GetEndRegularization();
		file<<" --optimizationIt "<<GetOptimizationIteration();
		file<<" --checkpointInter "<<GetCheckpointingInterval();

		if(GetProcrustesScalingState()){
			file<<" --RelativeWeight "<<Convert_Double_To_CharArray(GetRelativeWeighting());
		}
		if(GetProcrustesIntervalState()){
			file<<" --AdaptivityStrength "<<Convert_Double_To_CharArray(GetAdaptivityStrength());
		}
		if(GetAdaptivityStrengthState()){
			file<<" --ProcrustesInt "<<Convert_Double_To_CharArray(GetProcrustesInterval());
		}
		if(GetRelativeWeightingState()){
			file<<" --ProcrustesScal "<<Convert_Double_To_CharArray(GetCheckpointingInterval());
		}
		file<<" ";
		file<<m_ParamFile<<" ";
		file<<GetOutputDirectory();

		file.close();
	}
}



//.....................................................................................
//                                Tools
//....................................................................................


char * ParticleModuleParameters::Convert_Double_To_CharArray(double doubleVariable) {
	char *CharArrayOut;
	CharArrayOut = new char [512];
	std::string stringVariable;
	std::stringstream strStream;
	strStream << doubleVariable;
	stringVariable = strStream.str();
	strcpy(CharArrayOut,stringVariable.c_str());

	return CharArrayOut;
}

double ParticleModuleParameters::string_to_double( const std::string& s )
 {
   std::istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 } 

void ParticleModuleParameters::Run(std::vector<const char*> args, bool TimeOn)
{		
	//itk sys parameters
	int length;
	time_t start,end;
	time (&start);

	double timeout = 0.05;
	int result;
	char* dataitk = NULL;

	itksysProcess* gp = itksysProcess_New();
	itksysProcess_SetCommand(gp, &*args.begin());
	itksysProcess_SetOption(gp,itksysProcess_Option_HideWindow,1);
	itksysProcess_Execute(gp);
	while(int Value = itksysProcess_WaitForData(gp,&dataitk,&length,&timeout)) // wait for 1s
	{
		if ( ((Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR)) && dataitk[0]=='D' )
		{
			std::strstream st;
			for(int i=0;i<length;i++) 	
			{
				st<<dataitk[i];
			}
			std::string dim=st.str();
		}
			if(TimeOn){
				time (&end);
				cout<<"(processing since "<<difftime (end,start)<<" seconds) \r"<<flush;
				timeout = 0.05; 
			}  	
	}
	itksysProcess_WaitForExit(gp, 0);
	result = 1;
	switch(itksysProcess_GetState(gp))
	{
		case itksysProcess_State_Exited:
		{
			result = itksysProcess_GetExitValue(gp);
		} break;
		case itksysProcess_State_Error:
		{
			std::cerr<<"Error: Could not run " << args[0]<<":\n";
			std::cerr<<itksysProcess_GetErrorString(gp)<<"\n";
			std::cout<<"Error: Could not run " << args[0]<<":\n";
			std::cout<<itksysProcess_GetErrorString(gp)<<"\n";
		} break;
		case itksysProcess_State_Exception:
		{
			std::cerr<<"Error: "<<args[0]<<" terminated with an exception: "<<itksysProcess_GetExceptionString(gp)<<"\n";
			std::cout<<"Error: "<<args[0]<<" terminated with an exception: "<<itksysProcess_GetExceptionString(gp)<<"\n";
		} break;
		case itksysProcess_State_Starting:
		case itksysProcess_State_Executing:
		case itksysProcess_State_Expired:
		case itksysProcess_State_Killed:
		{
		// Should not get here.
		std::cerr<<"Unexpected ending state after running "<<args[0]<<std::endl;
		std::cout<<"Unexpected ending state after running "<<args[0]<<std::endl;
		} break;
	}
	itksysProcess_Delete(gp);  
}

//.....................................................................................
//                                Initializations
//....................................................................................


//create the sub directories
void ParticleModuleParameters::OrganizeOutputDirectory()
{
	std::string lptsbefore,lptsafter,preprocess,corresp,meta,mrml,trans,meshscale;
	lptsbefore.append(GetOutputDirectory());
	lptsafter.append(GetOutputDirectory());
	preprocess.append(GetOutputDirectory());
	corresp.append(GetOutputDirectory());
	meta.append(GetOutputDirectory());
	mrml.append(GetOutputDirectory());
	trans.append(GetOutputDirectory());
	meshscale.append(GetOutputDirectory());
	lptsbefore.append("/Initialization_Particles");
	lptsafter.append("/Corresponding_Particles");
	preprocess.append("/PreProcessing");
	corresp.append("/Corresponding_Meshes");
	meta.append("/InputMeshes_scale");
	mrml.append("/MRML");
	trans.append("/MRML/TransformFiles");
	meshscale.append("/Corresponding_Meshes_scale");

	itksys::SystemTools::MakeDirectory(lptsbefore.c_str());
	itksys::SystemTools::MakeDirectory(lptsafter.c_str());
	itksys::SystemTools::MakeDirectory(preprocess.c_str());
	itksys::SystemTools::MakeDirectory(corresp.c_str());
	itksys::SystemTools::MakeDirectory(meta.c_str());
	itksys::SystemTools::MakeDirectory(mrml.c_str());
	itksys::SystemTools::MakeDirectory(trans.c_str());
	itksys::SystemTools::MakeDirectory(meshscale.c_str());
}

//read the input file to know the .vtk usefull for this running
int ParticleModuleParameters::ReadFile(const char *_FileName)
{
	int column=0;
	int nbLine=0;
	int problem =0;
	vector <vector<string> > List;
	size_t found,found2;

	//initialisation;
	nbPoints=0;
	
	std::ifstream file(_FileName, std::ios::in);
	string Line;
	if(file)
	{
		while(getline(file,Line))
		{
			istringstream s(Line);
			string value,value2;

			if(!Line.empty())
			{
				while(getline(s,value,','))
					{

					if(nbLine!=0){
						if(m_ColumnMesh==column) 
						{
							found=value.find_last_of(".");
							found2=value.find_last_of("/");
							//name of the vtk in the Data vector
							Data.push_back(value);
							value2.append(value.substr(found2+1));
							value2.resize(value2.size()-4);
							//name of the files witout exentions in Data2
							Data2.push_back(value2);
						}
					}

					column++;	
				}
			}		column=0;
		nbLine++;
		}
	}

 	else{
                std::cout << "problem with the input file" << std::endl;
		problem=1;}

	return problem;
		
}


//.....................................................................................
//                                Scale Factor
//....................................................................................


void ParticleModuleParameters::ScaleMeta(std::string vtkinput,std::string vtkoutput , float factorX, float factorY, float factorZ)
{

	std::vector<const char*> args2;
	std::string factor;
	factor.append(Convert_Double_To_CharArray(factorX));
	factor.append(",");
	factor.append(Convert_Double_To_CharArray(factorY));
	factor.append(",");
	factor.append(Convert_Double_To_CharArray(factorZ));


	args2.push_back("MeshMath");
	args2.push_back(vtkinput.c_str());
	args2.push_back(vtkoutput.c_str());
	args2.push_back("-scaleMesh");
	args2.push_back(factor.c_str());

	args2.push_back(0);

	Run(args2,0);
}


//.....................................................................................
//                                Create LPTS files
//....................................................................................


void ParticleModuleParameters::SetLptsFilesName(std::string VTKName)
{
	size_t found;
	std::string lptsName,lptsNameFinal;
	found=VTKName.find_last_of("/");
	lptsNameFinal= VTKName.substr (found);
	lptsNameFinal.insert(0,"/Initialization_Particles");
	lptsNameFinal.insert(0,GetOutputDirectory());
	lptsNameFinal.erase(lptsNameFinal.end()-3,lptsNameFinal.end());
	lptsNameFinal.append("lpts");
	LPTS.push_back(lptsNameFinal);
}

void ParticleModuleParameters::SetVTKFilesName(std::string VTKName)
{
	std::string vtkoutput;
	vtkoutput.append(GetOutputDirectory());
	vtkoutput.append("/InputMeshes_scale/");
	vtkoutput.append( VTKName);
	vtkoutput.append("_scale.vtk");
	VTKPrePross.push_back(vtkoutput);
}

//create the .lpts files from the vtk files
void ParticleModuleParameters::CreateLptsFiles()
{
	
	for(unsigned int  i=0;i<Data.size();i++)
	{	
		SetVTKFilesName(Data2[i]);
		SetLptsFilesName(VTKPrePross.back());
		ScaleMeta(Data[i],VTKPrePross.back(),1/ GetEnforcedSpaceX(),1/ GetEnforcedSpaceY(),1/ GetEnforcedSpaceZ());
		VTKToLpts(VTKPrePross.back(),LPTS.back() );	
	}
	std::cout<<" "<<std::endl;

}

void ParticleModuleParameters::VTKToLpts(std::string VTKName, std::string LptsName)  
{
	int Pointwritten;
	Pointwritten=0;

	if(nbPoints ==0)
	{
		vtkPolyDataReader *meshin = vtkPolyDataReader::New();
		meshin->SetFileName(VTKName.c_str());
		meshin->Update();
		nbPoints=meshin->GetOutput()->GetNumberOfPoints();
	}

	std::ofstream lptsfile(LptsName.c_str(), std::ios::out | std::ios::trunc);

	vtkPolyDataReader *vtkreader = vtkPolyDataReader::New();
	vtkreader->SetFileName(VTKName.c_str());
	vtkreader->Update();
	double x[3];
	vtkPoints * PointVTK = vtkPoints::New();
	PointVTK=vtkreader->GetOutput()->GetPoints();
	for (int PointId = 0; PointId < (vtkreader->GetOutput()->GetNumberOfPoints()); PointId++)
	{
		PointVTK->GetPoint(PointId,x);
		lptsfile<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
	}
	lptsfile.close();
}


//.....................................................................................
//                                PreProcessing
//....................................................................................

void ParticleModuleParameters::CreatePreprocessFiles()
{
	std::cout<<"Creation of the CorrespondencePreprocessing.params file: "<<std::endl;
	preprocess_file.append(GetOutputDirectory());
	preprocess_file.append("/PreProcessing/CorrespondencePreprocessing.params");
	std::ofstream write(preprocess_file.c_str(), std::ios::out);
	if(write){

			write<<"NUMBER_OF_SHAPES="<<VTKPrePross.size()<<std::endl;
			write<<"NUMBER_OF_ATTRIBUTES=0"<<std::endl;
			write<<"TEMPLATE_INDEX=1"<<std::endl;
			write<<"WEIGHTS= "<<std::endl;
			for( unsigned int i=0; i<VTKPrePross.size();i++)
			{
				write<<VTKPrePross.at(i)<<std::endl;  
			}
			write.close();
	}
	std::cout<<preprocess_file<<std::endl;
	std::cout<<" "<<std::endl;
}

void ParticleModuleParameters::MakePPDirectory()
{
	PPOutputDirectory.append(GetOutputDirectory() );
	PPOutputDirectory.append("/PreProcessing");
}


void ParticleModuleParameters::CommandlineParticleCorrespondencePreprocessing()
{

	MakePPDirectory();

	args.push_back("ParticleCorrespondencePreprocessing" );
	args.push_back("--parameterFileName" );
	args.push_back(preprocess_file.c_str() );
	args.push_back("--voxelSize" );
	args.push_back("1" );
	args.push_back("--smoothing" );
	args.push_back(Convert_Double_To_CharArray(GetSmoomting()));
	args.push_back("--outputDirectory" );
	args.push_back(PPOutputDirectory.c_str());
	args.push_back(0);

	std::cout<<"ParticleCorrespondencePreprocessing command line: "<<std::endl;
	for( unsigned int k =0; k<args.size()-1;k++)
	{std::cout<<args.at(k)<<" ";}
	std::cout<<" "<<std::endl;std::cout<<" "<<std::endl;

}

void ParticleModuleParameters::RunParticleCorrespondencePreprocessing()
{
	CommandlineParticleCorrespondencePreprocessing();
	
	
	std::cout<<"....Run ParticleCorrespondencePreprocessing.... "<<std::endl;
	Run(args,1);	

}

//.....................................................................................
//                                Correspondence
//....................................................................................

//find all the .mha and stock there names
void ParticleModuleParameters::FindMhaFiles()
{
	itksys::Glob globMHAFile;
	std::string pathmha="/PreProcessing/*DistanceMap.mha";
	std::string pathFile;
	pathFile=GetOutputDirectory()+pathmha;
	globMHAFile.FindFiles(pathFile);
	mhaFile=globMHAFile.GetFiles();
}

void ParticleModuleParameters::CreateCorrespondenceFiles()
{
	std::cout<<" "<<std::endl;std::cout<<" "<<std::endl;
	std::cout<<"Creation of the Correspondence.params files: "<<std::endl;

	FindMhaFiles();

	Correspondence_file.append(GetOutputDirectory());

	Correspondence_file.append("/Corresponding_Particles/Correspondence.params");
	std::ofstream writemha(Correspondence_file.c_str(), std::ios::out);
	if(writemha){
;
			writemha<<"(inputs"<<std::endl;
			for( unsigned int m=0; m<mhaFile.size();m++)
			{
				writemha<<"\""<<mhaFile[m]<<"\""<<std::endl;
			}
			writemha<<")"<<std::endl;

			writemha<<"(point_files "<<std::endl;
			for( unsigned int i=0; i<LPTS.size();i++)
			{
				writemha<<"\""<<LPTS.at(i)<<"\""<<std::endl;
			}
			writemha<<")"<<std::endl;

			writemha<<"(number_of_particles ";
			writemha<<Convert_Double_To_CharArray(nbPoints);  
			writemha<<") "<<std::endl;
			writemha<<"(starting_regularization ";
			writemha<<Convert_Double_To_CharArray(GetStartRegularization());
			writemha<<") "<<std::endl;
			writemha<<"(ending_regularization ";
			writemha<<Convert_Double_To_CharArray(GetEndRegularization());
			writemha<<") "<<std::endl;
			writemha<<"(optimization_iterations ";
			writemha<<Convert_Double_To_CharArray(GetOptimizationIteration());
			writemha<<") "<<std::endl;
			writemha<<"(checkpointing_interval ";
			writemha<<Convert_Double_To_CharArray(GetCheckpointingInterval());
			writemha<<") "<<std::endl;

			if(GetProcrustesScalingState()){
				writemha<<"(relative_weighting ";
				writemha<<Convert_Double_To_CharArray(GetRelativeWeighting());
				writemha<<") "<<std::endl;
			}
			if(GetProcrustesIntervalState()){
				writemha<<"(adaptivity_strength ";
				writemha<<Convert_Double_To_CharArray(GetAdaptivityStrength());
				writemha<<") "<<std::endl;
			}
			if(GetAdaptivityStrengthState()){
				writemha<<"(procrustes_interval ";
				writemha<<Convert_Double_To_CharArray(GetProcrustesInterval());
				writemha<<") "<<std::endl;
			}
			if(GetRelativeWeightingState()){
				writemha<<"(procrustes_scaling ";
				writemha<<Convert_Double_To_CharArray(GetCheckpointingInterval());
				writemha<<") "<<std::endl;
			}


			std::string output;
			output.append(GetOutputDirectory());
			output.append("/Corresponding_Particles/corresp");
			writemha<<"(output_points_prefix \"";
			writemha<<output;
			writemha<<"\")"<<std::endl;
			writemha.close();

		}
	std::cout<<Correspondence_file<<std::endl;
	std::cout<<" "<<std::endl;
}

void ParticleModuleParameters::CommandlineShapeWorksRun()
{

	args.clear();
	args.push_back("ShapeWorksRun" );
	args.push_back(Correspondence_file.c_str() );
	args.push_back(0);

	std::cout<<"ShapeWorksRun command line: "<<std::endl;
	for( unsigned int n =0; n<args.size()-1;n++)
	{std::cout<<args.at(n);}
	std::cout<<" "<<std::endl;
}

void ParticleModuleParameters::RunShapeWorksRun()
{
	CommandlineShapeWorksRun();

	std::cout<<" "<<std::endl;
	std::cout<<"....Run ShapeWorksRun.... "<<std::endl;

	Run(args,1);

}



//.....................................................................................
//                                PostProcessing
//....................................................................................


void ParticleModuleParameters::CreateVTKFiles()
{
	std::cout<<" "<<std::endl;
	std::cout<<"....Post processing : creation of the meshes and of the mrml scene.... "<<std::endl;

	for(int j=1;j<3;j++)

	{
		if(j==0){FindPostLptsFiles("/Corresponding_Particles/corresp*.mode",PostLptsWptsFile);}
		if(j==2){FindPostLptsFiles("/Corresponding_Particles/corresp*.lpts",PostLptsWptsFile);}
		if(j==1){FindPostLptsFiles("/Corresponding_Particles/corresp*.wpts",PostLptsWptsFile);}

		for( unsigned int i =0; i<PostLptsWptsFile.size();i++)
		{
			int result;
			std::string newname;
			size_t found,found2;
			found=PostLptsWptsFile[i].find_last_of("/\\");
			found2=PostLptsWptsFile[i].find_last_of(".");
			newname.append(PostLptsWptsFile[i].substr(0,found+1));
			newname.append(Data2[i]);
			if(j==0){newname.append(".d0.mode");}
			if(j==2){newname.append(".lpts");}
			if(j==1){newname.append(".wpts");
				WPTSFile.push_back(newname);}
			result= rename( PostLptsWptsFile[i].c_str() , newname.c_str() );
			if(j==2){

				LptsToVTK(newname,0); //create in Corresponding_Meshes
				LptsToVTK(newname,1); // create in Corresponding_Meshes_scale
				ScaleMeta(VTKCorrespFile.back(),VTKCorrespFile.back(),GetEnforcedSpaceX(),GetEnforcedSpaceY(),GetEnforcedSpaceZ());
				}
		}
	}

	CreatePsotProcessFiles();
}

//find all the .lpts and stock there names
void ParticleModuleParameters::FindPostLptsFiles(std::string pathFiles,std::vector<std::string> VectorFile)
{
	PostLptsWptsFile.clear();
	itksys::Glob globPostLptsFile;
	std::string pathFile;
	pathFile=GetOutputDirectory()+pathFiles;
	globPostLptsFile.FindFiles(pathFile);
	PostLptsWptsFile=globPostLptsFile.GetFiles();

}

void ParticleModuleParameters::LptsToVTK(std::string lptsFile, bool scale)
{

	std::ifstream lptsfile(lptsFile.c_str(), std::ios::in); 

	std::string line,word;
	int compt =0;

	std::string namevtk = SetVTKName(lptsFile,scale);

	std::string stringtemplatevtk;
	stringtemplatevtk.append(Data.at(0));


	//double pt[3];
	vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
	int nbpt=0;

	// Create a set of points


	vtkSmartPointer<vtkPoints> pointsSource=vtkSmartPointer<vtkPoints>::New();
	
	while(getline(lptsfile, line)) 
		{	

			for(unsigned int i=0;i<line.length();i++)
			{
				if ((line.at(i)!=' ')&&(i!=line.length()-1))
				{	word=word+line.at(i);	}
	
				else
				{
					double pt[3];
					pt[compt]=string_to_double( word );
					compt++;
					word.clear();
					if(compt==3){
						compt=0;
						 pointsSource->InsertNextPoint( pt);
						}
				}

			}nbpt++;
		}

	lptsfile.close();	
	vtkPolyData *polyData = vtkPolyData::New();
	polyData->SetPoints(pointsSource);
	polyData->Update();

	//cells
	vtkPolyDataReader *meshin = vtkPolyDataReader::New(); 
	meshin->SetFileName(Data.at(0).c_str());
	meshin->Update();

	vtkCellArray *polys;
 	polys = meshin->GetOutput()->GetPolys();
	vtkSmartPointer<vtkCellArray >polysout = vtkSmartPointer<vtkCellArray>::New();

	int prim = 0;
	vtkIdType npts, *pts;

	
	for (polys->InitTraversal(); polys->GetNextCell(npts, pts); prim++)
	{
		polysout->InsertNextCell(npts,pts);
	}
	polyData->SetPolys(polysout);
	polyData->Update();
	SurfaceWriter->SetInput(polyData);
	SurfaceWriter->SetFileName(namevtk.c_str()); 
	SurfaceWriter->Update();
}


std::string ParticleModuleParameters::SetVTKName(std::string lptsName, bool scale)
{

	size_t found;
	std::string vtkName;
	found=lptsName.find_last_of("/");
	vtkName=lptsName.substr(found);
	if(scale==0){ vtkName.insert(0,"/Corresponding_Meshes");}
	if(scale==1){ vtkName.insert(0,"/Corresponding_Meshes_scale");}
	vtkName.insert(0,GetOutputDirectory());
	vtkName.erase(vtkName.end()-5,vtkName.end());
	if(scale==0){vtkName.append("_corr.vtk");
			VTKCorrespFile.push_back(vtkName);}
	if(scale==1){vtkName.append("_corr_scale.vtk");}
	return(vtkName);
}

void ParticleModuleParameters::CreatePsotProcessFiles()
{
	PP_file.append(GetOutputDirectory());

	PP_file.append("/Corresponding_Particles/ShapeWorksView.params");
	std::ofstream PPfile(PP_file.c_str(), std::ios::out);
	if(PPfile){

			PPfile<<"//ShapeWorksView.params"<<std::endl;
			PPfile<<"//Parameters for viewing correspondence point using ShapeWorksView"<<std::endl;
			PPfile<<"(point_files"<<std::endl;
			for(unsigned int i =0;i<WPTSFile.size();i++)
			{
				PPfile<<"\""<<WPTSFile.at(i)<<"\" "<<std::endl;
			}
			PPfile<<")"<<std::endl;
	}
}


//.....................................................................................
//                               MRML
//....................................................................................

void ParticleModuleParameters::CreateMrml()
{
	int nbShapesPerMRML= GetHorizontalGridPara() * GetVerticalGridPara();

	int DataNumber=VTKCorrespFile.size();



	int nbVTKlastMRML,nbMRML;
	if(DataNumber>nbShapesPerMRML){
		nbVTKlastMRML= DataNumber %nbShapesPerMRML;
		if(nbVTKlastMRML!=0){
			nbMRML= (DataNumber-nbVTKlastMRML)/nbShapesPerMRML + 1;}
		else{
			nbMRML= (DataNumber-nbVTKlastMRML)/nbShapesPerMRML;
		}
	}
	else{nbMRML=0;}


	for(int i=0; i<(nbMRML+1);i++)//+1 since the 1st is with all the datas
	{
		std::string mrmlfile;
		std::vector<std::string> transformfile ;
		std::vector<std::string> NameTrans;
		std::vector<std::string> Name;
		std::vector<std::string> Namevtk;
		std::vector<std::string> RelativePathToVTK;
		std::vector<std::string> NbTrans;
		std::vector<std::string> NbFidu;
		int first, last, DataNumberPerMRML;
		double dim0,dim1,dim2;
		int count_line, count_col,nbdisplay;
		count_line=1;count_col=0;nbdisplay=0;	
		dim0=0;dim1=0;dim2=0;

		std::vector<const char*> argsMRML;
		first=(i-1)*nbShapesPerMRML;
		last=(i*nbShapesPerMRML)-1;

		//create the mrml name
		mrmlfile.append(GetOutputDirectory());
		mrmlfile.append("/MRML/ParticleModuleMRMLscene");
		if(i==0)
		{
			mrmlfile.append(".mrml");
			DataNumberPerMRML=DataNumber;
		}
		else{
			mrmlfile.append("_");
			mrmlfile.append(Convert_Double_To_CharArray(first));
			mrmlfile.append("_");
			mrmlfile.append(Convert_Double_To_CharArray(last));
			mrmlfile.append(".mrml");
			DataNumberPerMRML=nbShapesPerMRML;
		}
		std::cout<<" "<<std::endl;
		std::cout<<"mrml scene : "<<mrmlfile<<std::endl;	
		argsMRML.push_back("CreateMRML");   
		argsMRML.push_back(mrmlfile.c_str());
		
		for(int k=0;k<DataNumberPerMRML;k++)
		{
			size_t found,found2;
			std::string relativeName,NameFile,TransName;
			vector<double>Dims;

			found=(VTKCorrespFile.at(k)).find_last_of("/\\");
			relativeName.append((VTKCorrespFile.at(k)).substr(found+1));
			Namevtk.push_back(relativeName);
			relativeName.insert(0,"../Corresponding_Meshes/");
			RelativePathToVTK.push_back(relativeName);
			found2=(Namevtk.back()).find_last_of(".\\");
			TransName.append((Namevtk.back()).substr(0,found2));
			Name.push_back(TransName);
			TransName.append("Trans");
			NameTrans.push_back(TransName);

			//know the orientation
			if(i==0) {SetImageDimensions((char *)(VTKCorrespFile.at(0)).c_str());}
			else{SetImageDimensions((char *)(VTKCorrespFile.at(first)).c_str());}
			Dims=GetImageDimensions();
	
			//calculation of the transforms
			if(GetDirectionToDisplay()=="XYZ")
			{
				dim0 = Dims[0] * count_col;
				dim1 = Dims[1] * count_line;
			}
			if(GetDirectionToDisplay()=="XZY")
			{
				dim0 = Dims[0] * count_col;
				dim2 = Dims[2] * count_line;
			}
			if(GetDirectionToDisplay()=="YXZ")
			{
				dim1 = Dims[1] * count_col;
				dim0 = Dims[0] * count_line;
			}
			if(GetDirectionToDisplay()=="YZX")
			{
				dim1 = Dims[1] * count_col;
				dim2 = Dims[2] * count_line;
			}
			if(GetDirectionToDisplay()=="ZXY")
			{
				dim2 = Dims[2] * count_col;
				dim0 = Dims[0] * count_line;
			}
			if(GetDirectionToDisplay()=="ZYX")
			{
				dim2 = Dims[2] * count_col;
				dim1 = Dims[1] * count_line;
			}
	
			//create the transform file.
			std::string tmp_transformfile ;  
			tmp_transformfile.append("./TransformFiles/transform");
			if(i==0) {tmp_transformfile.append(Convert_Double_To_CharArray(k));}
			else{tmp_transformfile.append(Convert_Double_To_CharArray(DataNumber+first+k));}
			tmp_transformfile.append(".txt");
			transformfile.push_back(tmp_transformfile);
			argsMRML.push_back("-t" );
			argsMRML.push_back("-f" );
			argsMRML.push_back((transformfile.back()).c_str());   

			argsMRML.push_back("-n" );
			argsMRML.push_back((NameTrans.back()).c_str());
			argsMRML.push_back("-l" );
			std::string tmp_NbTrans;  
			if(GetDirectionToDisplay()=="XYZ")
			{
				tmp_NbTrans = "1,0,0,0,1,0,0,0,1,";
				tmp_NbTrans += Convert_Double_To_CharArray(dim0) ;
				tmp_NbTrans += "," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim1) ;
				tmp_NbTrans += ",0" ;
			}
			if(GetDirectionToDisplay()=="XZY")
			{

				tmp_NbTrans = "1,0,0,0,1,0,0,0,1,";
				tmp_NbTrans += Convert_Double_To_CharArray(dim0) ;
				tmp_NbTrans += ",0," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim2*-1) ;
			}
			if(GetDirectionToDisplay()=="YXZ")
			{
				tmp_NbTrans = "1,0,0,0,1,0,0,0,1,";
				tmp_NbTrans += Convert_Double_To_CharArray(dim0) ;
				tmp_NbTrans += "," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim1) ;
				tmp_NbTrans += ",0" ;
			}
			if(GetDirectionToDisplay()=="YZX")
			{
				tmp_NbTrans = "1,0,0,0,1,0,0,0,1,0," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim1) ;
				tmp_NbTrans += "," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim2*-1) ;
			}
			if(GetDirectionToDisplay()=="ZXY")
			{
				tmp_NbTrans = "1,0,0,0,1,0,0,0,1,";
				tmp_NbTrans += Convert_Double_To_CharArray(dim0) ;
				tmp_NbTrans += ",0," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim2*-1) ;
			}
			if(GetDirectionToDisplay()=="ZYX")
			{
				tmp_NbTrans = "1,0,0,0,1,0,0,0,1,0," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim1) ;
				tmp_NbTrans += "," ;
				tmp_NbTrans += Convert_Double_To_CharArray(dim2*-1) ;
			}
			NbTrans.push_back(tmp_NbTrans);
			argsMRML.push_back((NbTrans.back()).c_str());

		
			//add shape
			argsMRML.push_back("-m" ); argsMRML.push_back("-f" ); argsMRML.push_back((RelativePathToVTK.back()).c_str()); argsMRML.push_back("-n");  argsMRML.push_back(Namevtk.back().c_str()); 


			//link shape and transform
			argsMRML.push_back("-p" );
			argsMRML.push_back((NameTrans.back()).c_str());

			//add fiducial
			argsMRML.push_back("-q"); argsMRML.push_back("-id"); 
			argsMRML.push_back((Name.back()).c_str()); 
			argsMRML.push_back("-lbl"); 
			argsMRML.push_back((Name.back()).c_str());
			argsMRML.push_back("-sc");
			argsMRML.push_back("0,0,0");
			argsMRML.push_back("-pos"); 

			std::string tmp_NbFidu;
			if(GetDirectionToDisplay()=="XYZ")
			{
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[3]+dim0+Dims[0]/2));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[6]+dim1));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[7]));

			}
			if(GetDirectionToDisplay()=="XZY")
			{
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[3]+dim0+Dims[0]/2));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[5]));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[8]+dim2));
			}
			if(GetDirectionToDisplay()=="YXZ")
			{
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[4]+dim0));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[5]+dim1+Dims[1]/2));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[7]));
			}
			if(GetDirectionToDisplay()=="YZX")
			{
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[3]));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[5]+dim1+Dims[1]/2));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[8]+dim2));
			}
			if(GetDirectionToDisplay()=="ZXY")
			{
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[4]+dim0));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[5]));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[7]+dim2+Dims[3]/2));
			}
			if(GetDirectionToDisplay()=="ZYX")
			{
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[3]));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[6]+dim1));
				tmp_NbFidu.append(",");
				tmp_NbFidu.append(Convert_Double_To_CharArray(Dims[7]+dim2+Dims[3]/2));
			}
			NbFidu.push_back(tmp_NbFidu);
			argsMRML.push_back((NbFidu.back()).c_str());

			nbdisplay++;

			//how many .vtk per line
			if(i==0){
				if(nbdisplay!=10) {count_col++;}
				else{
				count_col=0;
				count_line++;
				nbdisplay=0;}
			}
			else{
				if(nbdisplay!=GetHorizontalGridPara()) {count_col++;}
				else{
				count_col=0;
				count_line++;
				nbdisplay=0;}
			}


		}

		argsMRML.push_back(0);
		Run(argsMRML,0);
	}
}



// Set the dimension of the first object to display it in a MRML scene
void ParticleModuleParameters::SetImageDimensions(char *filename)
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
	m_Dims.push_back(minCoord[0]);//3
	m_Dims.push_back(maxCoord[0]);
	m_Dims.push_back(minCoord[1]);//5
	m_Dims.push_back(maxCoord[1]);
	m_Dims.push_back(minCoord[2]);//7
	m_Dims.push_back(maxCoord[2]);

}

vector <double> ParticleModuleParameters::GetImageDimensions()
{
	return m_Dims;
}

std::string ParticleModuleParameters::GetDirectionToDisplay()
{
	return m_directionToDisplay;
}









