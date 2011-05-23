#include "ShapeAnalysisModuleComputation.h"
#include <itksys/SystemTools.hxx>

ShapeAnalysisModuleComputation::ShapeAnalysisModuleComputation()
: Parameters()
{}

ShapeAnalysisModuleComputation::~ShapeAnalysisModuleComputation()
{}

// Compute Shape Analysis
void ShapeAnalysisModuleComputation::Computation()
{

	std::cout<<"\n\nComputing ShapeAnalysisModule..."<<std::endl<<std::endl;

	SetBMSShapeAnalysisModuleFile(false);
	SetBMSShapeAnalysisModuleMRMLFile(false);

	m_nbHorizontal=GetHorizontalGridPara();
	m_nbVertical=GetVerticalGridPara();
	nbShapesPerMRML= m_nbHorizontal * m_nbVertical;
	int nummrml=-1;//if -1 you will have one all the shapes 

	SetAllFilesName();
	OverWrite();
	WriteBMSShapeAnalysisModuleFile();	
	ExecuteBatchMake(GetBMSShapeAnalysisModuleFile());std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;

	if(GetTemplateMState()==true)
	{
		ComputationMean();

		WriteBMSShapeAnalysisModuleFile2();
		ExecuteBatchMake(GetBMSShapeAnalysisModuleFile2());
	}
	
	std::cout<<"\n\nExecute Meshmath..."<<std::endl<<std::endl;
	//execute MeshMath external application
	for(int i=0;i<GetDataNumber();i++)
	{
		ExecuteMeshMath(i,"phi",0);
		ExecuteMeshMath(i,"theta",0);
	}
	
	ExecuteMeshMathTemplate();

	// Delete the transform file;
	for(int type=0; type<3;type++)
	{ DeleteTransformsFolders(type);}
			
	int DataNumber=GetDataNumber();
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
	for(int i=0;i<nbMRML+1;i++){ //+1 since the 1st is with all the datas
		if(nummrml!=-1){SetFilesNameMRML(nummrml);}
		WriteBMSMRMLScene(nummrml);		
		nummrml++;
	}	

 /*	ExecuteBatchMake(GetBMSShapeAnalysisModuleMRMLFile()); 
	SetBMSShapeAnalysisModuleFile(true);*/
	/*std::cout<<"Modify output csv"<<std::endl;
	ModifyCSV();

	std::cout<<"Computing ShapeAnalysisModule: Done!"<<std::endl<<std::endl;*/

	//Particles 
	if(GetParticlesState())
	{
		RunParticlesModule();
		std::cout<<"Modify output csv"<<std::endl;
		ModifyCSV(1);
		for(int i=0;i<GetDataNumber();i++)
		{
			ExecuteMeshMath(i,"phi",1);
			ExecuteMeshMath(i,"theta",1);
		}
		std::cout<<"MRML"<<std::endl;
		CreateMrmlParticle();
		

	}
	else {ModifyCSV(0);}
  
  	return;
}

//Execute MeshMath to write a KWM scalar field (1D) into a PolyData Field Data Scalar to visualize in Slicer3
void ShapeAnalysisModuleComputation::ExecuteMeshMath(int numData, char * scalar, bool particule)
{
	int end;
	if(particule ==0) {end=3;}//execute MeshMath for each volume file: Original, Ellalign, Procalign
	else {end=1;}


	for(int j=0;j<end;j++)
	{
	
		std::vector<const char*> args;  
		char* data = NULL;
		int length;
		double timeout = 0.05;
		int result;
		char *fileType=NULL;
		
		args.push_back("MeshMath");
		
		if(particule ==0)
		{
			if (GetTemplateMState()==true) 
			{
				if(j==0) fileType=GetAllSurfmeanSPHARMFiles(numData);
				else if (j==1) fileType=GetAllSurfmeanSPHARMellalignFiles(numData);
				else if (j==2) fileType=GetAllSurfmeanSPHARMprocalignFiles(numData);
			}
			else
			{
				if(j==0) fileType=GetAllSurfSPHARMFiles(numData);
				else if (j==1) fileType=GetAllSurfSPHARMellalignFiles(numData);
				else if (j==2) fileType=GetAllSurfSPHARMprocalignFiles(numData);
			}
		}
		else
		{	std::string tmp;
			if (GetTemplateMState()==true) 
			{
			}
			else
			{
				fileType=GetPostCorrespondenceFiles(numData);
			}
		}
	
		args.push_back(fileType);
		args.push_back(fileType);
		args.push_back("-KWMtoPolyData");
		
		if(scalar=="phi")
		{  
			args.push_back(GetAllPhiFiles(0));  
			args.push_back("Color_Map_Phi");
		}
		if(scalar=="theta")
		{ 
			args.push_back(GetAllThetaFiles(0));
			args.push_back("Color_Map_Theta");
		}
		
		args.push_back(0);
		
		
		// Run the application
		itksysProcess* gp = itksysProcess_New();
		itksysProcess_SetCommand(gp, &*args.begin());
		itksysProcess_SetOption(gp,itksysProcess_Option_HideWindow,1);
		itksysProcess_Execute(gp);
	
	
		while(int Value = itksysProcess_WaitForData(gp,&data,&length,&timeout)) // wait for 1s
		{
		if ( ((Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR)) && data[0]=='D' )
		{
			strstream st;
			for(int i=0;i<length;i++) 	
			{
				st<<data[i];
			}
			string dim=st.str();
			istringstream s(dim);
			string value;
					
			while(getline(s,value,' '))
			{	
				m_Dims.push_back((atoi(value.c_str()))/2);		
			}
			
			m_Dims.erase(m_Dims.begin());
		}
			timeout = 0.05;   	
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
}


void ShapeAnalysisModuleComputation::ExecuteMeshMathTemplate()
{
	std::cout<<"meshmath template"<<std::endl;

	itksys::Glob globTemplateList;
	std::vector<std::string>  list_template;
	std::string pathFile =GetOutputDirectory();
	std::string path;
	path="/Template/*.vtk";
	pathFile=pathFile+path;
	globTemplateList.FindFiles(pathFile);
	list_template=globTemplateList.GetFiles();

	for(int color_map =0;  color_map<2; color_map++)
	{
		for(unsigned int  i=0;i<list_template.size();i++)
		{
			std::vector<const char*> args;  
			char* data = NULL;
			int length;
			double timeout = 0.05;
			int result;
			
			args.push_back("MeshMath");
			
				
			args.push_back(list_template[i].c_str());
			args.push_back(list_template[i].c_str());
			args.push_back("-KWMtoPolyData");
			
			if( color_map==0)
			{  
				args.push_back(GetAllPhiFiles(0));  
				args.push_back("Color_Map_Phi");
			}
			if( color_map==1)
			{ 
				args.push_back(GetAllThetaFiles(0));
				args.push_back("Color_Map_Theta");
			}
			
			args.push_back(0);
			
			
			// Run the application
			itksysProcess* gp = itksysProcess_New();
			itksysProcess_SetCommand(gp, &*args.begin());
			itksysProcess_SetOption(gp,itksysProcess_Option_HideWindow,1);
			itksysProcess_Execute(gp);
		
		
			while(int Value = itksysProcess_WaitForData(gp,&data,&length,&timeout)) // wait for 1s
			{
			if ( ((Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR)) && data[0]=='D' )
			{
				strstream st;
				for(int i=0;i<length;i++) 	
				{
					st<<data[i];
				}
				string dim=st.str();
				istringstream s(dim);
				string value;
						
				while(getline(s,value,' '))
				{	
					m_Dims.push_back((atoi(value.c_str()))/2);		
				}
				
				m_Dims.erase(m_Dims.begin());
			}
				timeout = 0.05;   	
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
		
	}
}

// Create BMS File to compute the SPHARM pipeline
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleFile(bool changeDirectory)
{
	std::strcpy(m_BMSShapeAnalysisModuleFile, GetOutputDirectory());
	std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;

	if(changeDirectory==false)
	{	std::strcat(m_BMSShapeAnalysisModuleFile, "/");
std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;}

	else {std::strcat(m_BMSShapeAnalysisModuleFile, "/BatchMake_Scripts/");
std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;}

std::cout<<GetRandomizeInputs()<<std::endl;

	if (GetRandomizeInputs())
	{
		int pID=getpid();
		sprintf(m_BMSShapeAnalysisModuleFile,"ShapeAnalysisModule_id%d.bms",pID);
		std::cout << " " << std::endl;std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;
	}	
	else{
			
		std::strcat(m_BMSShapeAnalysisModuleFile, "ShapeAnalysisModule.bms");std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;}

	std::cout<<GetBMSShapeAnalysisModuleFile()<<std::endl;
	return;  
}

 // Get BMS File
char *ShapeAnalysisModuleComputation::GetBMSShapeAnalysisModuleFile()
{
  	return m_BMSShapeAnalysisModuleFile;
}


// Create BMS for mean computation File
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleFile2(bool changeDirectory)
{
	std::strcpy(m_BMSShapeAnalysisModuleFile, GetOutputDirectory());
	
	if(changeDirectory==false)
		std::strcat(m_BMSShapeAnalysisModuleFile, "/");

	else std::strcat(m_BMSShapeAnalysisModuleFile, "/BatchMake_Scripts/");

	std::strcat(m_BMSShapeAnalysisModuleFile, "ShapeAnalysisModule_MeanAsTemplate.bms");

	return;  
}


char * Convert_Double_To_CharArray(double doubleVariable) {
	char *CharArrayOut;
	CharArrayOut = new char [512];
	std::string stringVariable;
	std::stringstream strStream;
	strStream << doubleVariable;
	stringVariable = strStream.str();
	strcpy(CharArrayOut,stringVariable.c_str());

	return CharArrayOut;
}

// Get BMS File for mean computation
char *ShapeAnalysisModuleComputation::GetBMSShapeAnalysisModuleFile2()
{
  	return m_BMSShapeAnalysisModuleFile;
}

// Create a BMS file to create a MRML scene
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleMRMLFile(bool changeDirectory)
{
	std::strcpy(m_BMSShapeAnalysisModuleMRLMFile, GetOutputDirectory());
	std::strcat(m_BMSShapeAnalysisModuleMRLMFile, "/BatchMake_Scripts/");
	std::strcat(m_BMSShapeAnalysisModuleMRLMFile, "ShapeAnalysisModuleMRML.bms");
}

char * ShapeAnalysisModuleComputation::GetBMSShapeAnalysisModuleMRMLFile()
{
	return m_BMSShapeAnalysisModuleMRLMFile;
}

//create the output file
void ShapeAnalysisModuleComputation::SetOuputFile()
{
	std::strcpy(m_OutputFile, GetOutputDirectory());
	std::strcat(m_OutputFile, "/OutputGroupFile/");
	std::strcat(m_OutputFile, "ShapeAnalysisModule_OutputFileVersion1.csv");
}

char * ShapeAnalysisModuleComputation::GetOutputFile()
{
	return m_OutputFile;
}

 // Execute BatchMake script
void ShapeAnalysisModuleComputation::ExecuteBatchMake(char *_Input)
{
  	std::cout<<"\tExecuting BatchMake..."<<std::endl;

	char * envpath = getenv("BatchmakeShapeAnalysisModule_Dir");
	std::string applicationPath;
	bm::ScriptParser m_Parser;

	// the module path exist
	if(chdir(GetModulePath())==0)
	{	
		//if bmm files are in the same directory as the Shape Analysis Module
		m_Parser.LoadWrappedApplication(GetModulePath());
	}

	else
	{
		// if the environment variable is set
		if (envpath) 
		{	
			applicationPath = std::string(envpath);
			m_Parser.LoadWrappedApplication(applicationPath.c_str());
		}

		else
		{
			std::cerr<<"The environment variable 'BatchmakeShapeAnalysisModule_Dir' needs to be set"<<std::endl;
			std::cerr<<"bash usage : export BatchmakeShapeAnalysisModule_Dir=<Batchmake ShapeAnalysisModule Directory>"<<std::endl;
			std::cerr<<"tcsh usage : setenv BatchmakeShapeAnalysisModule_Dir <Batchmake ShapeAnalysisModule Directory>"<<std::endl;
			exit(0);
		}
	}
		
	
	bool ComputationSuccess = m_Parser.Execute(_Input);
	if (!ComputationSuccess)
		cerr<<"\tExecuting BatchMake: Error!"<<endl;

 	 std::cout<<"\tExecuting BatchMake: Done!"<<std::endl<<endl;
  return;
}

// Write the BMS script to create a MRML file
void ShapeAnalysisModuleComputation::WriteBMSMRMLScene(int whichmrml)
{


	for(int count=0;count<3;count++) // None,ecalign,procallign
	{

		for(int nbcolormap=0;nbcolormap<3;nbcolormap++)
		{
			std::string mrmlfilePhi;
			std::string mrmlfileTheta;

	
			double dim0,dim1,dim2;
			int count_line, count_col,nbdisplay;
			count_line=1;count_col=0;nbdisplay=1;	
			dim0=0;dim1=0;dim2=0;
		
			int DataNumber;
			int Nb_Data=GetDataNumber();
			int first, last;	
			if( whichmrml==-1){ DataNumber=GetDataNumber();}
			else{
				first=whichmrml*nbShapesPerMRML;
				last=((whichmrml+1)*nbShapesPerMRML)-1;
			
				//to know if it's the last mrml
				if((Nb_Data>=first) && (Nb_Data<=last))
				{
					DataNumber= Nb_Data%nbShapesPerMRML;
				}
				else
				{
					DataNumber=nbShapesPerMRML;
				}
			}	
		
			std::vector<std::string> Name;
			std::vector<std::string> NameVTK;
			std::vector<std::string> NameTrans;
			std::vector<std::string> NbFidu; 
			std::vector<std::string> NameTemplate;
			std::vector<std::string> randomcolor;
			std::vector<std::string> transformfile ;
			std::vector<std::string> NbTrans;
			std::vector<std::string> NameColormap;
	
			//init all names
			for(int i=0;i<DataNumber;i++)
			{	
				std::string help_Name;
				if( whichmrml==-1) {help_Name.append(GetAllFilesName(i));}
				else{help_Name.append(GetListFiles(i));}
				if(GetTemplateMState()){help_Name.append("_pp_surf_tMeanSPHARM");}
				else {help_Name.append("_pp_surfSPHARM");}
				if(count==1){help_Name.append("_ellalign");}
				if(count==2){help_Name.append("_procalign");}
				Name.push_back(help_Name);
			
				std::string help_NameTrans;
				if( whichmrml==-1){help_NameTrans.append(GetAllFilesName(i));}
				else{help_NameTrans.append(GetListFiles(i));}
				if(GetTemplateMState()){help_NameTrans.append("_tMean");}
				if(count==1){help_NameTrans.append("_ellalign");}
				if(count==2){help_NameTrans.append("_procalign");}
				help_NameTrans.append("Trans");
				NameTrans.push_back(help_NameTrans);
				
				std::string help_NameVTK;
				if(count==1){help_NameVTK.append("../");}
				if(count==2){help_NameVTK.append("../");}
				help_NameVTK.append("../Mesh/SPHARM/");
				if( whichmrml==-1){help_NameVTK.append(GetAllFilesName(i));}
				else{help_NameVTK.append(GetListFiles(i));}
				if(GetTemplateMState()){help_NameVTK.append("_pp_surf_tMeanSPHARM");}
				else {help_NameVTK.append("_pp_surfSPHARM");}
				if(count==1){help_NameVTK.append("_ellalign");}
				if(count==2){help_NameVTK.append("_procalign");}
				help_NameVTK.append(".vtk");
				NameVTK.push_back(help_NameVTK);
			}
		
			//create the mrml name
			std::vector<const char*> args;
			std::string mrmlfile;
			if( whichmrml!=-1){
				mrmlfile.append(GetOutputDirectory());
				mrmlfile.append("/MRML");
				if(count==1){mrmlfile.append("/Ellalign");}
				if(count==2){mrmlfile.append("/Procalign");}
				if(nbcolormap==0){mrmlfile.append("/ShapeAnalysisModuleMRMLscene_Phi");}
				if(nbcolormap==1){mrmlfile.append("/ShapeAnalysisModuleMRMLscene_Theta");}
				if(nbcolormap==2){
					mrmlfile.append("/ShapeAnalysisModuleMRMLscene");
					mrmlfilePhi.append(mrmlfile);
					mrmlfilePhi.append("_Phi");
					mrmlfileTheta.append(mrmlfile);
					mrmlfileTheta.append("_Theta");	}
				if(GetTemplateMState()){
					mrmlfile.append("_tMean_");
					if(nbcolormap==2){
						mrmlfilePhi.append("_tMean_");
						mrmlfileTheta.append("_tMean_");
							}
				}
				mrmlfile.append(Convert_Double_To_CharArray(first));
				if(nbcolormap==2){ 
					mrmlfilePhi.append(Convert_Double_To_CharArray(first));
					mrmlfileTheta.append(Convert_Double_To_CharArray(first));}
				mrmlfile.append("_");
				if(nbcolormap==2){ 
					mrmlfilePhi.append("_");
					mrmlfileTheta.append("_");}
				if((Nb_Data>=first) && (Nb_Data<=last))
				{
					mrmlfile.append(Convert_Double_To_CharArray(first+Nb_Data%nbShapesPerMRML-1));
					if(nbcolormap==2){ 
						mrmlfilePhi.append(Convert_Double_To_CharArray(first+Nb_Data%nbShapesPerMRML-1));
						mrmlfileTheta.append(Convert_Double_To_CharArray(first+Nb_Data%nbShapesPerMRML-1));}
				}
				else
				{
					mrmlfile.append(Convert_Double_To_CharArray(last));
					if(nbcolormap==2){ 
						mrmlfilePhi.append(Convert_Double_To_CharArray(last));
						mrmlfileTheta.append(Convert_Double_To_CharArray(last));}
				}
				if(count==1){
						mrmlfile.append("_ellalign");
						if(nbcolormap==2){ 
							mrmlfilePhi.append("_ellalign");
							mrmlfileTheta.append("_ellalign");}}
				if(count==2){
					mrmlfile.append("_procalign");
						if(nbcolormap==2){ 
							mrmlfilePhi.append("_procalign");
							mrmlfileTheta.append("_procalign");}}
				mrmlfile.append(".mrml");
				if(nbcolormap==2){ 
					mrmlfilePhi.append(".mrml");
					mrmlfileTheta.append(".mrml");}
			}
			else
			{
				mrmlfile.append(GetOutputDirectory());
				mrmlfile.append("/MRML");
				if(count==1){mrmlfile.append("/Ellalign");}
				if(count==2){mrmlfile.append("/Procalign");}
				if(nbcolormap==0){mrmlfile.append("/ShapeAnalysisModuleMRMLscene_allVTK_Phi");}
				if(nbcolormap==1){mrmlfile.append("/ShapeAnalysisModuleMRMLscene_allVTK_Theta");}
				if(nbcolormap==2){mrmlfile.append("/ShapeAnalysisModuleMRMLscene_allVTK");}
				if(GetTemplateMState()){ mrmlfile.append("_tMean");}
				mrmlfile.append(".mrml");
			}
			
			std::cout<<"mrml "<<mrmlfile<<std::endl;	
			args.push_back("CreateMRML");   
			args.push_back(mrmlfile.c_str());

			//know the orientation
			std::string firstFile1;
			firstFile1.append(GetOutputDirectory());
			firstFile1.append("/Mesh/SPHARM/");
			if( whichmrml==-1){firstFile1.append(GetAllFilesName(0));}
			else{firstFile1.append(GetListFiles(0));}
			if(GetTemplateMState()){firstFile1.append("_pp_surf_tMeanSPHARM");}
			else {firstFile1.append("_pp_surfSPHARM");}
			if(count==1){firstFile1.append("_ellalign");}
			if(count==2){firstFile1.append("_procalign");}
			firstFile1.append(".vtk");
			vector<double>Dims;
			SetImageDimensions((char *)(firstFile1).c_str());
			Dims=GetImageDimensions();


			int NbDataFirstMRML;
			NbDataFirstMRML=GetDataNumber();

			for(int i=0;i<DataNumber;i++)
			{


				if((count ==0) ||(count ==1)|| ((count ==2) && (nbcolormap==2))){  //add the template   NB no template for procalign

					if(i==0 &&  whichmrml>=0){
						std::string tmp_NameTemplate;      
						tmp_NameTemplate.append(GetTemplate(count));
						NameTemplate.push_back(tmp_NameTemplate);

						//add shape
						args.push_back("-m" ); args.push_back("-f" ); args.push_back((NameTemplate.back()).c_str()); args.push_back("-n"); args.push_back("template"); 

						if(nbcolormap!=2){
							//add color map
							 
							std::string tmp_NameColormap;
							if(count==1){tmp_NameColormap.append("../");}
							if(count==2){tmp_NameColormap.append("../");}
							if(nbcolormap==0){tmp_NameColormap.append("../Mesh/SPHARM/customLUT_Color_Map_Phi.txt");}
							if(nbcolormap==1){tmp_NameColormap.append("../Mesh/SPHARM/customLUT_Color_Map_Theta.txt");}

							NameColormap.push_back(tmp_NameColormap);

							args.push_back("-as" );
							if(nbcolormap==0){args.push_back("Color_Map_Phi" );}
							if(nbcolormap==1){args.push_back("Color_Map_Theta" );}
							args.push_back("-cc" );
							args.push_back((NameColormap.back()).c_str() );}

				
						else{
							//random color
							std::string tmp_randomcolor;   
							double color;
							color = rand() % 255 + 1;
							tmp_randomcolor.append(Convert_Double_To_CharArray(color/255));
							tmp_randomcolor.append(",");
							color = rand() % 255 + 1;
							tmp_randomcolor.append(Convert_Double_To_CharArray(color/255));
							tmp_randomcolor.append(",");
							color = rand() % 255 + 1;
							tmp_randomcolor.append(Convert_Double_To_CharArray(color/255));
							randomcolor.push_back(tmp_randomcolor);
							args.push_back("-dc" );
							args.push_back((randomcolor.back()).c_str()); }

					}
				}
			

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
	
				if(nbcolormap!=2){
					//create the transform file.
					std::string tmp_transformfile ;  
					tmp_transformfile.append("./TransformFiles/transform");
					if(whichmrml==-1) {tmp_transformfile.append(Convert_Double_To_CharArray(i));}
					if(whichmrml!=-1){tmp_transformfile.append(Convert_Double_To_CharArray(NbDataFirstMRML+i));}
					tmp_transformfile.append(".txt");
					transformfile.push_back(tmp_transformfile);
	
					args.push_back("-t" );
					args.push_back("-f" );
					args.push_back((transformfile.back()).c_str());   
					args.push_back("-n" );
					args.push_back((NameTrans.at(i)).c_str());
					if ((whichmrml==-1 && nbcolormap==0) ||  (whichmrml==0 &&nbcolormap ==0 )){
						args.push_back("-l" );
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
						args.push_back((NbTrans.back()).c_str());
					}	
		
					//add shape
					args.push_back("-m" ); args.push_back("-f" ); args.push_back((NameVTK.at(i)).c_str()); args.push_back("-n");  args.push_back((Name.at(i)).c_str()); 
	
					//add color map
					if(nbcolormap!=2){
						std::string tmp_NameColormap;
						if(count==1){tmp_NameColormap.append("../");}
						if(count==2){tmp_NameColormap.append("../");}
						if(nbcolormap==0){tmp_NameColormap.append("../Mesh/SPHARM/customLUT_Color_Map_Phi.txt");}
						if(nbcolormap==1){tmp_NameColormap.append("../Mesh/SPHARM/customLUT_Color_Map_Theta.txt");}
						NameColormap.push_back(tmp_NameColormap);
						args.push_back("-as" );
						if(nbcolormap==0){args.push_back("Color_Map_Phi" );}
						if(nbcolormap==1){args.push_back("Color_Map_Theta" );}
						args.push_back("-cc" );
						args.push_back((NameColormap.back()).c_str() );
					}
					else{
						//random color
						std::string tmp_randomcolor;   
						double color;
						color = rand() % 255 + 1;
						tmp_randomcolor.append(Convert_Double_To_CharArray(color/255));
						tmp_randomcolor.append(",");
						color = rand() % 255 + 1;
						tmp_randomcolor.append(Convert_Double_To_CharArray(color/255));
						tmp_randomcolor.append(",");
						color = rand() % 255 + 1;
						tmp_randomcolor.append(Convert_Double_To_CharArray(color/255));
						randomcolor.push_back(tmp_randomcolor);
						args.push_back("-dc" );args.push_back((randomcolor.back()).c_str()); }
			
					//link shape and transform
					args.push_back("-p" );
					args.push_back((NameTrans.at(i)).c_str());
	
					//how many .vtk per line
					if(whichmrml==-1){
						if(nbdisplay!=10) {count_col++;}
						else{
						count_col=0;
						count_line++;
						nbdisplay=0;}
					}
					else{
						if(nbdisplay!=m_nbHorizontal) {count_col++;}
						else{
						count_col=0;
						count_line++;
						nbdisplay=0;}
					}
			
					//add fiducial
					args.push_back("-q"); args.push_back("-id"); 
					args.push_back((Name.at(i)).c_str()); 
					args.push_back("-lbl"); 
					args.push_back((Name.at(i)).c_str());
					args.push_back("-sc");
					args.push_back("0,0,0");
					args.push_back("-pos"); 

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
					args.push_back((NbFidu.back()).c_str());
	
					nbdisplay++;
				}//end if !=no colormap
		}//enf for each data

			
		args.push_back(0);


		//write files with the commande line using to create the mrml
		
		char fileCommanLine[512];
		std::strcpy (fileCommanLine,"/biomed-resimg/conte_projects/CONTE_NEO/Data/Vent_Shape/SPHARM/SPHARM_10_15_Slicer/LeftVentLong2/MRML/commandline");
		std::strcat (fileCommanLine,Convert_Double_To_CharArray(whichmrml));
		std::strcat (fileCommanLine,"_");
		std::strcat (fileCommanLine,Convert_Double_To_CharArray(count));
		std::strcat (fileCommanLine,"_");
		std::strcat (fileCommanLine,Convert_Double_To_CharArray(nbcolormap));
		std::strcat (fileCommanLine,".txt");
		std::cout<<fileCommanLine<<std::endl;
		std::ofstream file(fileCommanLine);
	
		if(file)
		{
			for(int k=0; k<args.size();k++)
				{file<<args.at(k)<<" ";}
			file.close();
		}


		//itk sys parameters
			int length;
			double timeout = 0.05;
			int result;
			char* dataitk = NULL;
		std::cout<<""<<std::endl;
		
		// Run the application
		itksysProcess* gp = itksysProcess_New();
		itksysProcess_SetCommand(gp, &*args.begin());
		itksysProcess_SetOption(gp,itksysProcess_Option_HideWindow,1);
		itksysProcess_Execute(gp);
		
		while(int Value = itksysProcess_WaitForData(gp,&dataitk,&length,&timeout)) // wait for 1s
		{
		if ( ((Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR)) && dataitk[0]=='D' )
		{
			strstream st;
			for(int i=0;i<length;i++) 	
			{
				st<<dataitk[i];
			}
			string dim=st.str();
		}
			timeout = 0.05;   	
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
			
		if (whichmrml!=-1){

		if(nbcolormap==2){ //merge the .mrml  

			ModifyMRML(mrmlfile ,mrmlfilePhi , mrmlfileTheta );
		}
		}
		}//end nbcolormap
	
	}//end for count
}


 
 // Write Batchmake script to compute Shape Analysis
void ShapeAnalysisModuleComputation::WriteBMSShapeAnalysisModuleFile()
{

	std::ofstream BMSShapeAnalysisModuleFile(GetBMSShapeAnalysisModuleFile());
	vector<string> OutputFileHeaders=GetOutputFileHeaders();
	int DataNumber;
	
	BMSShapeAnalysisModuleFile<<"#---------------------------- Shape Analysis ----------------------------"<<std::endl;
	BMSShapeAnalysisModuleFile<<"#------------------------------------------------------------------------"<<std::endl;
	BMSShapeAnalysisModuleFile<<"#------------------------------------------------------------------------"<<std::endl;
	BMSShapeAnalysisModuleFile<<"#------------------------------------------------------------------------"<<std::endl;


        for( int i = 0; i < GetDataNumber(); i++)
		m_permutations.push_back(0);
        //get here the permutation vector
        GetRandomNum(1, GetDataNumber());
	//for( int i = 0; i < GetDataNumber(); i++)
	//	std::cout << m_permutations[i] << std::endl;
	
	if (GetRandomizeInputs())
	{
		for(unsigned int i=0;i<OutputFileHeaders.size();i++)
		{
			BMSShapeAnalysisModuleFile<<" Set(Header"<<i<<" '"<<GetNthDataListValue(m_permutations[0],i)<<"')"<<std::endl;
		
		}
	
		BMSShapeAnalysisModuleFile<<"set (OrigCasesList '"<<GetNthDataListValue(m_permutations[0],GetColumnVolumeFile())<<"')"<<std::endl;
	}
	else
	{
	//Set the headers for the output file
		for(unsigned int i=0;i<OutputFileHeaders.size();i++)
		{
			BMSShapeAnalysisModuleFile<<" Set(Header"<<i<<" '"<<GetNthDataListValue(1,i)<<"')"<<std::endl;
		
		}
	
		BMSShapeAnalysisModuleFile<<"set (OrigCasesList '"<<GetNthDataListValue(1,GetColumnVolumeFile())<<"')"<<std::endl;
	}
	
	std::cout<<"Number of Datas: "<<GetDataNumber()<<std::endl;
	
	if (GetRandomizeInputs())
	{
		for (DataNumber = 1; DataNumber <GetDataNumber(); DataNumber++)
		{
		
			for(unsigned int i=0;i<OutputFileHeaders.size();i++)
			{
				BMSShapeAnalysisModuleFile<<"set (Header"<<i<<" ${Header"<<i<<"} '" <<GetNthDataListValue(m_permutations[DataNumber],i)<<"')"<<std::endl;
			
			}
			BMSShapeAnalysisModuleFile<<"set (OrigCasesList ${OrigCasesList} '"<<GetNthDataListValue(m_permutations[DataNumber],GetColumnVolumeFile())<<"')"<<std::endl;
		}
		BMSShapeAnalysisModuleFile<<""<<std::endl;
	}
	else
	{
		for (DataNumber = 2; DataNumber <=GetDataNumber(); DataNumber++)
		{
		
			for(unsigned int i=0;i<OutputFileHeaders.size();i++)
			{
				BMSShapeAnalysisModuleFile<<"set (Header"<<i<<" ${Header"<<i<<"} '" <<GetNthDataListValue(DataNumber,i)<<"')"<<std::endl;
			
			}
			BMSShapeAnalysisModuleFile<<"set (OrigCasesList ${OrigCasesList} '"<<GetNthDataListValue(DataNumber,GetColumnVolumeFile())<<"')"<<std::endl;
		}
		BMSShapeAnalysisModuleFile<<""<<std::endl;
	}
	
	BMSShapeAnalysisModuleFile<<"#Create directories"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/BatchMake_Scripts/)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Mesh/)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML/Procalign)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML/Ellalign)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML/Ellalign/TransformFiles)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML/Procalign/TransformFiles)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML/TransformFiles)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Mesh/PostProcess/)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Mesh/SPHARM/)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Template/)"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/OutputGroupFile/)"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/EulerFiles/)"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"set(BMSdir '"<<GetOutputDirectory()<<"/BatchMake_Scripts')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(datadir '"<<GetOutputDirectory()<<"/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(PPdir '"<<GetOutputDirectory()<<"/Mesh/PostProcess/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(SPHARMdir '"<<GetOutputDirectory()<<"/Mesh/SPHARM/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(tdir '"<<GetOutputDirectory()<<"/Template/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(Outputdir '"<<GetOutputDirectory()<<"/OutputGroupFile/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(Eulerdir '"<<GetOutputDirectory()<<"/EulerFiles/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(transformdir '"<<GetOutputDirectory()<<"/MRML/TransformFiles')"<<std::endl;

/*
BMSShapeAnalysisModuleFile<<" ListFileInDir(Transformfiles ${transformdir} *.txt)"<<std::endl;
BMSShapeAnalysisModuleFile<<"FileExists(fileCreated 'transform0.txt')"<<std::endl;
BMSShapeAnalysisModuleFile<<"If( ${fileCreated} == 1 )"<<std::endl;
BMSShapeAnalysisModuleFile<<"ForEach(files ${Transformfiles})"<<std::endl;
BMSShapeAnalysisModuleFile<<"DeleteFile(${files})"<<std::endl;
BMSShapeAnalysisModuleFile<<"EndForEach(files)"<<std::endl;
BMSShapeAnalysisModuleFile<<"EndIf(${fileCreated})"<<std::endl;*/

	BMSShapeAnalysisModuleFile<<"echo()"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"#Create OutputFile"<<std::endl;
	BMSShapeAnalysisModuleFile<<"Set(OutputFile "<<GetOutputDirectory()<<"/OutputGroupFile/ShapeAnalysisModule_OutputFileVersion1.csv)"<<std::endl;
	SetOuputFile();
	BMSShapeAnalysisModuleFile<<"Set(OutputFile "<<GetOutputFile()<<")"<<std::endl;

	BMSShapeAnalysisModuleFile<<" WriteFile(${OutputFile} '"<<OutputFileHeaders[0]<<",')"<<std::endl;
	for(unsigned int i=1;i<OutputFileHeaders.size();i++)
	{
		BMSShapeAnalysisModuleFile<<" appendFile(${OutputFile} '"<<OutputFileHeaders[i]<<",')"<<std::endl;
	}
	
	//Write headers of the output file
	BMSShapeAnalysisModuleFile<<" appendFile(${OutputFile} ' Post Processed Segmentation, Parameterization of Original Surface, SPHARM Surface in Original Space, SPHARM Coefficient in Original Space, SPHARM Surface in Ellipsoid Aligned Space, SPHARM Coefficient in Ellipsoid Aligned Space, SPHARM Surface in Procaligned Space\\n')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"echo()"<<std::endl;
	
	
	BMSShapeAnalysisModuleFile<<"  Set(count 0)"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<""<<std::endl;
	BMSShapeAnalysisModuleFile<<"ForEach(case ${OrigCasesList})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  #Extract basename"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  GetFilename(basename ${case} NAME_WITHOUT_EXTENSION)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo(Case: ${case})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testSeg ${PPdir} *${basename}*_pp"<<GetVolumeFileExtension()<<")"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testGen1 ${SPHARMdir} *${basename}*para.vtk*) "<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testGen2 ${SPHARMdir} *${basename}*surf.vtk*) "<<std::endl;
	BMSShapeAnalysisModuleFile<<"  set(testGen ${testGen1} ${testGen2}) "<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testPara ${SPHARMdir} *${basename}*SPHARM*)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testPara2 ${SPHARMdir} *${basename}*MeanSPHARM*)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testtemp ${tdir}  *SPHARM*)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  listFileInDir(testtemp2 ${tdir} *${basename}*SPHARM*)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  "<<std::endl;
	BMSShapeAnalysisModuleFile<<"  set(ppcase ${PPdir}${basename}_pp.gipl.gz)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  "<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  #Post Processing"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo('Doing Post Processing')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  if(${testSeg} == '')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetApp(Seg @SegPostProcessCLP)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Seg.fileName ${case})"<<std::endl;

	if (GetGaussianFilteringState()==true)
	{
		BMSShapeAnalysisModuleFile<<"    SetAppOption(Seg.gaussianOn 1)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(Seg.variance_vect.variance_vect "<<GetVarianceBoxX()<<","<<GetVarianceBoxY()<<","<<GetVarianceBoxZ()<<")"<<std::endl;
	}
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Seg.outfileName ${ppcase})"<<std::endl;
	if (GetLabelState()==true)
	{
		BMSShapeAnalysisModuleFile<<"    SetAppOption(Seg.label.label "<<GetLabel()<<")"<<std::endl;
	}
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Seg.spacing_vect.spacing_vect "<<GetEnforcedSpaceX()<<","<<GetEnforcedSpaceY()<<","<<GetEnforcedSpaceZ()<<")"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    Run(output ${Seg} error)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    echo('SegPostProcess Error:' ${error}')"<<std::endl;
	//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  endif(${testSeg})"<<std::endl;
	BMSShapeAnalysisModuleFile<<""<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  #GenParaMesh"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo('Doing GenParaMesh')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  Set(value 0)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  foreach(data ${testGen})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    Inc(${value} 1)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    Int(${value})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  EndForEach(data)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  if (${value} < 2)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetApp(Gen @GenParaMeshCLP)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.infile ${ppcase})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.outParaName ${SPHARMdir}${basename}_pp_para.vtk)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.outSurfName ${SPHARMdir}${basename}_pp_surf.vtk)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.numIterations.numIterations "<<GetNumberOfIterations()<<")"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.label.label 1)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.EulerFile 1)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.outEulerName 1)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    SetAppOption(Gen.outEulerName.outEulerName ${Eulerdir}${basename}_euler.txt)"<<std::endl;



	BMSShapeAnalysisModuleFile<<"    Run(output ${Gen} error)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"    echo('GenParaMesh Error:' ${error})"<<std::endl;
	//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  endif(${value})"<<std::endl;
	BMSShapeAnalysisModuleFile<<""<<std::endl;

	BMSShapeAnalysisModuleFile<<"  #ParaToSPHARMMesh"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo('Doing ParaToSPHARMMesh')"<<std::endl;

	
		BMSShapeAnalysisModuleFile<<"  Set(value 0)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  foreach(data ${testtemp})"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    Inc(${value} 1)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    Int(${value})"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  EndForEach(data)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  if (${value} < 4)"<<std::endl;

//TODO create templte

		BMSShapeAnalysisModuleFile<<"  #Create Template"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    echo('Creating Template')"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetApp(ParaT @ParaToSPHARMMeshCLP)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.outbase ${tdir}${basename}_pp_surf)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.subdivLevel 1)" <<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.subdivLevel.subdivLevel "<<GetSubdivLevel()<<")"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.spharmDegree 1)" <<std::endl;
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.spharmDegree.spharmDegree "<<GetSPHARMDegree()<<")"<<std::endl;
		BMSShapeAnalysisModuleFile<<"      SetAppOption(ParaT.finalFlipIndex 1)" <<std::endl;
		if (GetFinalFlipN()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex "<<0<<")"<<std::endl;}
		if (GetFinalFlipX()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<4<<")"<<std::endl;}
		if (GetFinalFlipY()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<5<<")"<<std::endl;}
		if (GetFinalFlipZ()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<7<<")"<<std::endl;}
		if (GetFinalFlipXY()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<1<<")"<<std::endl;}
		if (GetFinalFlipYZ()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<2<<")"<<std::endl;}
		if (GetFinalFlipXZ()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<3<<")"<<std::endl;}
		if (GetFinalFlipXYZ()==1){
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  "<<6<<")"<<std::endl;}

				BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *SPHARM.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      listFileInDir(flip ${tdir} *SPHARM.coef)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(regTemplate ${reg})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(flipTemplate ${flip})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
				//BMSShapeAnalysisModuleFile<<"      echo ('regTemplate: '${regTemplate})"<<std::endl;
				//BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate: '${flipTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;

if (GetParaOut1State()==true)
		{
		BMSShapeAnalysisModuleFile<<"    SetAppOption(ParaT.writePara 1)"<<std::endl;
		}
		BMSShapeAnalysisModuleFile<<"    Run(output ${ParaT} error)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    echo('ParaToSPHARMMesh: '${error})"<<std::endl;
		//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
		BMSShapeAnalysisModuleFile<<" DeleteFile("<<GetOutputDirectory()<<"/Template/${basename}_pp_surf_paraPhiHalf.txt)"<<std::endl;
		BMSShapeAnalysisModuleFile<<" DeleteFile("<<GetOutputDirectory()<<"/Template/${basename}_pp_surf_paraMix.txt)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  endif(${value})"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  "<<std::endl;

if (GetRegTemplateState()==false && GetFlipTemplateState()==false)//TODO
	{		
		
		
		
		if (GetTemplateMState()==true) 
		{
			if (MeanTemplateExist()==0)
			{
				BMSShapeAnalysisModuleFile<<"    Set(value 0)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    foreach(data ${testPara2})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Inc(${value} 1)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Int(${value})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    EndForEach(data)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    if (${value} < 4) "<<std::endl;
				BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *SPHARM.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      listFileInDir(flip ${tdir} *SPHARM.coef)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(regTemplate ${reg})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(flipTemplate ${flip})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      echo ('regTemplate: '${regTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate: '${flipTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      #The Template is the first data"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Set(value 0)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      foreach(data ${testPara})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        Inc(${value} 1)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        Int(${value})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      EndForEach(data)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      if (${value} < 5)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetApp(Para @ParaToSPHARMMeshCLP)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"	     SetAppOption(Para.flipTemplateFileOn 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.subdivLevel 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.subdivLevel.subdivLevel "<</*m_Parameters.*/GetSubdivLevel()<<")"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.spharmDegree 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.spharmDegree.spharmDegree "<</*m_Parameters.*/GetSPHARMDegree()<<")"<<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFileOn 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex 1)" <<std::endl;
				if (GetFinalFlipN()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<0<<")"<<std::endl;}
				if (GetFinalFlipX()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<4<<")"<<std::endl;}
				if (GetFinalFlipY()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<5<<")"<<std::endl;}
				if (GetFinalFlipZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<7<<")"<<std::endl;}
				if (GetFinalFlipXY()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<1<<")"<<std::endl;}
				if (GetFinalFlipYZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<2<<")"<<std::endl;}
				if (GetFinalFlipXZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<3<<")"<<std::endl;}
				if (GetFinalFlipXYZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<6<<")"<<std::endl;}
				
				if (GetParaOut1State()==true)
				{
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.writePara 1)"<<std::endl;
				}
				BMSShapeAnalysisModuleFile<<"        Run(output ${Para} error)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    echo('ParaToSPHARMMesh: '${error})"<<std::endl;
				//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      endif(${value})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    endif(${value})"<<std::endl;
				ComputeMean = 1;
			}
			if (MeanTemplateExist()!=0)
			{
				ComputeMean = 0;
				BMSShapeAnalysisModuleFile<<"    listFileInDir(deletelist ${SPHARMdir} ${basename}*SPHARM*) "<<std::endl;  
				BMSShapeAnalysisModuleFile<<"    Set(value 0)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    foreach(data ${testPara2})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Inc(${value} 1)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Int(${value})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    EndForEach(data)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    if (${value} < 4) "<<std::endl;
				BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *meanAll.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(regTemplate ${reg})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      listFileInDir(flip ${tdir} *SPHARM.coef)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(flipTemplate ${flip})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      echo()"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      echo ('regTemplate: '${regTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate: '${flipTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      echo()"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      #The Template is the Mean"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      listFileInDir(testPara ${SPHARMdir} ${basename}*SPHARM*)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetApp(Para @ParaToSPHARMMeshCLP)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf_tMean)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"	   SetAppOption(Para.flipTemplateFileOn 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.subdivLevel 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.subdivLevel.subdivLevel "<</*m_Parameters.*/GetSubdivLevel()<<")"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.spharmDegree 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.spharmDegree.spharmDegree "<</*m_Parameters.*/GetSPHARMDegree()<<")"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.regTemplateFileOn 1)" <<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex 1)" <<std::endl;
				if (GetFinalFlipN()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<0<<")"<<std::endl;}
				if (GetFinalFlipX()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<4<<")"<<std::endl;}
				if (GetFinalFlipY()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<5<<")"<<std::endl;}
				if (GetFinalFlipZ()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<7<<")"<<std::endl;}
				if (GetFinalFlipXY()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<1<<")"<<std::endl;}
				if (GetFinalFlipYZ()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<2<<")"<<std::endl;}
				if (GetFinalFlipXZ()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<3<<")"<<std::endl;}
				if (GetFinalFlipXYZ()==1){
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<6<<")"<<std::endl;}
			
				if (GetParaOut1State()==true)
				{
					BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.writePara 1)"<<std::endl;
				}
				BMSShapeAnalysisModuleFile<<"      Run(output ${Para} error)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    echo('ParaToSPHARMMesh: '${error})"<<std::endl;
				//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    endif(${value})"<<std::endl;
			}
		}
		if (GetTemplateMState()==false)
		{
			BMSShapeAnalysisModuleFile<<"    Set(value 0)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"    foreach(data ${testPara2})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      Inc(${value} 1)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      Int(${value})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"    EndForEach(data)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"    if (${value} < 4) "<<std::endl;
			BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *SPHARM.vtk)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      listFileInDir(flip ${tdir} *SPHARM.coef)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      set(regTemplate ${reg})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      set(flipTemplate ${flip})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      echo ('regTemplate: '${regTemplate})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate: '${flipTemplate})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      #The Template is the first data"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      Set(value 0)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      foreach(data ${testPara})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        Inc(${value} 1)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        Int(${value})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      EndForEach(data)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      if (${value} < 4)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetApp(Para @ParaToSPHARMMeshCLP)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"	     SetAppOption(Para.flipTemplateFileOn 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.subdivLevel 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.subdivLevel.subdivLevel "<<GetSubdivLevel()<<")"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.spharmDegree 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.spharmDegree.spharmDegree "<<GetSPHARMDegree()<<")"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFileOn 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex 1)" <<std::endl;
			if (GetFinalFlipN()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<0<<")"<<std::endl;}
			if (GetFinalFlipX()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<4<<")"<<std::endl;}
			if (GetFinalFlipY()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<5<<")"<<std::endl;}
			if (GetFinalFlipZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<7<<")"<<std::endl;}
			if (GetFinalFlipXY()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<1<<")"<<std::endl;}
			if (GetFinalFlipYZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<2<<")"<<std::endl;}
			if (GetFinalFlipXZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<3<<")"<<std::endl;}
			if (GetFinalFlipXYZ()==1){
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<6<<")"<<std::endl;}
			
			if (GetParaOut1State()==true)
			{
				BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.writePara 1)"<<std::endl;
			}
			BMSShapeAnalysisModuleFile<<"        Run(output ${Para} error)"<<std::endl;
			BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
			BMSShapeAnalysisModuleFile<<"    echo('ParaToSPHARMMesh: '${error})"<<std::endl;
			//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
			BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"      endif(${value})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"    endif(${value})"<<std::endl;
		}
 	}


	if ((GetRegTemplateState()==true && GetFlipTemplateState()==true) ||(GetRegTemplateState()==true && GetFlipTemplateState()==false)||(GetRegTemplateState()==false && GetFlipTemplateState()==true))

	{


BMSShapeAnalysisModuleFile<<"    Set(value 0)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    foreach(data ${testPara2})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Inc(${value} 1)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      Int(${value})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    EndForEach(data)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"    if (${value} < 4) "<<std::endl;
if(GetTemplateMState()==true &&  MeanTemplateExist()!=0 )
{BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *meanAll.vtk)"<<std::endl;}
else {BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *SPHARM.vtk)"<<std::endl;}
				BMSShapeAnalysisModuleFile<<"      listFileInDir(flip ${tdir} *SPHARM.coef)"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(regTemplate ${reg})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"      set(flipTemplate ${flip})"<<std::endl;
				BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
if (GetRegTemplateState()==false  )
				{BMSShapeAnalysisModuleFile<<"      echo ('regTemplate: '${regTemplate})"<<std::endl;}
if (GetFlipTemplateState()==false) 
				{BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate: '${flipTemplate})"<<std::endl;}
				BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;



		BMSShapeAnalysisModuleFile<<"  #The Template is choosed by the user"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
		if (GetRegTemplateState()==true)
			{BMSShapeAnalysisModuleFile<<"  echo ('regTemplate: '"<<GetRegTemplate()<<")"<<std::endl;}
		if (GetFlipTemplateState()==true)
			{BMSShapeAnalysisModuleFile<<"  echo ('flipTemplate: '"<<GetFlipTemplate()<<")"<<std::endl;}
		BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  Set(value 0)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  foreach(data ${testPara})"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    Inc(${value} 1)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    Int(${value})"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  EndForEach(data)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  if (${value} < 4)"<<std::endl;
		
		BMSShapeAnalysisModuleFile<<"SetApp(Para @ParaToSPHARMMeshCLP)" <<std::endl;
		BMSShapeAnalysisModuleFile<<"SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)" <<std::endl;
BMSShapeAnalysisModuleFile<<"SetAppOption(Para.inParaFile  ${SPHARMdir}${basename}_pp_para.vtk)" <<std::endl;

		if(MeanTemplateExist()!=0 &&  MeanTemplateExist()!=0 ) 
			{ BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf_tMean)"<<std::endl;}
		else
			{BMSShapeAnalysisModuleFile<<"SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf)" <<std::endl;}
		BMSShapeAnalysisModuleFile<<"SetAppOption(Para.subdivLevel 1)" <<std::endl;
		BMSShapeAnalysisModuleFile<<"SetAppOption(Para.subdivLevel.subdivLevel "<<GetSubdivLevel()<<")" <<std::endl;
		BMSShapeAnalysisModuleFile<<"SetAppOption(Para.spharmDegree 1)" <<std::endl;
		BMSShapeAnalysisModuleFile<<"SetAppOption(Para.spharmDegree.spharmDegree "<<GetSPHARMDegree()<<")" <<std::endl;
		if (GetFlipTemplateState()==true){
			BMSShapeAnalysisModuleFile<<"SetAppOption(Para.flipTemplateFileOn 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"SetAppOption(Para.flipTemplateFile.flipTemplateFile "<<GetFlipTemplate()<<")" <<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFileOn 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"<<std::endl;}
		if (GetRegTemplateState()==true){
			BMSShapeAnalysisModuleFile<<"SetAppOption(Para.flipTemplateFileOn 1)" <<std::endl;
			BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})"<<std::endl;
			BMSShapeAnalysisModuleFile<<"SetAppOption(Para.regTemplateFileOn 1)" <<std::endl;
			if(MeanTemplateExist()!=0 &&  MeanTemplateExist()!=0){BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"<<std::endl;}
			else{BMSShapeAnalysisModuleFile<<"SetAppOption(Para.regTemplateFile.regTemplateFile "<<GetRegTemplate()<<")" <<std::endl;}}
		BMSShapeAnalysisModuleFile<<"SetAppOption(Para.finalFlipIndex 1)" <<std::endl;
		if (GetFinalFlipN()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<0<<")"<<std::endl;}
		if (GetFinalFlipX()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<4<<")"<<std::endl;}
		if (GetFinalFlipY()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<5<<")"<<std::endl;}
		if (GetFinalFlipZ()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<7<<")"<<std::endl;}
		if (GetFinalFlipXY()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<1<<")"<<std::endl;}
		if (GetFinalFlipYZ()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<2<<")"<<std::endl;}
		if (GetFinalFlipXZ()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<3<<")"<<std::endl;}
		if (GetFinalFlipXYZ()==1){
			BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<6<<")"<<std::endl;}



		if (GetParaOut1State()==true)	{
			BMSShapeAnalysisModuleFile<<"    SetAppOption(Para.writePara 1)"<<std::endl;}
		BMSShapeAnalysisModuleFile<<"    Run(output ${Para} error)"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
		BMSShapeAnalysisModuleFile<<"    echo('ParaToSPHARMMesh: '${error})"<<std::endl;
		//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;}

	
	BMSShapeAnalysisModuleFile<<"Set(paraPhiHalf "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surf_paraPhiHalf.txt)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"FileExists(testParaPhiHalf ${paraPhiHalf})"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"If( ${testParaPhiHalf} == 1 )"<<std::endl;
	BMSShapeAnalysisModuleFile<<"echo(${paraPhiHalf})"<<std::endl;
	BMSShapeAnalysisModuleFile<<" DeleteFile(${paraPhiHalf})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"EndIf(${testParaPhiHalf})"<<std::endl;


	
	BMSShapeAnalysisModuleFile<<"Set(paraMix "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surf_paraMix.txt)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"FileExists(testParaMix ${paraMix})"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"If( ${testParaMix} == 1 )"<<std::endl;
	BMSShapeAnalysisModuleFile<<" DeleteFile(${paraMix})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"EndIf(${testParaMix})"<<std::endl;
	//Write the Ouput File
	BMSShapeAnalysisModuleFile<<"  #Write Ouput File"<<std::endl;
	
	for(unsigned int i=0;i<OutputFileHeaders.size();i++)
	{
		BMSShapeAnalysisModuleFile<<"  GetParam(header"<<i<<" ${Header"<<i<<"} ${count})"<<std::endl;
		BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${header"<<i<<"} ',' )"<<std::endl;
	}
	BMSShapeAnalysisModuleFile<<"  Inc(${count} 1)"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  Set(postProcessFile "<<GetOutputDirectory()<<"/Mesh/PostProcess/${basename}_pp"<<GetVolumeFileExtension()<<")"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  FileExists(testPostProcessFile ${postProcessFile})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  If(${testPostProcessFile} == 1)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${postProcessFile} ',')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  EndIf(${testPostProcessFile})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  If(${testPostProcessFile} == 0)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} 'none,')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  EndIf(${testPostProcessFile})"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  Set(originSurf "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_para.vtk)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${originSurf} ',')"<<std::endl;

	BMSShapeAnalysisModuleFile<<"  Set(surfSPHARM "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surfSPHARM.vtk)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${surfSPHARM} ',')"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  Set(surfSPHARMcoef "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surfSPHARM.coef)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${surfSPHARMcoef} ',')"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  Set(surfSPHARM_ellalign "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surfSPHARM_ellalign.vtk)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${surfSPHARM_ellalign} ',')"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  Set(surfSPHARM_ellalignCoef "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surfSPHARM_ellalign.coef)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${surfSPHARM_ellalignCoef} ',')"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"  Set(procalignFile "<<GetOutputDirectory()<<"/Mesh/SPHARM/${basename}_pp_surfSPHARM_procalign.vtk)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  FileExists(testProcalignFile ${procalignFile})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  If(${testProcalignFile} == 1)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} ${procalignFile} )"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  EndIf(${testProcalignFile})"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  If(${testProcalignFile} == 0)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} 'none')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  EndIf(${testProcalignFile})"<<std::endl;



	BMSShapeAnalysisModuleFile<<"  appendFile(${OutputFile} '\\n' )"<<std::endl;
	BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
	BMSShapeAnalysisModuleFile<<"EndForEach(case)"<<std::endl;
	BMSShapeAnalysisModuleFile<<""<<std::endl<<std::endl;
	
	/*BMSShapeAnalysisModuleFile<<" Set(BMSFile ${BMSdir}"<<"/ShapeAnalysisModule.bms"<<")"<<std::endl;
	BMSShapeAnalysisModuleFile<<" CopyFile("<<GetBMSShapeAnalysisModuleFile()<<" ${BMSdir})"<<std::endl;
	BMSShapeAnalysisModuleFile<<" DeleteFile("<<GetBMSShapeAnalysisModuleFile()<<")"<<std::endl;*/
	
	return;
}  

//Write BMS file for mean computation
void ShapeAnalysisModuleComputation::WriteBMSShapeAnalysisModuleFile2()
{

std::ofstream BMSShapeAnalysisModuleFile(GetBMSShapeAnalysisModuleFile2());

BMSShapeAnalysisModuleFile<<"#----- Shape Analysis Computing ParaToSPHARMMesh with mean as input -----"<<std::endl;
BMSShapeAnalysisModuleFile<<"#------------------------------------------------------------------------"<<std::endl;
BMSShapeAnalysisModuleFile<<"#------------------------------------------------------------------------"<<std::endl;
BMSShapeAnalysisModuleFile<<"#------------------------------------------------------------------------"<<std::endl;

BMSShapeAnalysisModuleFile<<"echo()"<<std::endl;
BMSShapeAnalysisModuleFile<<"echo('Computing with mean as a template')"<<std::endl;
BMSShapeAnalysisModuleFile<<"echo()"<<std::endl;
BMSShapeAnalysisModuleFile<<"set (OrigCasesList '"<<GetNthDataListValue(1,GetColumnVolumeFile())<<"')"<<std::endl;

int DataNumber;
for (DataNumber = 2; DataNumber <=GetDataNumber(); DataNumber++)
{
BMSShapeAnalysisModuleFile<<"set (OrigCasesList ${OrigCasesList} '"<<GetNthDataListValue(DataNumber,GetColumnVolumeFile())<<"')"<<std::endl;
}
BMSShapeAnalysisModuleFile<<""<<std::endl;

BMSShapeAnalysisModuleFile<<"set(SPHARMdir '"<<GetOutputDirectory()<<"/Mesh/SPHARM/')"<<std::endl;
BMSShapeAnalysisModuleFile<<"set(tdir '"<<GetOutputDirectory()<<"/Template/')"<<std::endl;

BMSShapeAnalysisModuleFile<<"ForEach(case ${OrigCasesList})"<<std::endl;
BMSShapeAnalysisModuleFile<<"  #Extract basename"<<std::endl;
BMSShapeAnalysisModuleFile<<"  GetFilename(basename ${case} NAME_WITHOUT_EXTENSION)"<<std::endl;
BMSShapeAnalysisModuleFile<<"  echo()"<<std::endl;
BMSShapeAnalysisModuleFile<<"  echo('Case: '${case})"<<std::endl;
BMSShapeAnalysisModuleFile<<"  "<<std::endl;
BMSShapeAnalysisModuleFile<<"  listFileInDir(testPara ${SPHARMdir} *${basename}*SPHARM*)"<<std::endl;
BMSShapeAnalysisModuleFile<<"  listFileInDir(testPara2 ${SPHARMdir} *${basename}*MeanSPHARM*)"<<std::endl;
BMSShapeAnalysisModuleFile<<"  listFileInDir(testtemp ${tdir}  *SPHARM*)"<<std::endl;
BMSShapeAnalysisModuleFile<<"  listFileInDir(testtemp2 ${tdir} *${basename}*SPHARM*)"<<std::endl;



BMSShapeAnalysisModuleFile<<"    if (${value} < 4) "<<std::endl;
BMSShapeAnalysisModuleFile<<"      listFileInDir(reg ${tdir} *meanAll.vtk)"<<std::endl;
BMSShapeAnalysisModuleFile<<"      set(regTemplate ${reg})"<<std::endl;
BMSShapeAnalysisModuleFile<<"      listFileInDir(flip ${tdir} *SPHARM.coef)"<<std::endl;
BMSShapeAnalysisModuleFile<<"      set(flipTemplate ${flip})"<<std::endl;
BMSShapeAnalysisModuleFile<<"      echo()"<<std::endl;
BMSShapeAnalysisModuleFile<<"      echo ('regTemplate : '${regTemplate})"<<std::endl;
BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate : '${flipTemplate})"<<std::endl;
BMSShapeAnalysisModuleFile<<"      echo()"<<std::endl;
BMSShapeAnalysisModuleFile<<"      #The Template is the Mean"<<std::endl;
BMSShapeAnalysisModuleFile<<"      listFileInDir(testPara ${SPHARMdir} ${basename}*SPHARM*)"<<std::endl;
BMSShapeAnalysisModuleFile<<"       SetApp(Para @ParaToSPHARMMeshCLP)"<<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"<<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"<<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf_tMean)"<<std::endl;
BMSShapeAnalysisModuleFile<<"	     SetAppOption(Para.flipTemplateFileOn 1)" <<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})"<<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.subdivLevel 1)" <<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.subdivLevel.subdivLevel "<<GetSubdivLevel()<<")"<<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.spharmDegree 1)" <<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.spharmDegree.spharmDegree "<<GetSPHARMDegree()<<")"<<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFileOn 1)" <<std::endl;
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"<<std::endl;
BMSShapeAnalysisModuleFile<<"SetAppOption(Para.finalFlipIndex 1)" <<std::endl;
if (GetFinalFlipN()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<0<<")"<<std::endl;}
if (GetFinalFlipX()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<4<<")"<<std::endl;}
if (GetFinalFlipY()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<5<<")"<<std::endl;}
if (GetFinalFlipZ()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<7<<")"<<std::endl;}
if (GetFinalFlipXY()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<1<<")"<<std::endl;}
if (GetFinalFlipYZ()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<2<<")"<<std::endl;}
if (GetFinalFlipXZ()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<3<<")"<<std::endl;}
if (GetFinalFlipXYZ()==1){
BMSShapeAnalysisModuleFile<<"        SetAppOption(Para.finalFlipIndex.finalFlipIndex "<<6<<")"<<std::endl;}
 
if (GetParaOut1State()==true)
 {
BMSShapeAnalysisModuleFile<<"      SetAppOption(Para.writePara 1)"<<std::endl;
 }
BMSShapeAnalysisModuleFile<<"      Run(output ${Para} error)"<<std::endl;
BMSShapeAnalysisModuleFile<<"  if(${error} != '')"<<std::endl;
BMSShapeAnalysisModuleFile<<"    echo('ParaToSPHARMMesh: '${error})"<<std::endl;
//BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
BMSShapeAnalysisModuleFile<<"  endif(${error})"<<std::endl;
BMSShapeAnalysisModuleFile<<"    endif(${value})"<<std::endl;

BMSShapeAnalysisModuleFile<<"EndForEach(${case})"<<std::endl;

BMSShapeAnalysisModuleFile<<"set(BMSdir '"<<GetOutputDirectory()<<"/BatchMake_Scripts')"<<std::endl;
BMSShapeAnalysisModuleFile<<" Set(BMSFile ${BMSdir}"<<"/ShapeAnalysisModule_MeanAsTemplate.bms"<<")"<<std::endl;
BMSShapeAnalysisModuleFile<<" CopyFile("<<GetBMSShapeAnalysisModuleFile()<<" ${BMSdir})"<<std::endl;
BMSShapeAnalysisModuleFile<<" DeleteFile("<<GetBMSShapeAnalysisModuleFile()<<")"<<std::endl;

}

//Overwrite (delete) the different files 
void ShapeAnalysisModuleComputation::OverWrite()
{
	int DataNumber=GetDataNumber();
	char **Base_Files; 
	Base_Files = new char *[DataNumber];
	

	for(int i=0;i<DataNumber;i++)
	{	
		Base_Files[i] = new char[512];
		std::strcpy(Base_Files[i],GetAllFilesName(i));
	}


	//Overwrite files created by SegPostProcess
	if( GetOverwriteSegPostProcess() )
	{
		char **PreProcessing_Files; 
		PreProcessing_Files = new char *[DataNumber];

	 	char *dir;
		dir=new  char[512] ;
		std::strcpy(dir, GetOutputDirectory());
		std::strcat(dir, "/Mesh/PostProcess/");

		for (int i = 0; i<DataNumber;i++)
		{
			PreProcessing_Files[i] = new char[512];  
			std::strcpy(PreProcessing_Files[i],dir);
			std::strcat(PreProcessing_Files[i],Base_Files[i]);
			std::strcat(PreProcessing_Files[i],"_pp.gipl.gz");
		}
		
		bool DirectoryEmpty=DirectoryIsEmpty(dir);
		if(!DirectoryEmpty)
		{
			cout<<"OverWrite SegPostProcess"<<endl;
			for (int k=0; k<DataNumber; k++)
			{	
				remove(PreProcessing_Files[k]);
				//cerr<<"Error deleting file "<<PreProcessing_Files[k]<<endl;
					
			}
		}
		else{ cerr<<"You have to compute SegPostProcess before overwriting!!"<<endl;}
	}
	
	//Overwrite files created by GenParaMesh
	if( GetOverwriteGenParaMesh() )
	{
		char **GenParaMesh_Files; 
		GenParaMesh_Files = new char *[2*DataNumber];

	 	char *dir;
		dir=new  char[512] ;
		std::strcpy(dir, GetOutputDirectory());
		std::strcat(dir, "/Mesh/SPHARM/");

		for (int i = 0; i<DataNumber;i++)
		{
			GenParaMesh_Files[i] = new char[512];  
			std::strcpy(GenParaMesh_Files[i],dir);
			std::strcat(GenParaMesh_Files[i],Base_Files[i]);
			std::strcat(GenParaMesh_Files[i],"_pp_para.vtk");

			GenParaMesh_Files[i+DataNumber] = new char[512];  
			std::strcpy(GenParaMesh_Files[i+DataNumber],dir);
			std::strcat(GenParaMesh_Files[i+DataNumber],Base_Files[i]);
			std::strcat(GenParaMesh_Files[i+DataNumber],"_pp_surf.vtk");
		}
		
		bool DirectoryEmpty=DirectoryIsEmpty(dir);
		if(!DirectoryEmpty)
		{
			cout<<"OverWrite GenParaMesh"<<endl;
			for (int k=0; k<2*DataNumber; k++)
			{
				remove(GenParaMesh_Files[k]);
				//cerr<<"Error deleting file "<<GenParaMesh_Files[k]<<endl;
			}
		}
		else{ cerr<<"You have to compute GenParaMesh before overwriting!!"<<endl;}
	}

	//Overwrite files created by ParaToSPHARMMesh
	if( GetOverwriteParaToSPHARMMesh() )
	{
		char **ParaToSPHARMMesh_Files; 
		ParaToSPHARMMesh_Files = new char *[7*DataNumber];

	 	char *dir;
		dir=new  char[512] ;
		std::strcpy(dir, GetOutputDirectory());
		std::strcat(dir, "/Mesh/SPHARM/");

		for (int i = 0; i<DataNumber;i++)
		{
			ParaToSPHARMMesh_Files[i] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i],dir);
			std::strcat(ParaToSPHARMMesh_Files[i],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i],"_pp_surfSPHARM.vtk");

			ParaToSPHARMMesh_Files[i+DataNumber] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i+DataNumber],dir);
			std::strcat(ParaToSPHARMMesh_Files[i+DataNumber],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i+DataNumber],"_pp_surfSPHARM_ellalign.vtk");

			ParaToSPHARMMesh_Files[i+2*DataNumber] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i+2*DataNumber],dir);
			std::strcat(ParaToSPHARMMesh_Files[i+2*DataNumber],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i+2*DataNumber],"_pp_surfSPHARM.coef");

			ParaToSPHARMMesh_Files[i+3*DataNumber] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i+3*DataNumber],dir);
			std::strcat(ParaToSPHARMMesh_Files[i+3*DataNumber],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i+3*DataNumber],"_pp_surfSPHARM_ellalign.coef");

			ParaToSPHARMMesh_Files[i+4*DataNumber] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i+4*DataNumber],dir);
			std::strcat(ParaToSPHARMMesh_Files[i+4*DataNumber],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i+4*DataNumber],"_pp_surfSPHARM_procalign.vtk");

			ParaToSPHARMMesh_Files[i+5*DataNumber] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i+5*DataNumber],dir);
			std::strcat(ParaToSPHARMMesh_Files[i+5*DataNumber],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i+5*DataNumber],"_pp_surf_paraPhi.txt");

			ParaToSPHARMMesh_Files[i+6*DataNumber] = new char[512];  
			std::strcpy(ParaToSPHARMMesh_Files[i+6*DataNumber],dir);
			std::strcat(ParaToSPHARMMesh_Files[i+6*DataNumber],Base_Files[i]);
			std::strcat(ParaToSPHARMMesh_Files[i+6*DataNumber],"_pp_surf_paraTheta.txt");

// 			ParaToSPHARMMesh_Files[i+7*DataNumber] = new char[512];  
// 			std::strcpy(ParaToSPHARMMesh_Files[i+7*DataNumber],dir);
// 			std::strcat(ParaToSPHARMMesh_Files[i+7*DataNumber],Base_Files[i]);
// 			std::strcat(ParaToSPHARMMesh_Files[i+7*DataNumber],"_pp_surf_tMeanSPHARM.vtk");
// 
// 			ParaToSPHARMMesh_Files[i+8*DataNumber] = new char[512];  
// 			std::strcpy(ParaToSPHARMMesh_Files[i+8*DataNumber],dir);
// 			std::strcat(ParaToSPHARMMesh_Files[i+8*DataNumber],Base_Files[i]);
// 			std::strcat(ParaToSPHARMMesh_Files[i+8*DataNumber],"_pp_surf_tMeanSPHARM_ellalign.vtk");
// 
// 			ParaToSPHARMMesh_Files[i+9*DataNumber] = new char[512];  
// 			std::strcpy(ParaToSPHARMMesh_Files[i+9*DataNumber],dir);
// 			std::strcat(ParaToSPHARMMesh_Files[i+9*DataNumber],Base_Files[i]);
// 			std::strcat(ParaToSPHARMMesh_Files[i+9*DataNumber],"_pp_surf_tMeanSPHARM.coef");

// 			ParaToSPHARMMesh_Files[i+10*DataNumber] = new char[512];  
// 			std::strcpy(ParaToSPHARMMesh_Files[i+10*DataNumber],dir);
// 			std::strcat(ParaToSPHARMMesh_Files[i+10*DataNumber],Base_Files[i]);
// 			std::strcat(ParaToSPHARMMesh_Files[i+10*DataNumber],"_pp_surf_tMeanSPHARM_ellalign.coef");
		
		}
		
		//check if the directory is empty
		bool DirectoryEmpty=DirectoryIsEmpty(dir);
		if(!DirectoryEmpty)
		{
			cout<<"OverWrite ParaToSPARHMMesh"<<endl;
			for (int k=0; k<7*DataNumber; k++)
			{
				remove(ParaToSPHARMMesh_Files[k]);
				//cerr<<"Error deleting file "<<ParaToSPHARMMesh_Files[k]<<endl;
			}
		}

		else{ cerr<<"You have to compute ParaToSPHARMMesh before overwriting!!"<<endl;}

		
		// Overwrite the template file;
		char *dirTemplate=NULL;
		dirTemplate=new  char[512] ;
		std::strcpy(dirTemplate, GetOutputDirectory());
		std::strcat(dirTemplate, "/Template/");

		bool templateDirectoryEmpty=DirectoryIsEmpty(dirTemplate);
		if(!templateDirectoryEmpty)
		{
			cout<<"OverWrite Templates"<<endl;
			int length=strlen(GetOutputDirectory());
			length=length+10;
	
			DIR *pdir = NULL;
			struct dirent *pent;
			pdir = opendir (dirTemplate);
			while ((pent=readdir(pdir)))
			{
				char *file=NULL;
				file= new char[512];
				strcpy(file,dirTemplate);
				strcat(file, pent->d_name);
				if(file[length]!='.')
				{	
					remove(file);
					//cerr<<"Error deleting file "<<file<<endl;
				}
			}
		}
	}
}

//Compute the mean of all the files
void ShapeAnalysisModuleComputation::ComputationMean()
{
	std::cout<<"Computing Mean"<<std::endl;
	int DataNumber=GetDataNumber();
	char **Base_Files; 
	Base_Files = new char *[DataNumber];

	int nbPoints=0;
	bool initialize=true;
vtkPolyData *poly= vtkPolyData::New();
	for(int i=0;i<DataNumber;i++)
	{	
		Base_Files[i] = new char[512];
		std::strcpy(Base_Files[i], "/");
		std::strcpy(Base_Files[i], GetOutputDirectory());
		std::strcat(Base_Files[i], "/Mesh/SPHARM/");
		std::strcat(Base_Files[i],GetAllFilesName(i));
		std::strcat(Base_Files[i],"_pp_surfSPHARM_procalign.vtk");

		//read a vtk file
		vtkPolyDataReader *meshin = vtkPolyDataReader::New();
		meshin->SetFileName(Base_Files[i]);
		
		try{
			meshin->Update();		
		}
		catch(...)
		{
			std::cout << "Cannot open file: " << Base_Files[i] << std::endl;
			return ;
		}
	
		
		poly=meshin->GetOutput();
	
		//Get number of points
		vtkIdType idNumPointsInFile=poly->GetNumberOfPoints();
		nbPoints=idNumPointsInFile;
	

		vtkPoints * pts;
		pts=poly->GetPoints();

		// add all values of the coordinates
		if(initialize==true)
		{	
			for( int i = 0; i < nbPoints; i++)
			{	
				m_meanX.push_back(0);
				m_meanY.push_back(0);
				m_meanZ.push_back(0);
			}
			initialize=false;
		}

		for( int i = 0; i < nbPoints; i++)
   		{	
			double *p;
			p=pts->GetPoint(i);
		
			m_meanX[i]=m_meanX[i]+p[0];
 			m_meanY[i]=m_meanY[i]+p[1];
 			m_meanZ[i]=m_meanZ[i]+p[2];		
		}
	}

	for(int i = 0; i < nbPoints; i++)
   	{	
		m_meanX[i]=m_meanX[i]/DataNumber;
		m_meanY[i]=m_meanY[i]/DataNumber;
		m_meanZ[i]=m_meanZ[i]/DataNumber;	
	}


	char m_PolyFile[512];

	std::strcpy(m_PolyFile,GetOutputDirectory());
	std::strcat(m_PolyFile,"/Template/");
	std::strcat(m_PolyFile,"polyAll.vtk");

	vtkCellArray* vtkcells =poly->GetPolys();
	vtkPolyData *polydata = vtkPolyData::New();
	polydata-> SetPolys(vtkcells);
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName (m_PolyFile);
	writer->SetInput ( polydata);
	writer->Write();


		WriteMeanFile(nbPoints);



}

//Write the mean file use as a template for ParaToSPHARMMesh
void ShapeAnalysisModuleComputation::WriteMeanFile(int nbPoints)
{
	char m_MeanFile[512];char m_PolyFile[512];
	std::strcpy(m_MeanFile,GetOutputDirectory());
	std::strcat(m_MeanFile,"/Template/");
	std::strcpy(m_PolyFile,m_MeanFile);
	std::strcat(m_MeanFile,"meanAll.vtk");
	std::strcat(m_PolyFile,"polyAll.vtk");

	std::string line;
	int nbline=0;

	std::ifstream tmpPolyFile(m_PolyFile);
	std::ofstream meanAll(m_MeanFile);
	meanAll<<"# vtk DataFile Version 3.0"<<std::endl;
	meanAll<<"vtk output"<<std::endl;
	meanAll<<"ASCII"<<std::endl;
	meanAll<<"DATASET POLYDATA"<<std::endl;
	meanAll<<"POINTS "<<nbPoints<<" float"<<std::endl;

	for(unsigned int i=0;i<m_meanX.size();i+=3)
	{
		meanAll<<m_meanX[i]<<" "<<m_meanY[i]<<" "<<m_meanZ[i]<<" "<<m_meanX[i+1]<<" "<<m_meanY[i+1]<<" "<<m_meanZ[i+1]<<" "<<m_meanX[i+2]<<" "<<m_meanY[i+2]<<" "<<m_meanZ[i+2]<<endl;
	}


	while(getline(tmpPolyFile, line)) 
	{	
		if(nbline>4)
		{
			meanAll<<line<<std::endl;
		}
		nbline++;
	}
	remove( m_PolyFile );


}

void ShapeAnalysisModuleComputation::GetRandomNum(int min, int max)
{
  //int min, max;

  srand(getpid()+time(NULL));

  int nofNum = max - min + 1;
  //std::cout << nofNum << " " << max << " " << min << std::endl;

  // Create a random permutation
  int index, temp;
  for (int i = 0; i < nofNum; i++)
    m_permutations[i] = i;
  
  for (int i = 1; i < nofNum; i++) {
    index = i + (rand() % (nofNum - i));
    temp = m_permutations[i];
    m_permutations[i] = m_permutations[index];
    m_permutations[index] = temp;
  }

  for (int i = 0; i < nofNum; i++)
    m_permutations[i]++;
  // random swaps

  //print perm
  //cout << "permutation : ";
  /*for (int i = 0; i < nofNum; i++)
    if (i == 0)
	{
		//m_permutations[i] = min + m_permutations[i];
		cout << min + m_permutations[i];
      	}
    else
      cout << " " << min + m_permutations[i];
 cout << endl;*/
}


void ShapeAnalysisModuleComputation::RunParticlesModule()
{
	std::cout<<" "<<std::endl;
	std::cout<<"Run the Particle Correspondence"<<std::endl;
	std::string particlesdirectory;
	std::string outputdir; outputdir.append(GetOutputDirectory());
	particlesdirectory.append(outputdir);
	particlesdirectory.append("/ParticleCorrespondence");
	itksys::SystemTools::MakeDirectory(particlesdirectory.c_str());

	std::string csvdirectory;
	csvdirectory.append(outputdir);
	csvdirectory.append("/OutputGroupFile/ShapeAnalysisModule_OutputFileVersion1.csv");


	//command line 
	std::vector<const char*> args;
	args.push_back("ParticleModule" );
	args.push_back("--columMeshFile" );
	if(GetUseProcalign()){args.push_back("9");}
	else{args.push_back("5");}//Original Space
	args.push_back("--sx" );
	args.push_back(Convert_Double_To_CharArray(GetEnforcedSpaceX()));
	args.push_back("--sy" );
	args.push_back(Convert_Double_To_CharArray(GetEnforcedSpaceY()));
	args.push_back("--sz" );
	args.push_back(Convert_Double_To_CharArray(GetEnforcedSpaceZ()));
	args.push_back(csvdirectory.c_str() );
	args.push_back(particlesdirectory.c_str());
	args.push_back(0);

	for( unsigned int k =0; k<args.size()-1;k++)
	{std::cout<<args.at(k)<<" ";}
	std::cout<<" "<<std::endl;std::cout<<" "<<std::endl;

	
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
			time (&end);
			cout<<"(processing since "<<difftime (end,start)<<" seconds) \r"<<flush;
			timeout = 0.05;   	
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



void ShapeAnalysisModuleComputation::CreateMrmlParticle()
{
	int nbShapesPerMRML= GetHorizontalGridPara() * GetVerticalGridPara();
	
	for( int colormap=0;colormap<2;colormap++)
	{

		int DataNumber=GetDataNumber();
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
			std::vector<std::string> NameColormap;
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
			mrmlfile.append("/ParticleCorrespondence/MRML/ParticleModuleMRMLscene_");
			if(colormap==0){mrmlfile.append("Phi");}
			else{mrmlfile.append("Theta");}
	
	
	
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
			std::cout<<"mrml "<<mrmlfile<<std::endl;	
			argsMRML.push_back("CreateMRML");   
			argsMRML.push_back(mrmlfile.c_str());
	
			
			for(int k=0;k<DataNumberPerMRML;k++)
			{
				size_t found,found2;
				std::string relativeName,NameFile,TransName;
				vector<double>Dims;

				std::string file;
				file.append ( (string) GetPostCorrespondenceFiles(k));
	
				found=file.find_last_of("/\\");
				relativeName.append(file.substr(found+1));
				
				Namevtk.push_back(relativeName);
				relativeName.insert(0,"../Corresponding_Meshes/");  

				RelativePathToVTK.push_back(relativeName);
				found2=(Namevtk.back()).find_last_of(".\\");

				TransName.append((Namevtk.back()).substr(0,found2));
				Name.push_back(TransName);

				TransName.append("Trans");
				NameTrans.push_back(TransName);

	
				std::string file0;
				file0.append ( (string) GetPostCorrespondenceFiles(0));
				SetImageDimensions((char *)file0.c_str());
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
				tmp_transformfile.append("./TransformFiles/transformCM");  
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
	
	
				//add color map
				std::string tmp_NameColormap;
				if(colormap==0){tmp_NameColormap.append("../../Mesh/SPHARM/customLUT_Color_Map_Phi.txt");}
				else{tmp_NameColormap.append("../../Mesh/SPHARM/customLUT_Color_Map_Theta.txt");}
				NameColormap.push_back(tmp_NameColormap);
				argsMRML.push_back("-as" );
				if(colormap==0){argsMRML.push_back("Color_Map_Phi" );}
				else{argsMRML.push_back("Color_Map_Theta" );}
				argsMRML.push_back("-cc" );
				argsMRML.push_back((NameColormap.back()).c_str() );
	
	
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
	
			/*for (int l =0; l<argsMRML.size();l++)
			{std::cout << argsMRML.at(l)<< " " ;}*/
	
			argsMRML.push_back(0);
			//itk sys parameters
			int length;
			time_t start,end;
			time (&start);
		
			double timeout = 0.05;
			int result;
			char* dataitk = NULL;
		
			itksysProcess* gp = itksysProcess_New();
			itksysProcess_SetCommand(gp, &*argsMRML.begin());
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
					std::cerr<<"Error: Could not run " << argsMRML[0]<<":\n";
					std::cerr<<itksysProcess_GetErrorString(gp)<<"\n";
					std::cout<<"Error: Could not run " << argsMRML[0]<<":\n";
					std::cout<<itksysProcess_GetErrorString(gp)<<"\n";
				} break;
				case itksysProcess_State_Exception:
				{
					std::cerr<<"Error: "<<argsMRML[0]<<" terminated with an exception: "<<itksysProcess_GetExceptionString(gp)<<"\n";
					std::cout<<"Error: "<<argsMRML[0]<<" terminated with an exception: "<<itksysProcess_GetExceptionString(gp)<<"\n";
				} break;
				case itksysProcess_State_Starting:
				case itksysProcess_State_Executing:
				case itksysProcess_State_Expired:
				case itksysProcess_State_Killed:
				{
				// Should not get here.
				std::cerr<<"Unexpected ending state after running "<<argsMRML[0]<<std::endl;
				std::cout<<"Unexpected ending state after running "<<argsMRML[0]<<std::endl;
				} break;
			}
			itksysProcess_Delete(gp);  
		}
	}
}

