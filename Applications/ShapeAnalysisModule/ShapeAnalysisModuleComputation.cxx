#include "ShapeAnalysisModuleComputation.h"

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
	SetAllFilesName();
	OverWrite();

	WriteBMSShapeAnalysisModuleFile();
 	ExecuteBatchMake(GetBMSShapeAnalysisModuleFile());

	if(GetTemplateMState()==true)
	{
		ComputationMean();

		WriteBMSShapeAnalysisModuleFile2();
		ExecuteBatchMake(GetBMSShapeAnalysisModuleFile2());
	}

	//execute MeshMath external application
	for(int i=0;i<GetDataNumber();i++)
	{
 		ExecuteMeshMath(i,"phi");
		ExecuteMeshMath(i,"theta");
	}

 	WriteBMSMRMLScene();
 	ExecuteBatchMake(GetBMSShapeAnalysisModuleMRMLFile()); 
	SetBMSShapeAnalysisModuleFile(true);

	ModifyCSV();


	std::cout<<"Computing ShapeAnalysisModule: Done!"<<std::endl<<std::endl;
  
  	return;
}

//Execute MeshMath to write a KWM scalar field (1D) into a PolyData Field Data Scalar to visualize in Slicer3
void ShapeAnalysisModuleComputation::ExecuteMeshMath(int numData, char * scalar)
{

//execute MeshMath for each volume file: Original, Ellalign, Procalign
	for(int j=0;j<3;j++)
	{
	
		std::vector<const char*> args;  
		char* data = NULL;
		int length;
		double timeout = 0.05;
		int result;
		char *fileType=NULL;
		
		args.push_back("MeshMath");
		
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
	
		args.push_back(fileType);
		args.push_back(fileType);
		args.push_back("-KWMtoPolyData");
		
		if(scalar=="phi")
		{  
			args.push_back(GetAllPhiFiles(0));  
			args.push_back("Color Map Phi");
		}
		if(scalar=="theta")
		{ 
			args.push_back(GetAllThetaFiles(0));
			args.push_back("Color Map Theta");
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


// Create BMS File to compute the SPHARM pipeline
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleFile(bool changeDirectory)
{
	std::strcpy(m_BMSShapeAnalysisModuleFile, GetOutputDirectory());
	
	if(changeDirectory==false)
		std::strcat(m_BMSShapeAnalysisModuleFile, "/");

	else std::strcat(m_BMSShapeAnalysisModuleFile, "/BatchMake_Scripts/");

	if (GetRandomizeInputs())
	{
		int pID=getpid();
		sprintf(m_BMSShapeAnalysisModuleFile,"ShapeAnalysisModule_id%d.bms",pID);
		std::cout << " " << std::endl;
	}	
	else
		std::strcat(m_BMSShapeAnalysisModuleFile, "ShapeAnalysisModule.bms");

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
	std::strcat(m_OutputFile, "/Output/");
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
void ShapeAnalysisModuleComputation::WriteBMSMRMLScene()
{
 std::cout<<"\tWritting Mrml"<<std::endl<<endl;
	int SnapShotNumber = SetNbSnapShot();

	SetAllFilesName();
	std::ofstream BMSShapeAnalysisModuleMRML(m_BMSShapeAnalysisModuleMRLMFile);


	BMSShapeAnalysisModuleMRML<<"  echo(\t'Writing MRML scripts...')"<<std::endl;
	for(int count=0;count<4;count++)
	{
		BMSShapeAnalysisModuleMRML<<"  # Script to create a MRML scene"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(countSnapNb 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  MakeDirectory("<<GetOutputDirectory()<<"/MRML)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" Set(FileDir '"<<GetOutputDirectory()<<"/Mesh/SPHARM')"<<std::endl;
	
	// create a MRML file with all the original SPHARM files
	if(count==0)
	{
		char *firstFile1;
		firstFile1 = new char[512];
		std::strcpy(firstFile1,GetAllSurfSPHARMFiles(0));
		BMSShapeAnalysisModuleMRML<<"  echo()"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  echo(\t'Writing 1st MRML scripts...')"<<std::endl;
		vector<double>Dims;
		SetImageDimensions(firstFile1);
		Dims=GetImageDimensions();
		BMSShapeAnalysisModuleMRML<<" ListFileInDir(MRMLfiles ${FileDir} *SPHARM.vtk)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" GetListSize( nbshape MRMLfiles)"<<std::endl;
	//	BMSShapeAnalysisModuleMRML<<"  Set( MRMLfiles "<< GetListFiles() <<" )"<<std::endl; 
		BMSShapeAnalysisModuleMRML<<"  Set(MRMLScene "<<GetOutputDirectory()<<"/MRML/ShapeAnalysisModuleMRMLScene.mrml)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim0 "<<Dims[0]<<")"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim1 "<<Dims[1]<<")"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim2 "<<Dims[2]<<")"<<std::endl;

	}
	//create a MRML files with all the SPHARM ellalign files
	else if(count==1)
	{
		char *firstFile2;
		firstFile2 = new char[512];
		std::strcpy(firstFile2,GetAllSurfSPHARMellalignFiles(0));
		BMSShapeAnalysisModuleMRML<<"  echo()"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  echo(\t'Writing 2nd MRML scripts...')"<<std::endl;
		vector<double>Dims;
		SetImageDimensions(firstFile2);
		Dims=GetImageDimensions();
		BMSShapeAnalysisModuleMRML<<" ListFileInDir(MRMLfiles ${FileDir} *SPHARM_ellalign.vtk)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" GetListSize( nbshape MRMLfiles)"<<std::endl;
	//	BMSShapeAnalysisModuleMRML<<"  Set( MRMLfiles "<< GetListFiles_ellalign() <<" )"<<std::endl; 
		BMSShapeAnalysisModuleMRML<<"  Set( MRMLScene "<<GetOutputDirectory()<<"/MRML/ShapeAnalysisModuleMRMLScene_ellalign.mrml)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim0 "<<Dims[0]<<")"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim1 "<<Dims[1]<<")"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim2 "<<Dims[2]<<")"<<std::endl;
	}
	//create a MRML file with all the procalign files
	else if(count==2)
	{
		char *firstFile3;
		firstFile3 = new char[512];
		std::strcpy(firstFile3,GetAllSurfSPHARMprocalignFiles(0));
		//std::strcat(firstFile3,"_pp_surfSPHARM_procalign.vtk");
		BMSShapeAnalysisModuleMRML<<"  echo()"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  echo(\t'Writing 3rd MRML scripts...')"<<std::endl;
		vector<double>Dims;
		SetImageDimensions(firstFile3);
		Dims=GetImageDimensions();
		BMSShapeAnalysisModuleMRML<<" ListFileInDir(MRMLfiles ${FileDir} *SPHARM_procalign.vtk)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" GetListSize( nbshape MRMLfiles)"<<std::endl;
	//	BMSShapeAnalysisModuleMRML<<"  Set( MRMLfiles "<< GetListFiles_procalign() <<" )"<<std::endl; 
		BMSShapeAnalysisModuleMRML<<"  Set( MRMLScene "<<GetOutputDirectory()<<"/MRML/ShapeAnalysisModuleMRMLScene_procalign.mrml)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim0 "<<Dims[0]<<")"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim1 "<<Dims[1]<<")"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(dim2 "<<Dims[2]<<")"<<std::endl;
	}
	BMSShapeAnalysisModuleMRML<<"  WriteFile(${MRMLScene} '<MRML userTags=\"\">\\n')"<<std::endl;

	BMSShapeAnalysisModuleMRML<<"  Set(countsnap 1)"<<std::endl;
		
	BMSShapeAnalysisModuleMRML<<"  Set(countj 0)"<<std::endl; //to know which scne is writting
	BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform 0)"<<endl;

	for(int j=0;j<SnapShotNumber*2+3;j++)
	{

		if(j==1)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SceneSnapshot\\n id=\"vtkMRMLSceneSnapshotNode1\" name=\"Color Map Phi\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\">\\n')"<<std::endl;
		}
		else if(j==2)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SceneSnapshot\\n id=\"vtkMRMLSceneSnapshotNode1\" name=\"Color Map Theta\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\">\\n')"<<std::endl;
		}
		else if(j>2)
		{
			if(j%2==0)
				{
					BMSShapeAnalysisModuleMRML<<"  Set(countnbsnapmorethan23 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SceneSnapshot\\n id=\"vtkMRMLSceneSnapshotNode1\" name=\"Color Map Theta'${countsnap}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\">\\n')"<<std::endl;
					
				}
			else
			{	if(j>3){BMSShapeAnalysisModuleMRML<<"  Inc(${countsnap} 1)"<<std::endl;}
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SceneSnapshot\\n id=\"vtkMRMLSceneSnapshotNode1\" name=\"Color Map Phi'${countsnap}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\">\\n')"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countnbsnapmorethan23 0)"<<std::endl;
	
			}
		}
		
		if(j<3){BMSShapeAnalysisModuleMRML<<"   Set(nbLine 1)"<<endl;}


/*	BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Selection\\n id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"NULL\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>\\n <Interaction\\n id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>\\n <Layout\\n id=\"vtkMRMLLayoutNode1\" name=\"vtkMRMLLayoutNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentViewArrangement=\"2\" guiPanelVisibility=\"1\" bottomPanelVisibility =\"1\" guiPanelLR=\"0\" numberOfCompareViewRows=\"0\" numberOfCompareViewColumns=\"0\" numberOfLightboxRows=\"1\" numberOfLightboxColumns=\"1\" mainPanelSize=\"400\" secondaryPanelSize=\"400\"></Layout>\\n <View\\n id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"true\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"true\" backgroundColor=\"0.70196 0.70196 0.90588\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>\\n <Camera\\n id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" position=\"-726.497 88.984 13.4559\" focalPoint=\"0 0 0\" viewUp=\"0 0 1\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>\\n ')"<<std::endl;
	BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<TGParameters\\n id=\"vtkMRMLChangeTrackerNode1\" name=\"vtkMRMLChangeTrackerNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ROIMin=\"-1 -1 -1\" ROIMax=\"-1 -1 -1\" SegmentThresholdMin=\"-1\" SegmentThresholdMax=\"-1\" Analysis_Intensity_Flag=\"0\" Analysis_Deformable_Flag=\"0\" UseITK=\"1\"></TGParameters>\n <VolumeRenderingSelection\\n id=\"vtkMRMLVolumeRenderingSelectionNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" activeVolumeRenderingID=\"NULL\"></VolumeRenderingSelection>\\n <Slice\\n' id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"387.5 250 1\" dimensions=\"496 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n <SliceComposite\\n id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n <Slice\\n id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n <SliceComposite\\n id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n '<Slice\\n id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"0 0 1 0 -1 0 0 0 0 1 0 0 0 0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n <SliceComposite\\n id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n <Crosshair\\n id=\"vtkMRMLCrosshairNode1\" name=\"vtkMRMLCrosshairNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" crosshairMode=\"NoCrosshair\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" crosshairRAS=\"0 0 0\"></Crosshair>\\n <ClipModels\\n id=\"vtkMRMLClipModelsNode1\" name=\"vtkMRMLClipModelsNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" clipType=\"0\" redSliceClipState=\"0\" yellowSliceClipState=\"0\" greenSliceClipState=\"0\"></ClipModels>\\n <ScriptedModule\\n id=\"vtkMRMLScriptedModuleNode1\" name=\"vtkMRMLScriptedModuleNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\"></ScriptedModule>\\n')"<<std::endl;*/
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Selection\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSelectionNode1\" name=\"vtkMRMLSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" secondaryVolumeID=\"NULL\" activeLabelVolumeID=\"NULL\" activeFiducialListID=\"NULL\" activeROIListID=\"NULL\" activeCameraID=\"NULL\" activeViewID=\"NULL\" activeLayoutID=\"vtkMRMLLayoutNode1\"></Selection>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Interaction\\n' )"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLInteractionNode1\" name=\"vtkMRMLInteractionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentInteractionMode=\"ViewTransform\" lastInteractionMode=\"ViewTransform\"></Interaction>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Layout\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLLayoutNode1\" name=\"vtkMRMLLayoutNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" currentViewArrangement=\"2\" guiPanelVisibility=\"1\" bottomPanelVisibility =\"1\" guiPanelLR=\"0\" numberOfCompareViewRows=\"0\" numberOfCompareViewColumns=\"0\" numberOfLightboxRows=\"1\" numberOfLightboxColumns=\"1\" mainPanelSize=\"400\" secondaryPanelSize=\"400\"></Layout>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<View\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLViewNode1\" name=\"vtkMRMLViewNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" active=\"false\" fieldOfView=\"200\" letterSize=\"0.05\" boxVisible=\"true\" fiducialsVisible=\"true\" fiducialLabelsVisible=\"true\" axisLabelsVisible=\"true\" backgroundColor=\"0.70196 0.70196 0.90588\" animationMode=\"Off\" viewAxisMode=\"LookFrom\" spinDegrees=\"2\" spinMs=\"5\" spinDirection=\"YawLeft\" rotateDegrees=\"5\" rockLength=\"200\" rockCount=\"0\" stereoType=\"NoStereo\" renderMode=\"Perspective\"></View>\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Camera\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLCameraNode1\" name=\"vtkMRMLCameraNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" position=\"-726.497 88.984 13.4559\" focalPoint=\"0 0 0\" viewUp=\"0 0 1\" parallelProjection=\"false\" parallelScale=\"1\" active=\"false\"></Camera>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<TGParameters\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLChangeTrackerNode1\" name=\"vtkMRMLChangeTrackerNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ROIMin=\"-1 -1 -1\" ROIMax=\"-1 -1 -1\" SegmentThresholdMin=\"-1\" SegmentThresholdMax=\"-1\" Analysis_Intensity_Flag=\"0\" Analysis_Deformable_Flag=\"0\" UseITK=\"1\"></TGParameters>\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<VolumeRenderingSelection\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLVolumeRenderingSelectionNode1\" name=\"vtkMRMLVolumeRenderingSelectionNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" activeVolumeID=\"NULL\" activeVolumeRenderingID=\"NULL\"></VolumeRenderingSelection>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Slice\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSliceNode1\" name=\"Green\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"387.5 250 1\" dimensions=\"496 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 0 1 0 0 1 0 0 0 0 0 1\" layoutName=\"Green\" orientation=\"Coronal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SliceComposite\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSliceCompositeNode1\" name=\"vtkMRMLSliceCompositeNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Green\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Slice\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSliceNode2\" name=\"Red\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"-1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1\" layoutName=\"Red\" orientation=\"Axial\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SliceComposite\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSliceCompositeNode2\" name=\"vtkMRMLSliceCompositeNode2\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Red\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Slice\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSliceNode3\" name=\"Yellow\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fieldOfView=\"386.719 250 1\" dimensions=\"495 320 1\" activeSlice=\"0\" layoutGridRows=\"1\" layoutGridColumns=\"1\" sliceToRAS=\"0 0 1 0 -1 0 0 0 0 1 0 0 0 0 0 1\" layoutName=\"Yellow\" orientation=\"Sagittal\" jumpMode=\"1\" sliceVisibility=\"false\" widgetVisibility=\"false\" useLabelOutline=\"false\" sliceSpacingMode=\"0\" prescribedSliceSpacing=\"1 1 1\"></Slice>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<SliceComposite\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLSliceCompositeNode3\" name=\"vtkMRMLSliceCompositeNode3\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" backgroundVolumeID=\"\" foregroundVolumeID=\"\" labelVolumeID=\"\" compositing=\"0\" labelOpacity=\"1\" linkedControl=\"0\" foregroundGrid=\"0\" backgroundGrid=\"0\" labelGrid=\"1\" fiducialVisibility=\"1\" fiducialLabelVisibility=\"1\" sliceIntersectionVisibility=\"0\" layoutName=\"Yellow\" annotationMode=\"All\" doPropagateVolumeSelection=\"1\"></SliceComposite>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Crosshair\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLCrosshairNode1\" name=\"vtkMRMLCrosshairNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" crosshairMode=\"NoCrosshair\" crosshairBehavior=\"Normal\" crosshairThickness=\"Fine\" crosshairRAS=\"0 0 0\"></Crosshair>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ClipModels\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLClipModelsNode1\" name=\"vtkMRMLClipModelsNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" clipType=\"0\" redSliceClipState=\"0\" yellowSliceClipState=\"0\" greenSliceClipState=\"0\"></ClipModels>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ScriptedModule\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLScriptedModuleNode1\" name=\"vtkMRMLScriptedModuleNode1\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" ModuleName =\"Editor\" parameter0= \"label 1\"></ScriptedModule>\\n')"<<std::endl;
	
		BMSShapeAnalysisModuleMRML<<"  Set(count 1)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(res 1)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(resPerSnap 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(Dim1 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(count3 1)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(countdata 0)"<<std::endl;
		


		BMSShapeAnalysisModuleMRML<<"  Set(countY 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(countZ 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(countX 0)"<<std::endl;
	
	
	
		
		BMSShapeAnalysisModuleMRML<<"  Set(column 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(columndisplay 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(linedisplay 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Set(nbLine 0)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"   Set(fiducialname 0)"<<endl;
		BMSShapeAnalysisModuleMRML<<"  ForEach(files ${MRMLfiles})"<<std::endl;

		if(j>2){
		//template displayed on every snapshot
			BMSShapeAnalysisModuleMRML<<" set(firstdata 0)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" If(${firstdata} == ${countdata})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"set(tdir '"<<GetOutputDirectory()<<"/Template/')"<<std::endl;
			if(count == 0){BMSShapeAnalysisModuleMRML<<"  listFileInDir(basename ${tdir} *_pp_surfSPHARM.vtk)"<<std::endl;}
			if(count == 1){BMSShapeAnalysisModuleMRML<<"  listFileInDir(basename ${tdir} *_pp_surfSPHARM_ellalign.vtk)"<<std::endl;}
			if(count == 2){BMSShapeAnalysisModuleMRML<<"  listFileInDir(basename ${tdir} *_pp_surfSPHARM_procalign.vtk)"<<std::endl;}
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelStorage\\n id=\"vtkMRMLModelStorageNode'${count}'\" name=\"vtkMRMLModelStorageNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<GetOutputDirectory()<<"/Template/'${data}'\" useCompression=\"1\" readState=\"0\" writeState=\"0\"></ModelStorage>\\n')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Randomize(Rvalue uniform 0 1) "<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Randomize(Gvalue uniform 0 1) "<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Randomize(Bvalue uniform 0 1) "<<std::endl;
	
			if(j%2==0) //if even: want to write a Theta snapshot
			{
				
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"Color Map Theta\" ></ModelDisplay>\\n')"<<std::endl;
				
			}
			else//if odd: want to write a Phi snapshot
			{
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"Color Map Phi\" ></ModelDisplay>\\n')"<<std::endl;
			}
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Model\\n id=\"vtkMRMLModelNode'${count}'\" name=\"'${tdir}''${basename}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\"  storageNodeRef=\"vtkMRMLModelStorageNode'${count}'\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode'${count}'\"></Model>\\n')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Inc(${count} 1)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" EndIf(${firstdata})"<<std::endl;
	
			BMSShapeAnalysisModuleMRML<<"  Inc(${countdatapersnap} 1)"<<std::endl;
		}
	
		BMSShapeAnalysisModuleMRML<<"  GetParam(data ${files} '${count2}' )"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelStorage\\n id=\"vtkMRMLModelStorageNode'${count}'\" name=\"vtkMRMLModelStorageNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" fileName=\""<<GetOutputDirectory()<<"/Mesh/SPHARM/'${data}'\" useCompression=\"1\" readState=\"0\" writeState=\"0\"></ModelStorage>\\n')"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Randomize(Rvalue uniform 0 1) "<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Randomize(Gvalue uniform 0 1) "<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Randomize(Bvalue uniform 0 1) "<<std::endl;
	
		if(j==0)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"false\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"\" ></ModelDisplay>\\n')"<<std::endl;
		}
		else if(j==1)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"Color Map Phi\" ></ModelDisplay>\\n')"<<std::endl;
		}
		else if(j==2)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"Color Map Theta\" ></ModelDisplay>\\n')"<<std::endl;
		}

		else if(j>2)
		{
			BMSShapeAnalysisModuleMRML<<" Math( infj ${countSnapNb} * 2 )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" Math( infj ${infj} + 3 )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" Math( supj ${countSnapNb} + 1)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" Math( supj ${supj} * 2)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" Math( supj ${supj} + 2)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" Math( infsnap ${countSnapNb} * 24 )"<<std::endl;
	
			BMSShapeAnalysisModuleMRML<<" Math( supsnap ${countSnapNb} + 1)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" Math( supsnap ${supsnap} * 24)"<<std::endl;

			if(j>4){
				BMSShapeAnalysisModuleMRML<<" Math( infsnap ${infsnap} + ${countsnap} )"<<std::endl;BMSShapeAnalysisModuleMRML<<" Math( supsnap ${supsnap} + ${countsnap})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<" Math( infsnap ${infsnap} - 1 )"<<std::endl;BMSShapeAnalysisModuleMRML<<" Math( supsnap ${supsnap} - 1)"<<std::endl;
			}	
	
			BMSShapeAnalysisModuleMRML<<" If(${infj} <= ${countj})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" If(${supj} >= ${countj})"<<std::endl;
	
			BMSShapeAnalysisModuleMRML<<" If(${infsnap} <= ${countdata})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" If(${supsnap} >= ${countdata})"<<std::endl;


			if(j%2==0) //if even: want to write a Theta snapshot
			{
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"Color Map Theta\" ></ModelDisplay>\\n')"<<std::endl;
				
			}
			else//if odd: want to write a Phi snapshot
			{
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<ModelDisplay\\n id=\"vtkMRMLModelDisplayNode'${count}'\" name=\"vtkMRMLModelDisplayNode'${count}'\" hideFromEditors=\"true\" selectable=\"true\" selected=\"false\" color=\"'${Rvalue}' '${Gvalue}' '${Bvalue}'\" selectedColor=\"1 0 0\" selectedAmbient=\"0.4\" ambient=\"0.18\" diffuse=\"0.82\" selectedSpecular=\"0.5\" specular=\"0.16\" power=\"1\" opacity=\"1\" visibility=\"true\" clipping=\"false\" sliceIntersectionVisibility=\"false\" backfaceCulling=\"true\" scalarVisibility=\"true\" vectorVisibility=\"false\" tensorVisibility=\"false\" autoScalarRange=\"true\" scalarRange=\"0 100\" colorNodeRef=\"vtkMRMLColorTableNodeFullRainbow\" activeScalarName=\"Color Map Phi\" ></ModelDisplay>\\n')"<<std::endl;
			}

			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Model\\n id=\"vtkMRMLModelNode'${count}'\" name=\""<<GetOutputDirectory()<<"/Mesh/SPHARM/'${data}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode'${count}'\" storageNodeRef=\"vtkMRMLModelStorageNode'${count}'\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode'${count}'\"></Model>\\n')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" EndIf(${supsnap})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" EndIf(${infsnap})"<<std::endl;
	

			if(j>3) //j=3 coordinate save in the parameters SnapCoordY and SnapCoordZ
			//j>3 just need to access to this parameters to write the <LinearTransform
			{
				BMSShapeAnalysisModuleMRML<<" set(infsnaptransform infsnap)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<" set(supsnaptransform supsnap)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<" Math( supsnaptransform ${supsnaptransform} + 1)"<<std::endl;
				
				if(j>4){BMSShapeAnalysisModuleMRML<<" Math( infsnaptransform ${infsnaptransform} + 1)"<<std::endl;}
			
				BMSShapeAnalysisModuleMRML<<" If(${infsnaptransform} <= ${res})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<" If(${supsnaptransform} >= ${res})"<<std::endl;

				if(j==4){BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	}
				if(j>4){BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 1)"<<std::endl;	}

				BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;

				BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
			
				if(j==4){BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;	}
				if(j>4){BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 1)"<<std::endl;	}
			
				BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
							
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;

				//new line
				if(GetDirectionToDisplay()=="ZYX")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(countZ  0)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;
				}
				if(GetDirectionToDisplay()=="YZX")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(countY  0)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				}
				if(GetDirectionToDisplay()=="YXZ")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(countY  0)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;
				}
				if(GetDirectionToDisplay()=="ZXY")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(countZ  0)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;
				}
				if(GetDirectionToDisplay()=="XZY")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(countX  0)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;
				}
				if(GetDirectionToDisplay()=="XYZ")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(countX  0)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;
				}

				//file the<LinearTransform
				if(GetDirectionToDisplay()=="XYZ")
				{

					BMSShapeAnalysisModuleMRML<<" If(${countX} == 4)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countX4_1 ${SnapCoordX} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countX4_0 ${SnapCoordX} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX4 ${countX4_0} - {countX4_1} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(X ${X} + {countX4} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countX} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${X}'  0 1 0 '${Y}'  0 0 1 0  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				}
	
				if(GetDirectionToDisplay()=="XZY")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countX} == 4)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countX4_1 ${SnapCoordX} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countX4_0 ${SnapCoordX} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countX4 ${countX4_0} - {countX4_1} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(X ${X} + {countX4} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countX} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${X}'  0 1 0 0  0 0 1 '${Z}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				}
	
				if(GetDirectionToDisplay()=="YXZ")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countY} == 4)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countY4_1 ${SnapCoordY} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countY4_0 ${SnapCoordY} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY4 ${countY4_0} - {countY4_1} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Y ${Y} + {countY4} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countX} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${X}'  0 1 0 '${Y}'  0 0 1 0  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				}
	
				if(GetDirectionToDisplay()=="YZX")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countY} == 4)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countY4_1 ${SnapCoordY} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countY4_0 ${SnapCoordY} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countY4 ${countY4_0} - {countY4_1} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Y ${Y} + {countY4} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 0  0 1 0 '${Y}'  0 0 1 '${Z}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;

 				}
	
				if(GetDirectionToDisplay()=="ZXY")
				{
					BMSShapeAnalysisModuleMRML<<" If(${countZ} == 4)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countZ4_1 ${SnapCoordZ} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countZ4_0 ${SnapCoordZ} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ4 ${countZ4_0} - {countZ4_1} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Z ${Z} + {countZ4} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countX} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${X}'  0 1 0 0  0 0 1 '${Z}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				}
			
				if(GetDirectionToDisplay()=="ZYX")
				{

					BMSShapeAnalysisModuleMRML<<" If(${countZ} == 4)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countZ4_1 ${SnapCoordZ} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(countZ4_0 ${SnapCoordZ} 1 )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(countZ4 ${countZ4_0} - {countZ4_1} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Z ${Z} + {countZ4} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ} )"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 0  0 1 0 '${Y}'  0 0 1 '${Z}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;

					BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				}

				BMSShapeAnalysisModuleMRML<<" EndIf(${supsnaptransform})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<" EndIf(${infsnaptransform})"<<std::endl;
			}//end if(j>3)
	
			BMSShapeAnalysisModuleMRML<<" EndIf(${supj})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" EndIf(${infj})"<<std::endl;
			
		}//end else if(j>2)

		if(j<3){
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<Model\\n id=\"vtkMRMLModelNode'${count}'\" name=\""<<GetOutputDirectory()<<"/Mesh/SPHARM/'${data}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" transformNodeRef=\"vtkMRMLLinearTransformNode'${count}'\" storageNodeRef=\"vtkMRMLModelStorageNode'${count}'\" userTags=\"\" displayNodeRef=\"vtkMRMLModelDisplayNode'${count}'\"></Model>\\n')"<<std::endl;
		}//end if(j<3)




		if(j<4){//display the scene and the 3 first snapshots
		//BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n')"<<std::endl;
	
			if(j==3){BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;//if 1shape of the snapshot
				BMSShapeAnalysisModuleMRML<<"  Set(res 01)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(count2 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(nbLine 0)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(count3 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;}		

			if(GetDirectionToDisplay()=="XYZ")
			{
				if(j==3){
	
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count2} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 '${Dim1}'  0 0 1 0  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
	
					BMSShapeAnalysisModuleMRML<<"   Set(max 26)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform ${count})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${max} > ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${SnapCoordX} ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${SnapCoordY} ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${max} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
	
			
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 5)"<<std::endl; //only 5 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
			
					
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					//BMSShapeAnalysisModuleMRML<<"  Inc(${nbLine} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
				}
	
				else
				{
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count2} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 '${Dim1}'  0 0 1 0  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 10)"<<std::endl; //only 10 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 10)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					if(j==0)
					{
						BMSShapeAnalysisModuleMRML<<"  If( ${count} == 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"  If( ${count} > 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY  ${fiducialCoordY} ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${fiducialCoordX} ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
					}
				}
			}//end if(GetDirectionToDisplay()=="XYZ")


			else if(GetDirectionToDisplay()=="XZY")
			{	
				if(j==3){
	
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count2} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 0  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
	
					BMSShapeAnalysisModuleMRML<<"   Set(max 26)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform ${count})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${max} > ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${SnapCoordX} ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${SnapCoordZ} ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${max} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
	
			
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 5)"<<std::endl; //only 5 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
			
					
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					}
	
				else
				{
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count2} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 0  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 10)"<<std::endl; //only 10 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 10)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					if(j==0)
					{
						BMSShapeAnalysisModuleMRML<<"  If( ${count} == 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"  If( ${count} > 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ  ${fiducialCoordZ} ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${fiducialCoordX} ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
					}
				}
			}//end else if(GetDirectionToDisplay()=="XZY")

			else if(GetDirectionToDisplay()=="YXZ")
			{	

				if(j==3){
	
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count2})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 '${Dim1}'  0 0 1 0  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
	
					BMSShapeAnalysisModuleMRML<<"   Set(max 26)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform ${count})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${max} > ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${SnapCoordX} ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${SnapCoordY} ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${max} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
	
			
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 5)"<<std::endl; //only 5 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
			
					
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
				}
	
				else
				{
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count2})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 '${Dim1}'  0 0 1 0  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 10)"<<std::endl; //only 10 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 10)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					if(j==0)
					{
						BMSShapeAnalysisModuleMRML<<"  If( ${count} == 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"  If( ${count} > 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY  ${fiducialCoordY} ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${fiducialCoordX} ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
					}
				}

			}//end else if(GetDirectionToDisplay()=="YXZ")


			else if(GetDirectionToDisplay()=="ZXY")
			{	
				if(j==3){
	
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count2} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 0  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
	
					BMSShapeAnalysisModuleMRML<<"   Set(max 26)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform ${count})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${max} > ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordX ${SnapCoordX} ${Dim0})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${SnapCoordZ} ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${max} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
	
			
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 5)"<<std::endl; //only 5 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
			
					
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
				}
		
				else
				{
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count2} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 '${Dim0}'  0 1 0 0  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 10)"<<std::endl; //only 10 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 10)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim0 ${dim0} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					if(j==0)
					{
						BMSShapeAnalysisModuleMRML<<"  If( ${count} == 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"  If( ${count} > 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ  ${fiducialCoordZ} ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim0} + ${dim0})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordX ${fiducialCoordX} ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
					}
				}
			}//end else if(GetDirectionToDisplay()=="ZXY")

			else if(GetDirectionToDisplay()=="ZYX")
			{	
	
				if(j==3){
	
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count2})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 0  0 1 0 '${Dim1}'  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
	
					BMSShapeAnalysisModuleMRML<<"   Set(max 26)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform ${count})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${max} > ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${SnapCoordY} ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${SnapCoordZ} ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${max} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
	
			
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 5)"<<std::endl; //only 5 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
			
					
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
				//	BMSShapeAnalysisModuleMRML<<"  Inc(${nbLine} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
				}
	
				else
				{
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count2})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 0  0 1 0 '${Dim1}'  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 10)"<<std::endl; //only 10 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 10)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					if(j==0)
					{
						BMSShapeAnalysisModuleMRML<<"  If( ${count} == 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"  If( ${count} > 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ  ${fiducialCoordZ} ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY ${fiducialCoordY} ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
					}
				}
			}//end else if(GetDirectionToDisplay()=="ZYX")

			if(GetDirectionToDisplay()=="YZX") 
			{	
				if(j==3){//j==3 save the coordinate in SnapCoordY SnapCoordZ
	
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;//if 1shape of the snapshot
					BMSShapeAnalysisModuleMRML<<"  Set(res 01)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(nbLine 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count3 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count2})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 0  0 1 0 '${Dim1}'  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(max 26)"<<endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} == 0)"<<std::endl;	
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapSaveTransform ${count})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${resPerSnap} > 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${max} > ${resPerSnap} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordZ ${SnapCoordZ} ${Dim2})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(SnapCoordY ${SnapCoordY} ${Dim1})"<<endl;
					BMSShapeAnalysisModuleMRML<<"   Set(fiducialname ${fiducialname} ${data})"<<endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${max} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${resPerSnap} )"<<std::endl;
	
			
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 5)"<<std::endl; //only 5 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
			
					
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 5)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Inc(${nbLine} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
				}
	
				else
				{
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim1 ${dim1} * ${count2})"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<LinearTransform\\n id=\"vtkMRMLLinearTransformNode'${count}'\" name=\"vtkMRMLLinearTransformNode'${count}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" userTags=\"\" matrixTransformToParent=\"1 0 0 0  0 1 0 '${Dim1}'  0 0 1 '${Dim2}'  0 0 0 1\"></LinearTransform>\\n')"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} != 10)"<<std::endl; //only 10 objects per line
					BMSShapeAnalysisModuleMRML<<"  Inc(${count2} 1)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  If( ${res} == 10)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Math(Dim2 ${dim2} * ${count3} )"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  Set(count2 0)"<<std::endl;//start the new line at the 1st position
					BMSShapeAnalysisModuleMRML<<"  Inc(${count3} 1)"<<std::endl;//new line
					BMSShapeAnalysisModuleMRML<<"  Set(res 0)"<<std::endl;
					BMSShapeAnalysisModuleMRML<<"  EndIf( ${res} )"<<std::endl;
					if(j==0)
					{
						BMSShapeAnalysisModuleMRML<<"  If( ${count} == 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"  If( ${count} > 1)"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp2 ${Dim1} + ${dim1})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordY  ${fiducialCoordY} ${temp2})"<<endl;
						BMSShapeAnalysisModuleMRML<<"   Math(temp ${Dim2} + ${dim2})"<<std::endl;
						BMSShapeAnalysisModuleMRML<<"   Set(fiducialCoordZ ${fiducialCoordZ} ${temp})"<<endl;
						BMSShapeAnalysisModuleMRML<<"  EndIf( ${count} )"<<std::endl;
		
					}
				}
	
			}//end if(GetDirectionToDisplay()=="YZX") 

		}//end if(j<4)


		BMSShapeAnalysisModuleMRML<<"  Inc(${count} 1)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" If(${infsnap} <= ${res})"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" If(${supsnap} >= ${res})"<<std::endl;
		
		BMSShapeAnalysisModuleMRML<<"  Inc(${resPerSnap} 1)"<<std::endl;
		
		BMSShapeAnalysisModuleMRML<<" EndIf(${supsnap})"<<std::endl;
		BMSShapeAnalysisModuleMRML<<" EndIf(${infsnap})"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Inc(${res} 1)"<<std::endl;
	
		
		
	
		BMSShapeAnalysisModuleMRML<<"  Inc(${countdata} 1)"<<std::endl;
		BMSShapeAnalysisModuleMRML<<"  Inc(${counter} 1)"<<std::endl;
	
	
	
		BMSShapeAnalysisModuleMRML<<"  EndForEach(files)"<<std::endl;

		if(j>3)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<FiducialList\\n')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Set(countsnapmax countsnap)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"   Math(countsnapmax ${countsnapmax} + 1)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLFiducialListNode'${countsnapmax}'\" name=\"label'${countsnap}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode'${countsnapmax}'\" userTags=\"\" symbolScale=\"0\" symbolType=\"1\" textScale=\"2\" visibility=\"1\" color=\"0 0 0\" selectedcolor=\"0.0117647 0.00784314 0.00784314\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" numberingScheme=\"0\" opacity=\"1\" fiducials=\"')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Set(countfidu 1)"<<std::endl;
			
			if(GetDirectionToDisplay()=="ZYX")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countY 0)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 1)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YZX")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countY 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YXZ")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countY 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countX 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="ZXY")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countX 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XZY")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countX 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XYZ")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countX 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countY 0)"<<std::endl;
			}





			BMSShapeAnalysisModuleMRML<<"ForEach(fiduname ${fiducialname})"<<std::endl;


			BMSShapeAnalysisModuleMRML<<"  If( ${countfidu} < 26)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  RegEx( result ${fiduname} '_pp_.*' REPLACE ' ' )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id '${fiduname}' labeltext '${result}' xyz ')"<<std::endl;

			if(GetDirectionToDisplay()=="ZYX")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 6)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countZ  1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_1 ${SnapCoordZ} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_0 ${SnapCoordZ} 0 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z5 ${Z5_1} - {Z5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z ${Z} + {Z5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;


				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YZX")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countY} == 6)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countY  1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_1 ${SnapCoordY} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_0 ${SnapCoordY} 0 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y5 ${Y5_1} - {Y5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y ${Y} + {Y5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ})"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YXZ")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countY} == 6)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countY  1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_1 ${SnapCoordY} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_0 ${SnapCoordY} 0 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y5 ${Y5_1} - {Y5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y ${Y} + {Y5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countX})"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="ZXY")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 6)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countZ  1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_1 ${SnapCoordZ} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_0 ${SnapCoordZ} 0 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z5 ${Z5_1} - {Z5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z ${Z} + {Z5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XZY")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countX} == 6)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countX  1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_1 ${SnapCoordX} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_0 ${SnapCoordX} 0 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X5 ${X5_1} - {X5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X ${X} + {X5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XYZ")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countX} == 6)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countX  1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_1 ${SnapCoordX} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_0 ${SnapCoordX} 0 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X5 ${X5_1} - {X5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X ${X} + {X5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
			}

	
			BMSShapeAnalysisModuleMRML<<" EndIf(${countfidu})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Inc(${countfidu} 1)"<<std::endl;
	
			BMSShapeAnalysisModuleMRML<<"EndForEach(fiduname)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '\"></FiducialList>\\n')"<<std::endl;
		}//end if(j>3)


		if(j==3)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<FiducialList\\n')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Set(countsnapmax countsnap)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"   Math(countsnapmax ${countsnapmax} + 1)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLFiducialListNode'${countsnapmax}'\" name=\"label'${countsnap}'\" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode'${countsnapmax}'\" userTags=\"\" symbolScale=\"0\" symbolType=\"1\" textScale=\"2\" visibility=\"1\" color=\"0 0 0\" selectedcolor=\"0.0117647 0.00784314 0.00784314\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" numberingScheme=\"0\" opacity=\"1\" fiducials=\"')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Set(countfidu 1)"<<std::endl;

			if(GetDirectionToDisplay()=="ZYX")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countY 0)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 1)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YZX")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countY 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YXZ")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countY 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countX 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="ZXY")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countX 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XZY")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countX 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countZ 0)"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XYZ")
			{
				BMSShapeAnalysisModuleMRML<<"  Set(countX 1)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Set(countY 0)"<<std::endl;
			}


			BMSShapeAnalysisModuleMRML<<"ForEach(fiduname ${fiducialname})"<<std::endl;

			BMSShapeAnalysisModuleMRML<<"  If( ${countfidu} < 6)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  Echo(${countY})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  RegEx( result ${fiduname} '_pp_.*' REPLACE ' ' )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id '${fiduname}' labeltext '${result}' xyz ')"<<std::endl;

			if(GetDirectionToDisplay()=="ZYX")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_1 ${SnapCoordZ} 2 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_0 ${SnapCoordZ} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z5 ${Z5_1} - {Z5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z ${Z} + {Z5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;


				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countZ  -1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YZX")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_1 ${SnapCoordY} 2 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_0 ${SnapCoordY} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y5 ${Y5_1} - {Y5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y ${Y} + {Y5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ})"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countY  -1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="YXZ")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_1 ${SnapCoordY} 2 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y5_0 ${SnapCoordY} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y5 ${Y5_1} - {Y5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Y ${Y} + {Y5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countX})"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countY  -1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="ZXY")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_1 ${SnapCoordZ} 2 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z5_0 ${SnapCoordZ} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z5 ${Z5_1} - {Z5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(Z ${Z} + {Z5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countZ  -1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XZY")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_1 ${SnapCoordX} 2 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_0 ${SnapCoordX} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X5 ${X5_1} - {X5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X ${X} + {X5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countX  -1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
			}

			else if(GetDirectionToDisplay()=="XYZ")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_1 ${SnapCoordX} 2 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X5_0 ${SnapCoordX} 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X5 ${X5_1} - {X5_0} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(X ${X} + {X5} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} < 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countX  -1)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
			}


			BMSShapeAnalysisModuleMRML<<" EndIf(${countfidu})"<<std::endl;

			//.......................................................................................................................................
		
			BMSShapeAnalysisModuleMRML<<"  If( ${countfidu} > 5)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" set(${twentyfive} 25)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" If(${twentyfive} < ${countfidu})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  RegEx( result ${fiduname} '_pp_.*' REPLACE ' ' )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id '${fiduname}' labeltext '${result}' xyz ')"<<std::endl;

			if(GetDirectionToDisplay()=="XYZ")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countX  0)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="XZY")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countX} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countX  0)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countX} )"<<std::endl;
				
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="YXZ")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countY  0)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="YZX")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countY} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countY  0)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countY} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countY} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countZ})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
				
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;

			}
	
			if(GetDirectionToDisplay()=="ZXY")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countX ${countX} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countZ  0)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${SnapCoordX} ${countfidu} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="ZYX")
			{
				BMSShapeAnalysisModuleMRML<<" If(${countZ} == 5)"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} / 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} + 1 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  Math(countY ${countY} * 5 )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"   Set(countZ  0)"<<endl;
				BMSShapeAnalysisModuleMRML<<"  EndIf( ${countZ} )"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${SnapCoordZ} ${countfidu})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${SnapCoordY} ${countfidu} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;

				BMSShapeAnalysisModuleMRML<<"  Math(countZ ${countZ} + 1 )"<<std::endl;
			}

			BMSShapeAnalysisModuleMRML<<" EndIf(${twentyfive})"<<std::endl;
			BMSShapeAnalysisModuleMRML<<" EndIf(${countfidu})"<<std::endl;

			BMSShapeAnalysisModuleMRML<<"  Inc(${countfidu} 1)"<<std::endl;

			BMSShapeAnalysisModuleMRML<<"EndForEach(fiduname)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '\"></FiducialList>\\n')"<<std::endl;
		}//end if(j==3)



		if(j<3)
		{
			BMSShapeAnalysisModuleMRML<<"  # Creating the fiducials )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '<FiducialList\\n')"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} 'id=\"vtkMRMLFiducialListNode1\" name=\"label: \" hideFromEditors=\"false\" selectable=\"true\" selected=\"false\" storageNodeRef=\"vtkMRMLFiducialListStorageNode1\" userTags=\"\" symbolScale=\"0\" symbolType=\"1\" textScale=\"2\" visibility=\"1\" color=\"0 0 0\" selectedcolor=\"0.0117647 0.00784314 0.00784314\" ambient=\"0\" diffuse=\"1\" specular=\"0\" power=\"1\" locked=\"0\" opacity=\"1\" fiducials=\"')"<<std::endl;
	
			BMSShapeAnalysisModuleMRML<<"  Set(counter 0)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  ForEach(files ${MRMLfiles})"<<std::endl;
		
			BMSShapeAnalysisModuleMRML<<"  GetParam(data ${files} '${counter}' )"<<std::endl;

			BMSShapeAnalysisModuleMRML<<"  RegEx( result ${data} '_pp_.*' REPLACE ' ' )"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ' id '${result}' labeltext '${result}' xyz ')"<<std::endl;
		
		
			if(GetDirectionToDisplay()=="XYZ")
			{
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${fiducialCoordX} ${counter})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${fiducialCoordY} ${counter} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="XZY")
			{
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${fiducialCoordX} ${counter})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${fiducialCoordZ} ${counter} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="YXZ")
			{
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${fiducialCoordX} ${counter})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${fiducialCoordY} ${counter} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' '${Y}' 15 orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="YZX")
			{
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${fiducialCoordZ} ${counter})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${fiducialCoordY} ${counter} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="ZXY")
			{
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${fiducialCoordZ} ${counter})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(X ${fiducialCoordX} ${counter} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} ''${X}' 15 '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
			}
		
			if(GetDirectionToDisplay()=="ZYX")
			{
				BMSShapeAnalysisModuleMRML<<"  GetParam(Z ${fiducialCoordZ} ${counter})"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  GetParam(Y ${fiducialCoordY} ${counter} )"<<std::endl;
				BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '15 '${Y}' '${Z}' orientationwxyz 0 0 0 1 selected 1 visibility 1 ')"<<std::endl;
			}
		
		
			BMSShapeAnalysisModuleMRML<<"  Inc(${counter} 1)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  EndForEach(files)"<<std::endl;
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '\"></FiducialList>\\n')"<<std::endl;
		}//end if(j<3)

	
		BMSShapeAnalysisModuleMRML<<"  Inc(${countj} 1)"<<std::endl;

		if(j>2){
		if(j%2==0){BMSShapeAnalysisModuleMRML<<"Inc(${countSnapNb} 1)"<<std::endl;}
			}

		if(j!=0)
		{
			BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '</SceneSnapshot>\\n')"<<std::endl;
		}
	}
	BMSShapeAnalysisModuleMRML<<"  AppendFile(${MRMLScene} '</MRML>\\n')"<<std::endl;
	}
	
	
	BMSShapeAnalysisModuleMRML<<" Set(BMSFile "<<GetOutputDirectory()<<"/BatchMake_Scripts/ShapeAnalysisModuleMRML.bms)"<<std::endl;
	BMSShapeAnalysisModuleMRML<<" CopyFile("<<GetBMSShapeAnalysisModuleMRMLFile()<<" ${BMSFile})"<<std::endl;
	//BMSShapeAnalysisModuleMRML<<" DeleteFile("<<GetBMSShapeAnalysisModuleMRMLFile()<<")"<<std::endl;

	BMSShapeAnalysisModuleMRML<<"  echo(\t'Writing MRML script: Done!')"<<std::endl;
	BMSShapeAnalysisModuleMRML<<"  echo()"<<std::endl;

	return;
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
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Mesh/PostProcess/)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Mesh/SPHARM/)"<<std::endl;
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Template/)"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/Output/)"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"MakeDirectory("<<GetOutputDirectory()<<"/EulerFiles/)"<<std::endl; 
	BMSShapeAnalysisModuleFile<<"set(BMSdir '"<<GetOutputDirectory()<<"/BatchMake_Scripts')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(datadir '"<<GetOutputDirectory()<<"/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(PPdir '"<<GetOutputDirectory()<<"/Mesh/PostProcess/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(SPHARMdir '"<<GetOutputDirectory()<<"/Mesh/SPHARM/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(tdir '"<<GetOutputDirectory()<<"/Template/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(Outputdir '"<<GetOutputDirectory()<<"/Output/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"set(Eulerdir '"<<GetOutputDirectory()<<"/EulerFiles/')"<<std::endl;
	BMSShapeAnalysisModuleFile<<"echo()"<<std::endl;
	
	BMSShapeAnalysisModuleFile<<"#Create OutputFile"<<std::endl;
	BMSShapeAnalysisModuleFile<<"Set(OutputFile "<<GetOutputDirectory()<<"/Output/ShapeAnalysisModule_OutputFileVersion1.csv)"<<std::endl;
	SetOuputFile();
	BMSShapeAnalysisModuleFile<<"Set(OutputFile "<<GetOutputFile()<<")"<<std::endl;

	BMSShapeAnalysisModuleFile<<" WriteFile(${OutputFile} '"<<OutputFileHeaders[0]<<",')"<<std::endl;
	for(unsigned int i=1;i<OutputFileHeaders.size();i++)
	{
		BMSShapeAnalysisModuleFile<<" appendFile(${OutputFile} '"<<OutputFileHeaders[i]<<",')"<<std::endl;
	}
	
	//Write headers of the output file
	BMSShapeAnalysisModuleFile<<" appendFile(${OutputFile} ' Post Processed Segmentation, Parametrisation of Original Surface, SPHARM Surface in Original Space, SPHARM Coefficient in Original Space, SPHARM Surface in Ellipsoid Aligned Space, SPHARM Coefficient in Ellipsoid Aligned Space, SPHARM Surface in Procaligned Space\\n')"<<std::endl;
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

