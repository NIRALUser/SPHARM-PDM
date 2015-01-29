#include "ShapeAnalysisModuleComputation.h"
#include <itksys/SystemTools.hxx>
#include <itksys/Directory.hxx>

ShapeAnalysisModuleComputation::ShapeAnalysisModuleComputation()
  : Parameters()
{
}

ShapeAnalysisModuleComputation::~ShapeAnalysisModuleComputation()
{
}

// Compute Shape Analysis
int ShapeAnalysisModuleComputation::Computation()
{

  std::cout << "\n\nComputing ShapeAnalysisModule..." << std::endl << std::endl;

  SetBMSShapeAnalysisModuleFile(false);
  SetBMSShapeAnalysisModuleMRMLFile(false);

  m_nbHorizontal = GetHorizontalGridPara();
  m_nbVertical = GetVerticalGridPara();
  this->m_nbShapesPerMRML = m_nbHorizontal * m_nbVertical;

  SetAllFilesName();
  OverWrite();
  WriteBMSShapeAnalysisModuleFile();
  ExecuteBatchMake(GetBMSShapeAnalysisModuleFile() ); 

  if( GetTemplateMState() == true )
    {
    ComputationMean();

    WriteBMSShapeAnalysisModuleFile2();
    ExecuteBatchMake(GetBMSShapeAnalysisModuleFile2() );
    }

  std::cout << "\n\nExecute Meshmath..." << std::endl << std::endl;
  // execute MeshMath external application
  for( int i = 0; i < GetDataNumber(); i++ )
    {
    ExecuteMeshMath(i, "phi", 0);
    ExecuteMeshMath(i, "theta", 0);
    }
#ifdef WIN32 
  //ExecuteMeshMathTemplate();
#else
  ExecuteMeshMathTemplate();
#endif

  // Particles
  if( GetParticlesState() )
    {
    RunParticlesModule();
    std::cout << "Modify output csv" << std::endl;
    ModifyCSV(1);
    for( int i = 0; i < GetDataNumber(); i++ )
      {
      ExecuteMeshMath(i, "phi", 1);
      ExecuteMeshMath(i, "theta", 1);
      }

    }
  else
    {
    ModifyCSV(0);
    }

  std::cout << "DONE COMPUTING SPHARM" << std::endl;

  return 0;
}

// Execute MeshMath to write a KWM scalar field (1D) into a PolyData Field Data Scalar to visualize in Slicer3
void ShapeAnalysisModuleComputation::ExecuteMeshMath(int numData, const char * scalar, bool particule)
{
  const std::string _scalar(scalar);
  int end;

  if( particule == 0 )
    {
    end = 3;
    }                       // execute MeshMath for each volume file: Original, Ellalign, Procalign
  else
    {
    end = 1;
    }

  for( int j = 0; j < end; j++ )
    {

    std::vector<const char *> args;
	std::string cmd ;
    char*                     data = NULL;
    int                       length;
    double                    timeout = 0.05;
    int                       result;
    char *                    fileType = NULL;


	args.push_back("MeshMath");

    if( particule == 0 )
      {
      if( GetTemplateMState() == true )
        {
        if( j == 0 )
          {
            fileType = GetAllSurfmeanSPHARMFiles(numData);
	    fileType = GetAllSurfSPHARMFiles(numData);
          }
        else if( j == 1 )
          {
            fileType = GetAllSurfmeanSPHARMellalignFiles(numData);
	    fileType = GetAllSurfSPHARMellalignFiles(numData);
          }
        else if( j == 2 )
          {
            fileType = GetAllSurfmeanSPHARMprocalignFiles(numData);
	    fileType = GetAllSurfSPHARMprocalignFiles(numData);
          }
        }
      else
        {
        if( j == 0 )
          {
            fileType = GetAllSurfSPHARMFiles(numData);
          }
        else if( j == 1 )
          {
            fileType = GetAllSurfSPHARMellalignFiles(numData);
          }
        else if( j == 2 )
          {
            fileType = GetAllSurfSPHARMprocalignFiles(numData);
          }
        }
      }
    else
      {
      std::string tmp;
      if( GetTemplateMState() == true )
        {
        }
      else
        {
        fileType = GetPostCorrespondenceFiles(numData);
        }
      }

    args.push_back(fileType);
    args.push_back(fileType);
    args.push_back("-KWMtoPolyData");
	std::string getAllPhiFilesval ;
	std::string getAllThetaFilesval ;
    if( _scalar == "phi" )
      {
		  getAllPhiFilesval = GetAllPhiFiles(numData) ;
		  args.push_back(getAllPhiFilesval.c_str() );
          args.push_back("Color_Map_Phi");
      }
    if( _scalar == "theta" )
      {
		 getAllThetaFilesval = GetAllThetaFiles(numData) ; 
		 args.push_back(getAllThetaFilesval.c_str() );
         args.push_back("Color_Map_Theta");
      }

    args.push_back(0);

    // Run the application
    itksysProcess* gp = itksysProcess_New();
    itksysProcess_SetCommand(gp, &*args.begin() );
	//itksysProcess_SetCommand(gp, &cmd_char );
    itksysProcess_SetOption(gp, itksysProcess_Option_HideWindow, 1);
    itksysProcess_Execute(gp);

    while( int Value = itksysProcess_WaitForData(gp, &data, &length, &timeout) ) // wait for 1s
      {
      if( ( (Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR) ) && data[0] == 'D' )
        {
        strstream st;
        for( int i = 0; i < length; i++ )
          {
          st << data[i];
          }
        string        dim = st.str();
        istringstream s(dim);
        string        value;

        while( getline(s, value, ' ') )
          {
          m_Dims.push_back( (atoi(value.c_str() ) ) / 2);
          }

        m_Dims.erase(m_Dims.begin() );
        }
      timeout = 0.05;
      }

    itksysProcess_WaitForExit(gp, 0);

    result = 1;

    switch( itksysProcess_GetState(gp) )
      {
      case itksysProcess_State_Exited:
        {
        result = itksysProcess_GetExitValue(gp);
        } break;
      case itksysProcess_State_Error:
        {
        std::cerr << "Error: Could not run " << args[0] << ":\n";
        std::cerr << itksysProcess_GetErrorString(gp) << "\n";
        std::cout << "Error: Could not run " << args[0] << ":\n";
        std::cout << itksysProcess_GetErrorString(gp) << "\n";
        } break;
      case itksysProcess_State_Exception:
        {
        std::cerr << "Error: " << args[0] << " terminated with an exception: "
                  << itksysProcess_GetExceptionString(gp) << "\n";
        std::cout << "Error: " << args[0] << " terminated with an exception: "
                  << itksysProcess_GetExceptionString(gp) << "\n";
        } break;
      case itksysProcess_State_Starting:
      case itksysProcess_State_Executing:
      case itksysProcess_State_Expired:
      case itksysProcess_State_Killed:
        {
        // Should not get here.
        std::cerr << "Unexpected ending state after running " << args[0] << std::endl;
        std::cout << "Unexpected ending state after running " << args[0] << std::endl;
        } break;
      }
    itksysProcess_Delete(gp);

    }
}

void ShapeAnalysisModuleComputation::ExecuteMeshMathTemplate()
{
  std::cout << "Mapping parameters in SPHARM template meshes using MeshMath" << std::endl;

  itksys::Glob             globTemplateList;
  std::vector<std::string> list_template;
  std::string              pathFile = GetOutputDirectory();
  std::string              path;
  //path = "/Template/*.vtk";
  path = "/Template/";
  pathFile = pathFile + path;
  std::cout << "Path for MeshMathTemplate: " << pathFile << std::endl;
  globTemplateList.FindFiles(pathFile);
  list_template = globTemplateList.GetFiles();
  for( int color_map = 0;  color_map < 2; color_map++ )
    {
    for( unsigned int  i = 0; i < list_template.size(); i++ )
      {
      std::vector<const char *> args;
      char*                     data = NULL;
      int                       length;
      double                    timeout = 0.05;
      int                       result;

      args.push_back("MeshMath");

      args.push_back(list_template[i].c_str() );
      args.push_back(list_template[i].c_str() );
      args.push_back("-KWMtoPolyData");

      if( color_map == 0 )
        {
        args.push_back(GetAllPhiFiles(0).c_str() );
        args.push_back("Color_Map_Phi");
        }
      if( color_map == 1 )
        {
        args.push_back(GetAllThetaFiles(0).c_str() );
        args.push_back("Color_Map_Theta");
        }
	  
      args.push_back(0);
	
      // Run the application
      itksysProcess* gp = itksysProcess_New();
      itksysProcess_SetCommand(gp, &*args.begin() );
      itksysProcess_SetOption(gp, itksysProcess_Option_HideWindow, 1);
      itksysProcess_Execute(gp);

      while( int Value = itksysProcess_WaitForData(gp, &data, &length, &timeout) ) // wait for 1s
        {
        if( ( (Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR) ) && data[0] == 'D' )
          {
          strstream st;
          for( int j = 0; j < length; ++j )
            {
            st << data[j];
            }
          string        dim = st.str();
          istringstream s(dim);
          string        value;

          while( getline(s, value, ' ') )
            {
            m_Dims.push_back( (atoi(value.c_str() ) ) / 2);
            }

          m_Dims.erase(m_Dims.begin() );
          }
        timeout = 0.05;
        }

      itksysProcess_WaitForExit(gp, 0);

      result = 1;

      switch( itksysProcess_GetState(gp) )
        {
        case itksysProcess_State_Exited:
          {
          result = itksysProcess_GetExitValue(gp);
          } break;
        case itksysProcess_State_Error:
          {
          std::cerr << "Error: Could not run " << args[0] << ":\n";
          std::cerr << itksysProcess_GetErrorString(gp) << "\n";
          std::cout << "Error: Could not run " << args[0] << ":\n";
          std::cout << itksysProcess_GetErrorString(gp) << "\n";
          } break;
        case itksysProcess_State_Exception:
          {
          std::cerr << "Error: " << args[0] << " terminated with an exception: "
                    << itksysProcess_GetExceptionString(gp) << "\n";
          std::cout << "Error: " << args[0] << " terminated with an exception: "
                    << itksysProcess_GetExceptionString(gp) << "\n";
          } break;
        case itksysProcess_State_Starting:
        case itksysProcess_State_Executing:
        case itksysProcess_State_Expired:
        case itksysProcess_State_Killed:
          {
          // Should not get here.
          std::cerr << "Unexpected ending state after running " << args[0] << std::endl;
          std::cout << "Unexpected ending state after running " << args[0] << std::endl;
          } break;
        }
      itksysProcess_Delete(gp);
      }

    }
}

// Create BMS File to compute the SPHARM pipeline
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleFile(bool changeDirectory)
{
  std::strcpy(m_BMSShapeAnalysisModuleFile, GetOutputDirectory() );
//  std::cout << GetBMSShapeAnalysisModuleFile() << std::endl;

 // if( changeDirectory == false )
  //  {
    std::strcat(m_BMSShapeAnalysisModuleFile, "/");
    std::cout << GetBMSShapeAnalysisModuleFile() << std::endl;
   // }

 /* else
    {
    std::strcat(m_BMSShapeAnalysisModuleFile, "/BatchMake_Scripts/");
    std::cout << GetBMSShapeAnalysisModuleFile() << std::endl;
    }

  if( GetRandomizeInputs() )
    {
    int pID = getpid();
    sprintf(m_BMSShapeAnalysisModuleFile, "ShapeAnalysisModule_id%d.bms", pID);
    std::cout << " " << std::endl; std::cout << GetBMSShapeAnalysisModuleFile() << std::endl;
    }
  else
    {*/

    std::strcat(m_BMSShapeAnalysisModuleFile,
                "ShapeAnalysisModule.bms"); 
   // std::cout << GetBMSShapeAnalysisModuleFile() << std::endl;
  //  }

  std::cout << GetBMSShapeAnalysisModuleFile() << std::endl;
  return;
}

// Get BMS File
char * ShapeAnalysisModuleComputation::GetBMSShapeAnalysisModuleFile()
{
  return m_BMSShapeAnalysisModuleFile;
}

// Create BMS for mean computation File
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleFile2(bool changeDirectory)
{
  std::strcpy(m_BMSShapeAnalysisModuleFile, GetOutputDirectory() );

  if( changeDirectory == false )
    {
    std::strcat(m_BMSShapeAnalysisModuleFile, "/");
    }

  else
    {
    std::strcat(m_BMSShapeAnalysisModuleFile, "/BatchMake_Scripts/");
    }

  std::strcat(m_BMSShapeAnalysisModuleFile, "ShapeAnalysisModule_MeanAsTemplate.bms");

  return;
}

char * Convert_Double_To_CharArray(double doubleVariable)
{
  char *CharArrayOut;

  CharArrayOut = new char[512];
  std::string       stringVariable;
  std::stringstream strStream;
  strStream << doubleVariable;
  stringVariable = strStream.str();
  strcpy(CharArrayOut, stringVariable.c_str() );

  return CharArrayOut;
}

// Get BMS File for mean computation
char * ShapeAnalysisModuleComputation::GetBMSShapeAnalysisModuleFile2()
{
  return m_BMSShapeAnalysisModuleFile;
}

// Create a BMS file to create a MRML scene
void ShapeAnalysisModuleComputation::SetBMSShapeAnalysisModuleMRMLFile(bool /* changeDirectory */)
{
  std::strcpy(m_BMSShapeAnalysisModuleMRLMFile, GetOutputDirectory() );
  std::strcat(m_BMSShapeAnalysisModuleMRLMFile, "/BatchMake_Scripts/");
  std::strcat(m_BMSShapeAnalysisModuleMRLMFile, "ShapeAnalysisModuleMRML.bms");
}

char * ShapeAnalysisModuleComputation::GetBMSShapeAnalysisModuleMRMLFile()
{
  return m_BMSShapeAnalysisModuleMRLMFile;
}

// create the output file
void ShapeAnalysisModuleComputation::SetOuputFile()
{
  std::strcpy(m_OutputFile, GetOutputDirectory() );
  std::strcat(m_OutputFile, "/OutputGroupFile/");
  std::strcat(m_OutputFile, "ShapeAnalysisModule_OutputFileVersion1.csv");
}

char * ShapeAnalysisModuleComputation::GetOutputFile()
{
  return m_OutputFile;
}

// Execute BatchMake script
void ShapeAnalysisModuleComputation::ExecuteBatchMake(const char *_Input)
{
  std::cout << "\tExecuting BatchMake..." << std::endl;

  char *           envpath = getenv("BatchmakeShapeAnalysisModule_Dir");
  //std::cout<<"envpath:"<<envpath<<std::endl;
  std::string      applicationPath;
  bm::ScriptParser m_Parser;

  // the module path exist
  if( itksys::SystemTools::ChangeDirectory(GetModulePath() ) == 0 )
    {
    // if bmm files are in the same directory as the Shape Analysis Module
    m_Parser.LoadWrappedApplication(GetModulePath() );
    }

  else
    {
    // if the environment variable is set
    if( envpath )
      {
      applicationPath = std::string(envpath);
      m_Parser.LoadWrappedApplication(applicationPath.c_str() );
      }

    else
      {
      std::cerr << "The environment variable 'BatchmakeShapeAnalysisModule_Dir' needs to be set" << std::endl;
      std::cerr << "bash usage : export BatchmakeShapeAnalysisModule_Dir=<Batchmake ShapeAnalysisModule Directory>"
                << std::endl;
      std::cerr << "tcsh usage : setenv BatchmakeShapeAnalysisModule_Dir <Batchmake ShapeAnalysisModule Directory>"
                << std::endl;
      exit(0);
      }
    }
  
  bool ComputationSuccess;
  if ( GetRandomizeInputs() )
  {
    std::string batchmake_file;
    batchmake_file.append(GetOutputDirectory());
    batchmake_file.append(_Input);
    ComputationSuccess = m_Parser.Execute(batchmake_file);
  }
  else
    ComputationSuccess = m_Parser.Execute(_Input);
  
  if( !ComputationSuccess )
    {
    cerr << "\tExecuting BatchMake: Error!" << endl;
    }

  std::cout << "\tExecuting BatchMake: Done!" << std::endl << endl;
  return;
}

// Write Batchmake script to compute Shape Analysis
void ShapeAnalysisModuleComputation::WriteBMSShapeAnalysisModuleFile()
{

  std::ofstream BMSShapeAnalysisModuleFile(GetBMSShapeAnalysisModuleFile() );

  vector<string> OutputFileHeaders = GetOutputFileHeaders();
  int            DataNumber;

  BMSShapeAnalysisModuleFile << "#---------------------------- Shape Analysis ----------------------------"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "#------------------------------------------------------------------------"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "#------------------------------------------------------------------------"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "#------------------------------------------------------------------------"
                             << std::endl;
  for( int i = 0; i < GetDataNumber(); i++ )
    {
    m_permutations.push_back(0);
    }
  // get here the permutation vector
  GetRandomNum(1, GetDataNumber() );
  // for( int i = 0; i < GetDataNumber(); i++)
  //	std::cout << m_permutations[i] << std::endl;

  if( GetRandomizeInputs() )
    {
    for( unsigned int i = 0; i < OutputFileHeaders.size(); i++ )
      {
      BMSShapeAnalysisModuleFile << " Set(Header" << i << " '"
                                 << GetNthDataListValue(m_permutations[0], i) << "')" << std::endl;

      }

    BMSShapeAnalysisModuleFile << "set (OrigCasesList '"
                               << GetNthDataListValue(m_permutations[0], GetColumnVolumeFile() ) << "')" << std::endl;
    }
  else
    {
    // Set the headers for the output file
    for( unsigned int i = 0; i < OutputFileHeaders.size(); i++ )
      {
      BMSShapeAnalysisModuleFile << " Set(Header" << i << " '" << GetNthDataListValue(1, i) << "')" << std::endl;

      }

    BMSShapeAnalysisModuleFile << "set (OrigCasesList '"
                               << GetNthDataListValue(1, GetColumnVolumeFile() ) << "')" << std::endl;
    }

  std::cout << "Number of Datas: " << GetDataNumber() << std::endl;

  if( GetRandomizeInputs() )
    {
    for( DataNumber = 1; DataNumber < GetDataNumber(); DataNumber++ )
      {
      for( unsigned int i = 0; i < OutputFileHeaders.size(); i++ )
        {
        BMSShapeAnalysisModuleFile << "set (Header" << i << " ${Header" << i << "} '" << GetNthDataListValue(
          m_permutations[DataNumber], i) << "')" << std::endl;

        }
      BMSShapeAnalysisModuleFile << "set (OrigCasesList ${OrigCasesList} '"
                                 << GetNthDataListValue(m_permutations[DataNumber],
                             GetColumnVolumeFile() ) << "')" << std::endl;
      }
    }
  else
    {
    for( DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
      {
      for( unsigned int i = 0; i < OutputFileHeaders.size(); i++ )
        {
        BMSShapeAnalysisModuleFile << "set (Header" << i << " ${Header" << i << "} '" << GetNthDataListValue(DataNumber,
                                                                                                             i)
                                   << "')" << std::endl;

        }
      BMSShapeAnalysisModuleFile << "set (OrigCasesList ${OrigCasesList} '" << GetNthDataListValue(DataNumber,
                                                                                                   GetColumnVolumeFile() )
                                 << "')" << std::endl;
      }
    }

  BMSShapeAnalysisModuleFile << "#Create directories" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/BatchMake_Scripts/)" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/Mesh/)" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/Mesh/PostProcess/)" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/Mesh/SPHARM/)" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/Template/)" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/OutputGroupFile/)" << std::endl;
  BMSShapeAnalysisModuleFile << "MakeDirectory(" << GetOutputDirectory() << "/EulerFiles/)" << std::endl;
  BMSShapeAnalysisModuleFile << "set(BMSdir '" << GetOutputDirectory() << "/BatchMake_Scripts')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(datadir '" << GetOutputDirectory() << "/')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(PPdir '" << GetOutputDirectory() << "/Mesh/PostProcess/')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(SPHARMdir '" << GetOutputDirectory() << "/Mesh/SPHARM/')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(tdir '" << GetOutputDirectory() << "/Template/')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(Outputdir '" << GetOutputDirectory() << "/OutputGroupFile/')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(Eulerdir '" << GetOutputDirectory() << "/EulerFiles/')" << std::endl;

  BMSShapeAnalysisModuleFile << "echo()" << std::endl;

  BMSShapeAnalysisModuleFile << "#Create OutputFile" << std::endl;
  BMSShapeAnalysisModuleFile << "Set(OutputFile " << GetOutputDirectory()
                             << "/OutputGroupFile/ShapeAnalysisModule_OutputFileVersion1.csv)" << std::endl;
  SetOuputFile();
  BMSShapeAnalysisModuleFile << "Set(OutputFile " << GetOutputFile() << ")" << std::endl;

  BMSShapeAnalysisModuleFile << " WriteFile(${OutputFile} '" << OutputFileHeaders[0] << ",')" << std::endl;
  for( unsigned int i = 1; i < OutputFileHeaders.size(); i++ )
    {
    BMSShapeAnalysisModuleFile << " appendFile(${OutputFile} '" << OutputFileHeaders[i] << ",')" << std::endl;
    }

  // Write headers of the output file
  BMSShapeAnalysisModuleFile
  <<
  " appendFile(${OutputFile} ' Post Processed Segmentation, Parameterization of Original Surface, SPHARM Surface in Original Space, SPHARM Coefficient in Original Space, SPHARM Surface in Ellipsoid Aligned Space, SPHARM Coefficient in Ellipsoid Aligned Space, SPHARM Surface in Procaligned Space\\n')"
  << std::endl;
  BMSShapeAnalysisModuleFile << "echo()" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(count 0)" << std::endl;

  BMSShapeAnalysisModuleFile << "ForEach(case ${OrigCasesList})" << std::endl;
  BMSShapeAnalysisModuleFile << "  #Extract basename" << std::endl;
  BMSShapeAnalysisModuleFile << "  GetFilename(basename ${case} NAME_WITHOUT_EXTENSION)" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo(Case: ${case})" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testSeg ${PPdir} *${basename}*_pp" << GetVolumeFileExtension()
                             << ")" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testGen1 ${SPHARMdir} *${basename}*para.vtk*) " << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testGen2 ${SPHARMdir} *${basename}*surf.vtk*) " << std::endl;
  BMSShapeAnalysisModuleFile << "  set(testGen ${testGen1} ${testGen2}) " << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testPara ${SPHARMdir} *${basename}*SPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testPara2 ${SPHARMdir} *${basename}*MeanSPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testtemp ${tdir}  *SPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testtemp2 ${tdir} *${basename}*SPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  set(ppcase ${PPdir}${basename}_pp.gipl.gz)" << std::endl;

  BMSShapeAnalysisModuleFile << "  #Post Processing" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo('Doing Post Processing')" << std::endl;
  BMSShapeAnalysisModuleFile << "  if(${testSeg} == '')" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetApp(Seg @SegPostProcessCLP)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.fileName ${case})" << std::endl;

  if( GetGaussianFilteringState() == true )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.gaussianOn 1)" << std::endl;
    BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.variance_vect.variance_vect " << GetVarianceBoxX() << ","
                               << GetVarianceBoxY() << "," << GetVarianceBoxZ() << ")" << std::endl;
    }
  BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.outfileName ${ppcase})" << std::endl;
  if( GetLabelState() == true )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.label.label " << GetLabel() << ")" << std::endl;
    }
  if( GetRescaleState() == true )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.scaleOn 1)" << std::endl;
    }

  BMSShapeAnalysisModuleFile << "    SetAppOption(Seg.spacing_vect.spacing_vect " << GetEnforcedSpaceX() << ","
                             << GetEnforcedSpaceY() << "," << GetEnforcedSpaceZ() << ")" << std::endl;

  BMSShapeAnalysisModuleFile << "    Run(output ${Seg} error)" << std::endl;
  BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
  BMSShapeAnalysisModuleFile << "    echo('SegPostProcess Error:' ${error})" << std::endl;
  // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${testSeg})" << std::endl;
  

  BMSShapeAnalysisModuleFile << "  #GenParaMesh" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo('Doing GenParaMesh')" << std::endl;
  BMSShapeAnalysisModuleFile << "  Set(value 0)" << std::endl;
  BMSShapeAnalysisModuleFile << "  foreach(data ${testGen})" << std::endl;
  BMSShapeAnalysisModuleFile << "    Inc(${value} 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    Int(${value})" << std::endl;
  BMSShapeAnalysisModuleFile << "  EndForEach(data)" << std::endl;
  BMSShapeAnalysisModuleFile << "  if (${value} < 2)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetApp(Gen @GenParaMeshCLP)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.infile ${ppcase})" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.outParaName ${SPHARMdir}${basename}_pp_para.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.outSurfName ${SPHARMdir}${basename}_pp_surf.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.numIterations.numIterations " << GetNumberOfIterations()
                             << ")" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.label.label 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.EulerFile 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.outEulerName 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Gen.outEulerName.outEulerName ${Eulerdir}${basename}_euler.txt)"
                             << std::endl;

  BMSShapeAnalysisModuleFile << "    Run(output ${Gen} error)" << std::endl;
  BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
  BMSShapeAnalysisModuleFile << "    echo('GenParaMesh Error:' ${error})" << std::endl;
  // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${value})" << std::endl;

  BMSShapeAnalysisModuleFile << "  #ParaToSPHARMMesh" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo('Doing ParaToSPHARMMesh')" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(value 0)" << std::endl;
  BMSShapeAnalysisModuleFile << "  foreach(data ${testtemp})" << std::endl;
  BMSShapeAnalysisModuleFile << "    Inc(${value} 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    Int(${value})" << std::endl;
  BMSShapeAnalysisModuleFile << "  EndForEach(data)" << std::endl;
  BMSShapeAnalysisModuleFile << "  if (${value} < 4)" << std::endl;

// TODO create templte

  BMSShapeAnalysisModuleFile << "  #Create Template" << std::endl;
  BMSShapeAnalysisModuleFile << "    echo('Creating Template')" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetApp(ParaT @ParaToSPHARMMeshCLP)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.outbase ${tdir}${basename}_pp_surf)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.subdivLevel 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.subdivLevel.subdivLevel " << GetSubdivLevel() << ")"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.spharmDegree 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.spharmDegree.spharmDegree " << GetSPHARMDegree() << ")"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.thetaIteration 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.thetaIteration.thetaIteration " << GetThetaIteration() << ")"
		  << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.phiIteration 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.phiIteration.phiIteration " << GetPhiIteration() << ")"
		  << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh.medialMesh " << GetPhiIteration() << ")"
		  << std::endl;
  BMSShapeAnalysisModuleFile << "      SetAppOption(ParaT.finalFlipIndex 1)" << std::endl;
  if( GetFinalFlipN() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex " << 0 << ")" << std::endl;
    }
  if( GetFinalFlipX() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 4 << ")" << std::endl;
    }
  if( GetFinalFlipY() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 5 << ")" << std::endl;
    }
  if( GetFinalFlipZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 7 << ")" << std::endl;
    }
  if( GetFinalFlipXY() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 1 << ")" << std::endl;
    }
  if( GetFinalFlipYZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 2 << ")" << std::endl;
    }
  if( GetFinalFlipXZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 3 << ")" << std::endl;
    }
  if( GetFinalFlipXYZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.finalFlipIndex.finalFlipIndex  " << 6 << ")" << std::endl;
    }

  BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *SPHARM.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "      listFileInDir(flip ${tdir} *SPHARM.coef)" << std::endl;
  BMSShapeAnalysisModuleFile << "      set(regTemplate ${reg})" << std::endl;
  BMSShapeAnalysisModuleFile << "      set(flipTemplate ${flip})" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
  // BMSShapeAnalysisModuleFile<<"      echo ('regTemplate: '${regTemplate})"<<std::endl;
  // BMSShapeAnalysisModuleFile<<"      echo ('flipTemplate: '${flipTemplate})"<<std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;

  if( GetParaOut1State() == true )
    {
    BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.writePara 1)" << std::endl;
    }
  BMSShapeAnalysisModuleFile << "    Run(output ${ParaT} error)" << std::endl;
  BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
  BMSShapeAnalysisModuleFile << "    echo('ParaToSPHARMMesh: '${error})" << std::endl;
  // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
  BMSShapeAnalysisModuleFile << " DeleteFile(" << GetOutputDirectory()
                             << "/Template/${basename}_pp_surf_paraPhiHalf.txt)" << std::endl;
  BMSShapeAnalysisModuleFile << " DeleteFile(" << GetOutputDirectory()
                             << "/Template/${basename}_pp_surf_paraMix.txt)" << std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${value})" << std::endl;

  if( GetRegTemplateState() == false && GetFlipTemplateState() == false ) // TODO
    {

    if( GetTemplateMState() == true )
      {
      if( MeanTemplateExist() == 0 )
        {
        BMSShapeAnalysisModuleFile << "    Set(value 0)" << std::endl;
        BMSShapeAnalysisModuleFile << "    foreach(data ${testPara2})" << std::endl;
        BMSShapeAnalysisModuleFile << "      Inc(${value} 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "      Int(${value})" << std::endl;
        BMSShapeAnalysisModuleFile << "    EndForEach(data)" << std::endl;
        BMSShapeAnalysisModuleFile << "    if (${value} < 4) " << std::endl;
        BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *SPHARM.vtk)" << std::endl;
        BMSShapeAnalysisModuleFile << "      listFileInDir(flip ${tdir} *SPHARM.coef)" << std::endl;
        BMSShapeAnalysisModuleFile << "      set(regTemplate ${reg})" << std::endl;
        BMSShapeAnalysisModuleFile << "      set(flipTemplate ${flip})" << std::endl;
        BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
        BMSShapeAnalysisModuleFile << "      echo ('regTemplate: '${regTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "      echo ('flipTemplate: '${flipTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
        BMSShapeAnalysisModuleFile << "      #The Template is the first data" << std::endl;
        BMSShapeAnalysisModuleFile << "      Set(value 0)" << std::endl;
        BMSShapeAnalysisModuleFile << "      foreach(data ${testPara})" << std::endl;
        BMSShapeAnalysisModuleFile << "        Inc(${value} 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "        Int(${value})" << std::endl;
        BMSShapeAnalysisModuleFile << "      EndForEach(data)" << std::endl;
        BMSShapeAnalysisModuleFile << "      if (${value} < 5)" << std::endl;
        BMSShapeAnalysisModuleFile << "        SetApp(Para @ParaToSPHARMMeshCLP)" << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"
                                   << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"
                                   << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf)"
                                   << std::endl;
        BMSShapeAnalysisModuleFile << "	     SetAppOption(Para.flipTemplateFileOn 1)" << std::endl;
        BMSShapeAnalysisModuleFile
        << "        SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.subdivLevel 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.subdivLevel.subdivLevel "
                                   << /*m_Parameters.*/ GetSubdivLevel() << ")" << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.spharmDegree 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.spharmDegree.spharmDegree "
                                   << /*m_Parameters.*/ GetSPHARMDegree() << ")" << std::endl;
		  BMSShapeAnalysisModuleFile << "    	  SetAppOption(Para.thetaIteration 1)" << std::endl;
		  BMSShapeAnalysisModuleFile << "    	  SetAppOption(Para.thetaIteration.thetaIteration " << GetThetaIteration() << ")"
				  << std::endl;
		  BMSShapeAnalysisModuleFile << "    	  SetAppOption(Para.phiIteration 1)" << std::endl;
		  BMSShapeAnalysisModuleFile << "    	  SetAppOption(Para.phiIteration.phiIteration " << GetPhiIteration() << ")"
				  << std::endl;
		  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh 1)" << std::endl;
		  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh.medialMesh " << GetPhiIteration() << ")"
				  << std::endl;
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.regTemplateFileOn 1)" << std::endl;
        BMSShapeAnalysisModuleFile
        << "        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex 1)" << std::endl;
        if( GetFinalFlipN() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 0 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipX() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 4 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipY() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 5 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 7 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipXY() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 1 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipYZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 2 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipXZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 3 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipXYZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 6 << ")"
                                     << std::endl;
          }

        if( GetParaOut1State() == true )
          {
          BMSShapeAnalysisModuleFile << "        SetAppOption(Para.writePara 1)" << std::endl;
          }
        BMSShapeAnalysisModuleFile << "        Run(output ${Para} error)" << std::endl;
        BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
        BMSShapeAnalysisModuleFile << "    echo('ParaToSPHARMMesh: '${error})" << std::endl;
        // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
        BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
        BMSShapeAnalysisModuleFile << "      endif(${value})" << std::endl;
        BMSShapeAnalysisModuleFile << "    endif(${value})" << std::endl;
        ComputeMean = 1;
        }
      if( MeanTemplateExist() != 0 )
        {
        ComputeMean = 0;
        BMSShapeAnalysisModuleFile << "    listFileInDir(deletelist ${SPHARMdir} ${basename}*SPHARM*) " << std::endl;
        BMSShapeAnalysisModuleFile << "    Set(value 0)" << std::endl;
        BMSShapeAnalysisModuleFile << "    foreach(data ${testPara2})" << std::endl;
        BMSShapeAnalysisModuleFile << "      Inc(${value} 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "      Int(${value})" << std::endl;
        BMSShapeAnalysisModuleFile << "    EndForEach(data)" << std::endl;
        BMSShapeAnalysisModuleFile << "    if (${value} < 4) " << std::endl;
        BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *meanAll.vtk)" << std::endl;
        BMSShapeAnalysisModuleFile << "      set(regTemplate ${reg})" << std::endl;
        BMSShapeAnalysisModuleFile << "      listFileInDir(flip ${tdir} *SPHARM.coef)" << std::endl;
        BMSShapeAnalysisModuleFile << "      set(flipTemplate ${flip})" << std::endl;
        BMSShapeAnalysisModuleFile << "      echo()" << std::endl;
        BMSShapeAnalysisModuleFile << "      echo ('regTemplate: '${regTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "      echo ('flipTemplate: '${flipTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "      echo()" << std::endl;
        BMSShapeAnalysisModuleFile << "      #The Template is the Mean" << std::endl;
        BMSShapeAnalysisModuleFile << "      listFileInDir(testPara ${SPHARMdir} ${basename}*SPHARM*)" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetApp(Para @ParaToSPHARMMeshCLP)" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"
                                   << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"
                                   << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf_tMean)"
                                   << std::endl;
        BMSShapeAnalysisModuleFile << "	   SetAppOption(Para.flipTemplateFileOn 1)" << std::endl;
        BMSShapeAnalysisModuleFile
        << "      SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.subdivLevel 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.subdivLevel.subdivLevel "
                                   << /*m_Parameters.*/ GetSubdivLevel() << ")" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.spharmDegree 1)" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.spharmDegree.spharmDegree "
                                   << /*m_Parameters.*/ GetSPHARMDegree() << ")" << std::endl;
		  BMSShapeAnalysisModuleFile << "      SetAppOption(Para.thetaIteration 1)" << std::endl;
		  BMSShapeAnalysisModuleFile << "      SetAppOption(Para.thetaIteration.thetaIteration " << GetThetaIteration() << ")"
				  << std::endl;
		  BMSShapeAnalysisModuleFile << "      SetAppOption(Para.phiIteration 1)" << std::endl;
		  BMSShapeAnalysisModuleFile << "      SetAppOption(Para.phiIteration.phiIteration " << GetPhiIteration() << ")"
				  << std::endl;
		  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh 1)" << std::endl;
		  BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh.medialMesh " << GetPhiIteration() << ")"
				  << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.regTemplateFileOn 1)" << std::endl;
        BMSShapeAnalysisModuleFile
        << "      SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})" << std::endl;
        BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex 1)" << std::endl;
        if( GetFinalFlipN() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 0 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipX() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 4 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipY() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 5 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 7 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipXY() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 1 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipYZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 2 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipXZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 3 << ")"
                                     << std::endl;
          }
        if( GetFinalFlipXYZ() == 1 )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 6 << ")"
                                     << std::endl;
          }

        if( GetParaOut1State() == true )
          {
          BMSShapeAnalysisModuleFile << "      SetAppOption(Para.writePara 1)" << std::endl;
          }
        BMSShapeAnalysisModuleFile << "      Run(output ${Para} error)" << std::endl;
        BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
        BMSShapeAnalysisModuleFile << "    echo('ParaToSPHARMMesh: '${error})" << std::endl;
        // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
        BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
        BMSShapeAnalysisModuleFile << "    endif(${value})" << std::endl;
        }
      }
    if( GetTemplateMState() == false )
      {
      BMSShapeAnalysisModuleFile << "    Set(value 0)" << std::endl;
      BMSShapeAnalysisModuleFile << "    foreach(data ${testPara2})" << std::endl;
      BMSShapeAnalysisModuleFile << "      Inc(${value} 1)" << std::endl;
      BMSShapeAnalysisModuleFile << "      Int(${value})" << std::endl;
      BMSShapeAnalysisModuleFile << "    EndForEach(data)" << std::endl;
      BMSShapeAnalysisModuleFile << "    if (${value} < 4) " << std::endl;
      BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *SPHARM.vtk)" << std::endl;
      BMSShapeAnalysisModuleFile << "      listFileInDir(flip ${tdir} *SPHARM.coef)" << std::endl;
      BMSShapeAnalysisModuleFile << "      set(regTemplate ${reg})" << std::endl;
      BMSShapeAnalysisModuleFile << "      set(flipTemplate ${flip})" << std::endl;
      BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
      BMSShapeAnalysisModuleFile << "      echo ('regTemplate: '${regTemplate})" << std::endl;
      BMSShapeAnalysisModuleFile << "      echo ('flipTemplate: '${flipTemplate})" << std::endl;
      BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
      BMSShapeAnalysisModuleFile << "      #The Template is the first data" << std::endl;
      BMSShapeAnalysisModuleFile << "      Set(value 0)" << std::endl;
      BMSShapeAnalysisModuleFile << "      foreach(data ${testPara})" << std::endl;
      BMSShapeAnalysisModuleFile << "        Inc(${value} 1)" << std::endl;
      BMSShapeAnalysisModuleFile << "        Int(${value})" << std::endl;
      BMSShapeAnalysisModuleFile << "      EndForEach(data)" << std::endl;
      BMSShapeAnalysisModuleFile << "      if (${value} < 4)" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetApp(Para @ParaToSPHARMMeshCLP)" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"
                                 << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"
                                 << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf)" << std::endl;
      BMSShapeAnalysisModuleFile << "	     SetAppOption(Para.flipTemplateFileOn 1)" << std::endl;
      BMSShapeAnalysisModuleFile
      << "        SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.subdivLevel 1)" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.subdivLevel.subdivLevel " << GetSubdivLevel() << ")"
                                 << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.spharmDegree 1)" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.spharmDegree.spharmDegree " << GetSPHARMDegree()
                                 << ")" << std::endl;
		BMSShapeAnalysisModuleFile << "        SetAppOption(Para.thetaIteration 1)" << std::endl;
		BMSShapeAnalysisModuleFile << "        SetAppOption(Para.thetaIteration.thetaIteration " << GetThetaIteration() << ")"
				<< std::endl;
		BMSShapeAnalysisModuleFile << "        SetAppOption(Para.phiIteration 1)" << std::endl;
		BMSShapeAnalysisModuleFile << "        SetAppOption(Para.phiIteration.phiIteration " << GetPhiIteration() << ")"
				<< std::endl;
		BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh 1)" << std::endl;
		BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh.medialMesh " << GetPhiIteration() << ")"
				<< std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.regTemplateFileOn 1)" << std::endl;
      BMSShapeAnalysisModuleFile
      << "        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex 1)" << std::endl;
      if( GetFinalFlipN() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 0 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipX() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 4 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipY() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 5 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipZ() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 7 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipXY() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 1 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipYZ() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 2 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipXZ() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 3 << ")"
                                   << std::endl;
        }
      if( GetFinalFlipXYZ() == 1 )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 6 << ")"
                                   << std::endl;
        }

      if( GetParaOut1State() == true )
        {
        BMSShapeAnalysisModuleFile << "        SetAppOption(Para.writePara 1)" << std::endl;
        }
      BMSShapeAnalysisModuleFile << "        Run(output ${Para} error)" << std::endl;
      BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
      BMSShapeAnalysisModuleFile << "    echo('ParaToSPHARMMesh: '${error})" << std::endl;
      // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
      BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
      BMSShapeAnalysisModuleFile << "      endif(${value})" << std::endl;
      BMSShapeAnalysisModuleFile << "    endif(${value})" << std::endl;
      }
    }

  if( (GetRegTemplateState() == true &&
       GetFlipTemplateState() == true) ||
      (GetRegTemplateState() == true &&
       GetFlipTemplateState() == false) || (GetRegTemplateState() == false && GetFlipTemplateState() == true) )

    {

    BMSShapeAnalysisModuleFile << "    Set(value 0)" << std::endl;
    BMSShapeAnalysisModuleFile << "    foreach(data ${testPara2})" << std::endl;
    BMSShapeAnalysisModuleFile << "      Inc(${value} 1)" << std::endl;
    BMSShapeAnalysisModuleFile << "      Int(${value})" << std::endl;
    BMSShapeAnalysisModuleFile << "    EndForEach(data)" << std::endl;
    BMSShapeAnalysisModuleFile << "    if (${value} < 4) " << std::endl;
    if( GetTemplateMState() == true &&  MeanTemplateExist() != 0 )
      {
      BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *meanAll.vtk)" << std::endl;
      }
    else
      {
      BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *SPHARM.vtk)" << std::endl;
      }
    BMSShapeAnalysisModuleFile << "      listFileInDir(flip ${tdir} *SPHARM.coef)" << std::endl;
    BMSShapeAnalysisModuleFile << "      set(regTemplate ${reg})" << std::endl;
    BMSShapeAnalysisModuleFile << "      set(flipTemplate ${flip})" << std::endl;
    BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
    if( GetRegTemplateState() == false  )
      {
      BMSShapeAnalysisModuleFile << "      echo ('regTemplate: '${regTemplate})" << std::endl;
      }
    if( GetFlipTemplateState() == false )
      {
      BMSShapeAnalysisModuleFile << "      echo ('flipTemplate: '${flipTemplate})" << std::endl;
      }
    BMSShapeAnalysisModuleFile << "  echo()" << std::endl;

    BMSShapeAnalysisModuleFile << "  #The Template is choosed by the user" << std::endl;
    BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
    if( GetRegTemplateState() == true )
      {
      BMSShapeAnalysisModuleFile << "  echo ('regTemplate: '" << GetRegTemplate() << ")" << std::endl;
      }
    if( GetFlipTemplateState() == true )
      {
      BMSShapeAnalysisModuleFile << "  echo ('flipTemplate: '" << GetFlipTemplate() << ")" << std::endl;
      }
    BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
    BMSShapeAnalysisModuleFile << "  Set(value 0)" << std::endl;
    BMSShapeAnalysisModuleFile << "  foreach(data ${testPara})" << std::endl;
    BMSShapeAnalysisModuleFile << "    Inc(${value} 1)" << std::endl;
    BMSShapeAnalysisModuleFile << "    Int(${value})" << std::endl;
    BMSShapeAnalysisModuleFile << "  EndForEach(data)" << std::endl;
    BMSShapeAnalysisModuleFile << "  if (${value} < 4)" << std::endl;

    BMSShapeAnalysisModuleFile << "SetApp(Para @ParaToSPHARMMeshCLP)" << std::endl;
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)" << std::endl;
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.inParaFile  ${SPHARMdir}${basename}_pp_para.vtk)" << std::endl;

    if( MeanTemplateExist() != 0 &&  MeanTemplateExist() != 0 )
      {
      BMSShapeAnalysisModuleFile << "      SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf_tMean)"
                                 << std::endl;
      }
    else
      {
      BMSShapeAnalysisModuleFile << "SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf)" << std::endl;
      }
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.subdivLevel 1)" << std::endl;
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.subdivLevel.subdivLevel " << GetSubdivLevel() << ")" << std::endl;
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.spharmDegree 1)" << std::endl;
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.spharmDegree.spharmDegree " << GetSPHARMDegree() << ")"
                               << std::endl;
	 BMSShapeAnalysisModuleFile << "SetAppOption(Para.thetaIteration 1)" << std::endl;
	 BMSShapeAnalysisModuleFile << "SetAppOption(Para.thetaIteration.thetaIteration " << GetThetaIteration() << ")"
			 << std::endl;
	 BMSShapeAnalysisModuleFile << "SetAppOption(Para.phiIteration 1)" << std::endl;
	 BMSShapeAnalysisModuleFile << "SetAppOption(Para.phiIteration.phiIteration " << GetPhiIteration() << ")"
			 << std::endl;
	 BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh 1)" << std::endl;
	 BMSShapeAnalysisModuleFile << "    SetAppOption(ParaT.medialMesh.medialMesh " << GetPhiIteration() << ")"
			 << std::endl;
    if( GetFlipTemplateState() == true )
      {
      BMSShapeAnalysisModuleFile << "SetAppOption(Para.flipTemplateFileOn 1)" << std::endl;
      BMSShapeAnalysisModuleFile << "SetAppOption(Para.flipTemplateFile.flipTemplateFile " << GetFlipTemplate()
                                 << ")" << std::endl;
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.regTemplateFileOn 1)" << std::endl;
      BMSShapeAnalysisModuleFile
      << "        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})" << std::endl;
      }
    if( GetRegTemplateState() == true )
      {
      BMSShapeAnalysisModuleFile << "SetAppOption(Para.flipTemplateFileOn 1)" << std::endl;
      BMSShapeAnalysisModuleFile
      << "      SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})" << std::endl;
      BMSShapeAnalysisModuleFile << "SetAppOption(Para.regTemplateFileOn 1)" << std::endl;
      if( MeanTemplateExist() != 0 &&  MeanTemplateExist() != 0 )
        {
        BMSShapeAnalysisModuleFile
        << "        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})" << std::endl;
        }
      else
        {
        BMSShapeAnalysisModuleFile << "SetAppOption(Para.regTemplateFile.regTemplateFile " << GetRegTemplate()
                                   << ")" << std::endl;
        }
      }
    BMSShapeAnalysisModuleFile << "SetAppOption(Para.finalFlipIndex 1)" << std::endl;
    if( GetFinalFlipN() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 0 << ")" << std::endl;
      }
    if( GetFinalFlipX() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 4 << ")" << std::endl;
      }
    if( GetFinalFlipY() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 5 << ")" << std::endl;
      }
    if( GetFinalFlipZ() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 7 << ")" << std::endl;
      }
    if( GetFinalFlipXY() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 1 << ")" << std::endl;
      }
    if( GetFinalFlipYZ() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 2 << ")" << std::endl;
      }
    if( GetFinalFlipXZ() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 3 << ")" << std::endl;
      }
    if( GetFinalFlipXYZ() == 1 )
      {
      BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 6 << ")" << std::endl;
      }

    if( GetParaOut1State() == true )
      {
      BMSShapeAnalysisModuleFile << "    SetAppOption(Para.writePara 1)" << std::endl;
      }
    BMSShapeAnalysisModuleFile << "    Run(output ${Para} error)" << std::endl;
    BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
    BMSShapeAnalysisModuleFile << "    echo('ParaToSPHARMMesh: '${error})" << std::endl;
    // BMSShapeAnalysisModuleFile<<"    exit()"<<std::endl;
    BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
    }

  BMSShapeAnalysisModuleFile << "Set(paraPhiHalf " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surf_paraPhiHalf.txt)" << std::endl;
  BMSShapeAnalysisModuleFile << "FileExists(testParaPhiHalf ${paraPhiHalf})" << std::endl;
  BMSShapeAnalysisModuleFile << "If( ${testParaPhiHalf} == 1 )" << std::endl;
  BMSShapeAnalysisModuleFile << "echo(${paraPhiHalf})" << std::endl;
  BMSShapeAnalysisModuleFile << " DeleteFile(${paraPhiHalf})" << std::endl;
  BMSShapeAnalysisModuleFile << "EndIf(${testParaPhiHalf})" << std::endl;

  BMSShapeAnalysisModuleFile << "Set(paraMix " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surf_paraMix.txt)" << std::endl;
  BMSShapeAnalysisModuleFile << "FileExists(testParaMix ${paraMix})" << std::endl;
  BMSShapeAnalysisModuleFile << "If( ${testParaMix} == 1 )" << std::endl;
  BMSShapeAnalysisModuleFile << " DeleteFile(${paraMix})" << std::endl;
  BMSShapeAnalysisModuleFile << "EndIf(${testParaMix})" << std::endl;
  // Write the Ouput File
  BMSShapeAnalysisModuleFile << "  #Write Ouput File" << std::endl;
  for( unsigned int i = 0; i < OutputFileHeaders.size(); i++ )
    {
    BMSShapeAnalysisModuleFile << "  GetParam(header" << i << " ${Header" << i << "} ${count})" << std::endl;
    BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${header" << i << "} ',' )" << std::endl;
    }
  BMSShapeAnalysisModuleFile << "  Inc(${count} 1)" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(postProcessFile " << GetOutputDirectory()
                             << "/Mesh/PostProcess/${basename}_pp" << GetVolumeFileExtension() << ")" << std::endl;
  BMSShapeAnalysisModuleFile << "  FileExists(testPostProcessFile ${postProcessFile})" << std::endl;
  BMSShapeAnalysisModuleFile << "  If(${testPostProcessFile} == 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${postProcessFile} ',')" << std::endl;
  BMSShapeAnalysisModuleFile << "  EndIf(${testPostProcessFile})" << std::endl;
  BMSShapeAnalysisModuleFile << "  If(${testPostProcessFile} == 0)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} 'none,')" << std::endl;
  BMSShapeAnalysisModuleFile << "  EndIf(${testPostProcessFile})" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(originSurf " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_para.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${originSurf} ',')" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(surfSPHARM " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surfSPHARM.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${surfSPHARM} ',')" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(surfSPHARMcoef " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surfSPHARM.coef)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${surfSPHARMcoef} ',')" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(surfSPHARM_ellalign " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surfSPHARM_ellalign.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${surfSPHARM_ellalign} ',')" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(surfSPHARM_ellalignCoef " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surfSPHARM_ellalign.coef)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${surfSPHARM_ellalignCoef} ',')" << std::endl;

  BMSShapeAnalysisModuleFile << "  Set(procalignFile " << GetOutputDirectory()
                             << "/Mesh/SPHARM/${basename}_pp_surfSPHARM_procalign.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "  FileExists(testProcalignFile ${procalignFile})" << std::endl;
  BMSShapeAnalysisModuleFile << "  If(${testProcalignFile} == 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} ${procalignFile} )" << std::endl;
  BMSShapeAnalysisModuleFile << "  EndIf(${testProcalignFile})" << std::endl;
  BMSShapeAnalysisModuleFile << "  If(${testProcalignFile} == 0)" << std::endl;
  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} 'none')" << std::endl;
  BMSShapeAnalysisModuleFile << "  EndIf(${testProcalignFile})" << std::endl;

  BMSShapeAnalysisModuleFile << "  appendFile(${OutputFile} '\\n' )" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "EndForEach(case)" << std::endl;

  /*BMSShapeAnalysisModuleFile<<" Set(BMSFile ${BMSdir}"<<"/ShapeAnalysisModule.bms"<<")"<<std::endl;
  BMSShapeAnalysisModuleFile<<" CopyFile("<<GetBMSShapeAnalysisModuleFile()<<" ${BMSdir})"<<std::endl;
  BMSShapeAnalysisModuleFile<<" DeleteFile("<<GetBMSShapeAnalysisModuleFile()<<")"<<std::endl;*/

  return;
}

// Write BMS file for mean computation
void ShapeAnalysisModuleComputation::WriteBMSShapeAnalysisModuleFile2()
{

  std::ofstream BMSShapeAnalysisModuleFile(GetBMSShapeAnalysisModuleFile2() );

  BMSShapeAnalysisModuleFile << "#----- Shape Analysis Computing ParaToSPHARMMesh with mean as input -----"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "#------------------------------------------------------------------------"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "#------------------------------------------------------------------------"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "#------------------------------------------------------------------------"
                             << std::endl;

  BMSShapeAnalysisModuleFile << "echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "echo('Computing with mean as a template')" << std::endl;
  BMSShapeAnalysisModuleFile << "echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "set (OrigCasesList '"
                             << GetNthDataListValue(1, GetColumnVolumeFile() ) << "')" << std::endl;

  int DataNumber;
  for( DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
    {
    BMSShapeAnalysisModuleFile << "set (OrigCasesList ${OrigCasesList} '" << GetNthDataListValue(DataNumber,
                                                                                                 GetColumnVolumeFile() )
                               << "')" << std::endl;
    }

  BMSShapeAnalysisModuleFile << "set(SPHARMdir '" << GetOutputDirectory() << "/Mesh/SPHARM/')" << std::endl;
  BMSShapeAnalysisModuleFile << "set(tdir '" << GetOutputDirectory() << "/Template/')" << std::endl;

  BMSShapeAnalysisModuleFile << "ForEach(case ${OrigCasesList})" << std::endl;
  BMSShapeAnalysisModuleFile << "  #Extract basename" << std::endl;
  BMSShapeAnalysisModuleFile << "  GetFilename(basename ${case} NAME_WITHOUT_EXTENSION)" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "  echo('Case: '${case})" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testPara ${SPHARMdir} *${basename}*SPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testPara2 ${SPHARMdir} *${basename}*MeanSPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testtemp ${tdir}  *SPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "  listFileInDir(testtemp2 ${tdir} *${basename}*SPHARM*)" << std::endl;

  BMSShapeAnalysisModuleFile << "    if (${value} < 4) " << std::endl;
  BMSShapeAnalysisModuleFile << "      listFileInDir(reg ${tdir} *meanAll.vtk)" << std::endl;
  BMSShapeAnalysisModuleFile << "      set(regTemplate ${reg})" << std::endl;
  BMSShapeAnalysisModuleFile << "      listFileInDir(flip ${tdir} *SPHARM.coef)" << std::endl;
  BMSShapeAnalysisModuleFile << "      set(flipTemplate ${flip})" << std::endl;
  BMSShapeAnalysisModuleFile << "      echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "      echo ('regTemplate : '${regTemplate})" << std::endl;
  BMSShapeAnalysisModuleFile << "      echo ('flipTemplate : '${flipTemplate})" << std::endl;
  BMSShapeAnalysisModuleFile << "      echo()" << std::endl;
  BMSShapeAnalysisModuleFile << "      #The Template is the Mean" << std::endl;
  BMSShapeAnalysisModuleFile << "      listFileInDir(testPara ${SPHARMdir} ${basename}*SPHARM*)" << std::endl;
  BMSShapeAnalysisModuleFile << "       SetApp(Para @ParaToSPHARMMeshCLP)" << std::endl;

  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.inSurfFile ${SPHARMdir}${basename}_pp_surf.vtk)"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.inParaFile ${SPHARMdir}${basename}_pp_para.vtk)"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.outbase  ${SPHARMdir}${basename}_pp_surf_tMean)"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "	     SetAppOption(Para.flipTemplateFileOn 1)" << std::endl;
  BMSShapeAnalysisModuleFile
  << "        SetAppOption(Para.flipTemplateFile.flipTemplateFile ${tdir}${flipTemplate})" << std::endl;

  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.subdivLevel 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.subdivLevel.subdivLevel " << GetSubdivLevel() << ")"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.spharmDegree 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.spharmDegree.spharmDegree " << GetSPHARMDegree() << ")"
                             << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.thetaIteration 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.thetaIteration.thetaIteration " << GetThetaIteration() << ")"
		  << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.phiIteration 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.phiIteration.phiIteration " << GetPhiIteration() << ")"
		  << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Para.medialMesh 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "    SetAppOption(Para.medialMesh.medialMesh " << GetPhiIteration() << ")"
		  << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.regTemplateFileOn 1)" << std::endl;
  BMSShapeAnalysisModuleFile << "        SetAppOption(Para.regTemplateFile.regTemplateFile ${tdir}${regTemplate})"
                             << std::endl;

  BMSShapeAnalysisModuleFile << "SetAppOption(Para.finalFlipIndex 1)" << std::endl;

  if( GetFinalFlipN() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 0 << ")" << std::endl;
    }
  if( GetFinalFlipX() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 4 << ")" << std::endl;
    }
  if( GetFinalFlipY() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 5 << ")" << std::endl;
    }
  if( GetFinalFlipZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 7 << ")" << std::endl;
    }
  if( GetFinalFlipXY() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 1 << ")" << std::endl;
    }
  if( GetFinalFlipYZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 2 << ")" << std::endl;
    }
  if( GetFinalFlipXZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 3 << ")" << std::endl;
    }
  if( GetFinalFlipXYZ() == 1 )
    {
    BMSShapeAnalysisModuleFile << "        SetAppOption(Para.finalFlipIndex.finalFlipIndex " << 6 << ")" << std::endl;
   }

  if( GetParaOut1State() == true )
    {
    BMSShapeAnalysisModuleFile << "      SetAppOption(Para.writePara 1)" << std::endl;
    }



  BMSShapeAnalysisModuleFile << "      Run(output ${Para} error)" << std::endl;
  BMSShapeAnalysisModuleFile << "  if(${error} != '')" << std::endl;
  BMSShapeAnalysisModuleFile << "    echo('ParaToSPHARMMesh: '${error})" << std::endl;
  BMSShapeAnalysisModuleFile << "  endif(${error})" << std::endl;
  BMSShapeAnalysisModuleFile << "    endif(${value})" << std::endl;

  BMSShapeAnalysisModuleFile << "EndForEach(${case})" << std::endl;

  BMSShapeAnalysisModuleFile << "set(BMSdir '" << GetOutputDirectory() << "/BatchMake_Scripts')" << std::endl;
  BMSShapeAnalysisModuleFile << " Set(BMSFile ${BMSdir}" << "/ShapeAnalysisModule_MeanAsTemplate.bms" << ")"
                             << std::endl;
  BMSShapeAnalysisModuleFile << " CopyFile(" << GetBMSShapeAnalysisModuleFile() << " ${BMSdir})" << std::endl;
  BMSShapeAnalysisModuleFile << " DeleteFile(" << GetBMSShapeAnalysisModuleFile() << ")" << std::endl;

}

// Overwrite (delete) the different files
void ShapeAnalysisModuleComputation::OverWrite()
{
  int     DataNumber = GetDataNumber();
  char * *Base_Files;

  Base_Files = new char *[DataNumber];
  for( int i = 0; i < DataNumber; i++ )
    {
    Base_Files[i] = new char[512];
    std::strcpy(Base_Files[i], GetAllFilesName(i).c_str() );
    }

  // Overwrite files created by SegPostProcess
  if( GetOverwriteSegPostProcess() )
    {
    char * *PreProcessing_Files;
    PreProcessing_Files = new char *[DataNumber];

    char *dir;
    dir = new  char[512];
    std::strcpy(dir, GetOutputDirectory() );
    std::strcat(dir, "/Mesh/PostProcess/");
    for( int i = 0; i < DataNumber; i++ )
      {
      PreProcessing_Files[i] = new char[512];
      std::strcpy(PreProcessing_Files[i], dir);
      std::strcat(PreProcessing_Files[i], Base_Files[i]);
      std::strcat(PreProcessing_Files[i], "_pp.gipl.gz");
      }

    bool DirectoryEmpty = DirectoryIsEmpty(dir);
    if( !DirectoryEmpty )
      {
      cout << "OverWrite SegPostProcess" << endl;
      for( int k = 0; k < DataNumber; k++ )
        {
        remove(PreProcessing_Files[k]);
        // cerr<<"Error deleting file "<<PreProcessing_Files[k]<<endl;

        }
      }
    else
      {
      cerr << "You have to compute SegPostProcess before overwriting!!" << endl;
      }
    }

  // Overwrite files created by GenParaMesh
  if( GetOverwriteGenParaMesh() )
    {
    char * *GenParaMesh_Files;
    GenParaMesh_Files = new char *[2 * DataNumber];

    char *dir;
    dir = new  char[512];
    std::strcpy(dir, GetOutputDirectory() );
    std::strcat(dir, "/Mesh/SPHARM/");
    for( int i = 0; i < DataNumber; i++ )
      {
      GenParaMesh_Files[i] = new char[512];
      std::strcpy(GenParaMesh_Files[i], dir);
      std::strcat(GenParaMesh_Files[i], Base_Files[i]);
      std::strcat(GenParaMesh_Files[i], "_pp_para.vtk");

      GenParaMesh_Files[i + DataNumber] = new char[512];
      std::strcpy(GenParaMesh_Files[i + DataNumber], dir);
      std::strcat(GenParaMesh_Files[i + DataNumber], Base_Files[i]);
      std::strcat(GenParaMesh_Files[i + DataNumber], "_pp_surf.vtk");
      }

    bool DirectoryEmpty = DirectoryIsEmpty(dir);
    if( !DirectoryEmpty )
      {
      cout << "OverWrite GenParaMesh" << endl;
      for( int k = 0; k < 2 * DataNumber; k++ )
        {
        remove(GenParaMesh_Files[k]);
        // cerr<<"Error deleting file "<<GenParaMesh_Files[k]<<endl;
        }
      }
    else
      {
      cerr << "You have to compute GenParaMesh before overwriting!!" << endl;
      }
    }

  // Overwrite files created by ParaToSPHARMMesh
  if( GetOverwriteParaToSPHARMMesh() )
    {
    char * *ParaToSPHARMMesh_Files;
    ParaToSPHARMMesh_Files = new char *[7 * DataNumber];

    char *dir;
    dir = new  char[512];
    std::strcpy(dir, GetOutputDirectory() );
    std::strcat(dir, "/Mesh/SPHARM/");
    for( int i = 0; i < DataNumber; i++ )
      {
      ParaToSPHARMMesh_Files[i] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i], dir);
      std::strcat(ParaToSPHARMMesh_Files[i], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i], "_pp_surfSPHARM.vtk");

      ParaToSPHARMMesh_Files[i + DataNumber] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i + DataNumber], dir);
      std::strcat(ParaToSPHARMMesh_Files[i + DataNumber], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i + DataNumber], "_pp_surfSPHARM_ellalign.vtk");

      ParaToSPHARMMesh_Files[i + 2 * DataNumber] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i + 2 * DataNumber], dir);
      std::strcat(ParaToSPHARMMesh_Files[i + 2 * DataNumber], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i + 2 * DataNumber], "_pp_surfSPHARM.coef");

      ParaToSPHARMMesh_Files[i + 3 * DataNumber] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i + 3 * DataNumber], dir);
      std::strcat(ParaToSPHARMMesh_Files[i + 3 * DataNumber], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i + 3 * DataNumber], "_pp_surfSPHARM_ellalign.coef");

      ParaToSPHARMMesh_Files[i + 4 * DataNumber] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i + 4 * DataNumber], dir);
      std::strcat(ParaToSPHARMMesh_Files[i + 4 * DataNumber], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i + 4 * DataNumber], "_pp_surfSPHARM_procalign.vtk");

      ParaToSPHARMMesh_Files[i + 5 * DataNumber] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i + 5 * DataNumber], dir);
      std::strcat(ParaToSPHARMMesh_Files[i + 5 * DataNumber], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i + 5 * DataNumber], "_pp_surf_paraPhi.txt");

      ParaToSPHARMMesh_Files[i + 6 * DataNumber] = new char[512];
      std::strcpy(ParaToSPHARMMesh_Files[i + 6 * DataNumber], dir);
      std::strcat(ParaToSPHARMMesh_Files[i + 6 * DataNumber], Base_Files[i]);
      std::strcat(ParaToSPHARMMesh_Files[i + 6 * DataNumber], "_pp_surf_paraTheta.txt");

      }

    // check if the directory is empty
    bool DirectoryEmpty = DirectoryIsEmpty(dir);
    if( !DirectoryEmpty )
      {
      cout << "OverWrite ParaToSPARHMMesh" << endl;
      for( int k = 0; k < 7 * DataNumber; k++ )
        {
        remove(ParaToSPHARMMesh_Files[k]);
        // cerr<<"Error deleting file "<<ParaToSPHARMMesh_Files[k]<<endl;
        }
      }

    else
      {
      cerr << "You have to compute ParaToSPHARMMesh before overwriting!!" << endl;
      }

    // Overwrite the template file;
    char *dirTemplate = NULL;
    dirTemplate = new  char[512];
    std::strcpy(dirTemplate, GetOutputDirectory() );
    std::strcat(dirTemplate, "/Template/");

    bool templateDirectoryEmpty = DirectoryIsEmpty(dirTemplate);
    if( !templateDirectoryEmpty )
      {
      cout << "OverWrite Templates" << endl;
      int length = strlen(GetOutputDirectory() );
      length = length + 10;

    /*  DIR *          pdir = NULL;
      struct dirent *pent;
      pdir = opendir(dirTemplate);
      while( (pent = readdir(pdir) ) )
        {
        char *file = NULL;
        file = new char[512];
        strcpy(file, dirTemplate);
        strcat(file, pent->d_name);
        if( file[length] != '.' )
          {
          remove(file);
          // cerr<<"Error deleting file "<<file<<endl;
          }
        }*/
  unsigned long nbFiles= itksys:: Directory::GetNumberOfFilesInDirectory(dirTemplate);
   unsigned long i =0;
   
    for(i=0 ; i< nbFiles-1 ; i++)
    {
      itksys:: Directory direc ; 
    
      if (!direc.Load (dirTemplate))
      {
       cerr << "Directory  cannot be read!" << endl;
      }
      const char * FileDir=direc.GetFile(i) ;
     
      char *file = NULL;
      file = new char[512];
      strcpy(file, dirTemplate);
      strcat(file,FileDir);
      if( file[length] != '.' )
        {
        remove(file);
        // cerr<<"Error deleting file "<<file<<endl;
        }
    }


      }
    }
}

// Compute the mean of all the files
void ShapeAnalysisModuleComputation::ComputationMean()
{
  std::cout << "Computing Mean" << std::endl;
  int     DataNumber = GetDataNumber();
  char * *Base_Files;
  Base_Files = new char *[DataNumber];

  int          nbPoints = 0;
  bool         initialize = true;
  vtkPolyData *poly = vtkPolyData::New();
  for( int i = 0; i < DataNumber; i++ )
    {
    Base_Files[i] = new char[512];
    std::strcpy(Base_Files[i], "/");
    std::strcpy(Base_Files[i], GetOutputDirectory() );
    std::strcat(Base_Files[i], "/Mesh/SPHARM/");
    std::strcat(Base_Files[i], GetAllFilesName(i).c_str() );
    std::strcat(Base_Files[i], "_pp_surfSPHARM_procalign.vtk");

    // read a vtk file
    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(Base_Files[i]);

    try
      {
      meshin->Update();
      }
    catch( ... )
      {
      std::cout << "Cannot open file: " << Base_Files[i] << std::endl;
      return;
      }

    poly = meshin->GetOutput();

    // Get number of points
    vtkIdType idNumPointsInFile = poly->GetNumberOfPoints();
    nbPoints = idNumPointsInFile;

    vtkPoints * pts;
    pts = poly->GetPoints();

    // add all values of the coordinates
    if( initialize == true )
      {
      for( int _i = 0; _i < nbPoints; ++_i )
        {
        m_meanX.push_back(0);
        m_meanY.push_back(0);
        m_meanZ.push_back(0);
        }
      initialize = false;
      }
    for( int j = 0; j < nbPoints; ++j )
      {
       double *p;

       p = pts->GetPoint(j);
       m_meanX[j] = m_meanX[j] + p[0];
       m_meanY[j] = m_meanY[j] + p[1];
       m_meanZ[j] = m_meanZ[j] + p[2];
      }
    }
  for( int i = 0; i < nbPoints; ++i )
    {
    m_meanX[i] = m_meanX[i] / DataNumber;
    m_meanY[i] = m_meanY[i] / DataNumber;
    m_meanZ[i] = m_meanZ[i] / DataNumber;
    }

  char m_PolyFile[512];

  std::strcpy(m_PolyFile, GetOutputDirectory() );
  std::strcat(m_PolyFile, "/Template/");
  std::strcat(m_PolyFile, "polyAll.vtk");

  vtkCellArray* vtkcells = poly->GetPolys();
  vtkPolyData * polydata = vtkPolyData::New();
  polydata->SetPolys(vtkcells);
  vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
  writer->SetFileName(m_PolyFile);
  #if VTK_MAJOR_VERSION > 5
  writer->SetInputData(polydata);
  #else
  writer->SetInput(polydata);
  #endif
  writer->Write();

  WriteMeanFile(nbPoints);
}

// Write the mean file use as a template for ParaToSPHARMMesh
void ShapeAnalysisModuleComputation::WriteMeanFile(int nbPoints)
{
  char m_MeanFile[512]; char m_PolyFile[512];

  std::strcpy(m_MeanFile, GetOutputDirectory() );
  std::strcat(m_MeanFile, "/Template/");
  std::strcpy(m_PolyFile, m_MeanFile);
  std::strcat(m_MeanFile, "meanAll.vtk");
  std::strcat(m_PolyFile, "polyAll.vtk");

  std::string line;
  int         nbline = 0;

  std::ifstream tmpPolyFile(m_PolyFile);
  std::ofstream meanAll;
  meanAll.open(m_MeanFile, ios::out);
  meanAll << "# vtk DataFile Version 3.0" << std::endl;
  meanAll << "vtk output" << std::endl;
  meanAll << "ASCII" << std::endl;
  meanAll << "DATASET POLYDATA" << std::endl;
  meanAll << "POINTS " << nbPoints << " float" << std::endl;
  for( unsigned int i = 0; i < m_meanX.size(); i += 3 )
    {

    meanAll << m_meanX[i] << " " << m_meanY[i] << " " << m_meanZ[i] << " "
            << m_meanX[i + 1] << " " << m_meanY[i + 1] << " "
            << m_meanZ[i + 1] << " " << m_meanX[i + 2] << " " << m_meanY[i + 2] << " " << m_meanZ[i + 2] << endl;
    }
    meanAll << std::endl;

  while( getline(tmpPolyFile, line) )
    {
    if( nbline > 4 )
      {
      meanAll << line << std::endl;
      }
    nbline++;
    }

  meanAll.close();

  remove( m_PolyFile );

}

void ShapeAnalysisModuleComputation::GetRandomNum(int min, int max)
{
  // int min, max;

  srand(getpid() + time(NULL) );

  int nofNum = max - min + 1;
  // std::cout << nofNum << " " << max << " " << min << std::endl;

  // Create a random permutation
  int index, temp;
  for( int i = 0; i < nofNum; ++i )
    {
    m_permutations[i] = i;
    }
  for( int i = 1; i < nofNum; ++i )
    {
    index = i + (rand() % (nofNum - i) );
    temp = m_permutations[i];
    m_permutations[i] = m_permutations[index];
    m_permutations[index] = temp;
    }
  for( int i = 0; i < nofNum; ++i )
    {
    m_permutations[i]++;
    }
  // random swaps

  // print perm
  // cout << "permutation : ";
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
  std::cout << " " << std::endl;
  std::cout << "Run the Particle Correspondence" << std::endl;
  std::string particlesdirectory;
  std::string outputdir; outputdir.append(GetOutputDirectory() );
  particlesdirectory.append(outputdir);
  particlesdirectory.append("/ParticleCorrespondence");
  itksys::SystemTools::MakeDirectory(particlesdirectory.c_str() );

  std::string csvdirectory;
  csvdirectory.append(outputdir);
  csvdirectory.append("/OutputGroupFile/ShapeAnalysisModule_OutputFileVersion1.csv");

  // command line
  std::vector<const char *> args;
  args.push_back("ParticleModule" );
  args.push_back("--columMeshFile" );
  if( GetUseProcalign() )
    {
    args.push_back("9");
    }
  else
    {
    args.push_back("5");
    }                       // Original Space
  args.push_back("--sx" );
  args.push_back(Convert_Double_To_CharArray(GetEnforcedSpaceX() ) );
  args.push_back("--sy" );
  args.push_back(Convert_Double_To_CharArray(GetEnforcedSpaceY() ) );
  args.push_back("--sz" );
  args.push_back(Convert_Double_To_CharArray(GetEnforcedSpaceZ() ) );
  args.push_back("--HorizontalGridPara" );
  args.push_back(Convert_Double_To_CharArray(GetHorizontalGridPara() ) );
  args.push_back("--VerticalGridPara" );
  args.push_back(Convert_Double_To_CharArray(GetVerticalGridPara() ) );
  args.push_back(csvdirectory.c_str() );
  args.push_back(particlesdirectory.c_str() );
  args.push_back(0);
  for( unsigned int k = 0; k < args.size() - 1; k++ )
    {
    std::cout << args.at(k) << " ";
    }
  std::cout << " " << std::endl; std::cout << " " << std::endl;

  // itk sys parameters
  int    length;
  time_t start, end;
  time(&start);

  double timeout = 0.05;
  int    result;
  char*  dataitk = NULL;

  itksysProcess* gp = itksysProcess_New();
  itksysProcess_SetCommand(gp, &*args.begin() );
  itksysProcess_SetOption(gp, itksysProcess_Option_HideWindow, 1);
  itksysProcess_Execute(gp);
  while( int Value = itksysProcess_WaitForData(gp, &dataitk, &length, &timeout) ) // wait for 1s
    {
    if( ( (Value == itksysProcess_Pipe_STDOUT) || (Value == itksysProcess_Pipe_STDERR) ) && dataitk[0] == 'D' )
      {
      std::strstream st;
      for( int i = 0; i < length; ++i )
        {
        st << dataitk[i];
        }
      std::string dim = st.str();
      }
    time(&end);
    cout << "(processing since " << difftime(end, start) << " seconds) \r" << flush;
    timeout = 0.05;
    }

  itksysProcess_WaitForExit(gp, 0);
  result = 1;

  switch( itksysProcess_GetState(gp) )
    {
    case itksysProcess_State_Exited:
      {
      result = itksysProcess_GetExitValue(gp);
      } break;
    case itksysProcess_State_Error:
      {
      std::cerr << "Error: Could not run " << args[0] << ":\n";
      std::cerr << itksysProcess_GetErrorString(gp) << "\n";
      std::cout << "Error: Could not run " << args[0] << ":\n";
      std::cout << itksysProcess_GetErrorString(gp) << "\n";
      } break;
    case itksysProcess_State_Exception:
      {
      std::cerr << "Error: " << args[0] << " terminated with an exception: " << itksysProcess_GetExceptionString(gp)
                << "\n";
      std::cout << "Error: " << args[0] << " terminated with an exception: " << itksysProcess_GetExceptionString(gp)
                << "\n";
      } break;
    case itksysProcess_State_Starting:
    case itksysProcess_State_Executing:
    case itksysProcess_State_Expired:
    case itksysProcess_State_Killed:
      {
      // Should not get here.
      std::cerr << "Unexpected ending state after running " << args[0] << std::endl;
      std::cout << "Unexpected ending state after running " << args[0] << std::endl;
      } break;
    }
  itksysProcess_Delete(gp);

}

void ShapeAnalysisModuleComputation::WriteComputationLog()
{
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "NUMBER OF DATA PROCESSED " << GetDataNumber() << std::endl;
  std::cout << "------------------------------------------" << std::endl;

  std::cout << std::endl; 
  std::cout << ">>>>>>>>>>>>>>>>  " << "SEGPOSTPROCESS SUMMARY RESULTS" << std::endl;
  std::cout << std::endl;

  std::vector<int> data_flags(GetDataNumber(),0);
  for( int DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
    {
    	//FOR ALL FILES IN THE INPUT FILE, DO THE PERTINENT CHECKINGS 
    	//std::cout << DataNumber << " - " << GetNthDataListValue(DataNumber,GetColumnVolumeFile() ) << std::endl;
    	
	char *aux2, *aux;
	string input_file = GetNthDataListValue(DataNumber,GetColumnVolumeFile());
        aux = new char [input_file.size()+1];
	strcpy (aux, input_file.c_str());

	aux2 = strtok(aux,".");
        string PreProcess_file = aux2;
    	PreProcess_file.append("_pp.gipl.gz");
	
	string PostProcess_path = GetOutputDirectory();
	PostProcess_path.append("Mesh/PostProcess/");
	string filename_found;

        if ( itksys::SystemTools::LocateFileInDir(PreProcess_file.c_str(),PostProcess_path.c_str(),filename_found,0) )
	{
		if (GetDebug())
			std::cout << "SegPostProcess file computed correctly : " << filename_found << std::endl;
	}
	else
	{
		std::cout << "SegPostProcess file missing : " << PreProcess_file << " CHECK INPUT DATA" << std::endl;
		data_flags[DataNumber-1] = 1;
	}
    }

  std::cout << std::endl; 
  std::cout << ">>>>>>>>>>>>>>>>  " << "GENPARAMESH SUMMARY RESULTS" << std::endl;
  std::cout << std::endl;

  for( int DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
    {
    	//FOR ALL FILES IN THE INPUT FILE, DO THE PERTINENT CHECKINGS 
    	//std::cout << DataNumber << " - " << GetNthDataListValue(DataNumber,GetColumnVolumeFile() ) << std::endl;
    	if ( data_flags[DataNumber-1]==0 )
	{
		char *aux2, *aux;
		string input_file = GetNthDataListValue(DataNumber,GetColumnVolumeFile());
		aux = new char [input_file.size()+1];
		strcpy (aux, input_file.c_str());
	
		aux2 = strtok(aux,".");
		string GenParaMesh_path = GetOutputDirectory();
		GenParaMesh_path.append("Mesh/SPHARM/");
		
		// Find surf.vtk files
		string Surf_file = aux2;
		Surf_file.append("_pp_surf.vtk");
		string aux3 = GenParaMesh_path;
		aux3.append(itksys::SystemTools::GetFilenameName(Surf_file));
		Surf_file = aux3;
		string filename_found;
	
		if ( itksys::SystemTools::LocateFileInDir(Surf_file.c_str(),GenParaMesh_path.c_str(),filename_found,0) )
		{
			if (GetDebug())
				std::cout << "GenParaMesh surface file computed correctly : " << filename_found << std::endl;
		}
		else
		{
			std::cout << "GenParaMesh surface file is missing : " << Surf_file << std::endl;
			data_flags[DataNumber-1]=1;
		}
	
		// Find para.vtk files
		string Para_file = aux2;
		Para_file.append("_pp_para.vtk");
		aux3 = GenParaMesh_path;
		aux3.append(itksys::SystemTools::GetFilenameName(Para_file));
		Surf_file = aux3;
	
		if ( itksys::SystemTools::LocateFileInDir(Para_file.c_str(),GenParaMesh_path.c_str(),filename_found,0) )
		{
			if (GetDebug())
				std::cout << "GenParaMesh parameterization file computed correctly : " << filename_found << std::endl;
		}
		else
		{
			std::cout << "GenParaMesh parameterization file is missing : " << Para_file << std::endl;
			string Euler_path = GetOutputDirectory();
			Euler_path.append("EulerFiles/");
			data_flags[DataNumber-1]=1;
		}
	}
    }

    std::cout << std::endl; 
    std::cout << ">>>>>>>>>>>>>>>>  " << "PARATOSPHARMMESH SUMMARY RESULTS" << std::endl;
    std::cout << std::endl;

    for( int DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
    {
    	//FOR ALL FILES IN THE INPUT FILE, DO THE PERTINENT CHECKINGS 
	char *aux2, *aux;
	string input_file = GetNthDataListValue(DataNumber,GetColumnVolumeFile());
	aux = new char [input_file.size()+1];
	strcpy (aux, input_file.c_str());

	aux2 = strtok(aux,".");
	string GenParaMesh_path = GetOutputDirectory();
	GenParaMesh_path.append("Mesh/SPHARM/");

	// Find surfSPHARM files
	string surfSPHARM_file = aux2;
	surfSPHARM_file.append("_pp_surfSPHARM.vtk");
	string aux3 = GenParaMesh_path;
	aux3.append(itksys::SystemTools::GetFilenameName(surfSPHARM_file));
	surfSPHARM_file = aux3;
	string filename_found;

	if ( itksys::SystemTools::LocateFileInDir(surfSPHARM_file.c_str(),GenParaMesh_path.c_str(),filename_found,0) )
	{
		if (GetDebug())
		{
			std::cout << "ParaToSPHARMMesh surfSPHARM file computed correctly : " << filename_found << std::endl;
		}

		vtkPolyDataReader *VTKreader = vtkPolyDataReader::New();
		vtkPolyData *surfaceMesh = vtkPolyData::New();
		VTKreader->SetFileName(filename_found.c_str() );
		VTKreader->Update();
		surfaceMesh = VTKreader->GetOutput();

		if (surfaceMesh->GetNumberOfPoints() == 0 )
		{
			std::cout << filename_found << " has been computed, but has no points... deleting!" << std::endl;
			std::vector<const char *> args;
			args.push_back("rm");
			args.push_back(filename_found.c_str());
			args.push_back(0);

			// Run the application
			itksysProcess* gp = itksysProcess_New();
			itksysProcess_SetCommand(gp, &*args.begin() );
			itksysProcess_SetOption(gp, itksysProcess_Option_HideWindow, 1);
			itksysProcess_Execute(gp);
			data_flags[DataNumber-1]=1;
		} 
	}
	else
	{
		std::cout << "One or more surfSPHARM file(s) missing : " << surfSPHARM_file << std::endl;
		data_flags[DataNumber-1]=1;
	}
    }

    for( int DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
    {
    	//FOR ALL FILES IN THE INPUT FILE, DO THE PERTINENT CHECKINGS 
	char *aux2, *aux;
	string input_file = GetNthDataListValue(DataNumber,GetColumnVolumeFile());
	aux = new char [input_file.size()+1];
	strcpy (aux, input_file.c_str());

	aux2 = strtok(aux,".");
	string GenParaMesh_path = GetOutputDirectory();
	GenParaMesh_path.append("Mesh/SPHARM/");

	// Find surfSPHARM files
	string surfSPHARM_file = aux2;
	surfSPHARM_file.append("_pp_surfSPHARM_procalign.vtk");
	string aux3 = GenParaMesh_path;
	aux3.append(itksys::SystemTools::GetFilenameName(surfSPHARM_file));
	surfSPHARM_file = aux3;
	string filename_found;

	if ( itksys::SystemTools::LocateFileInDir(surfSPHARM_file.c_str(),GenParaMesh_path.c_str(),filename_found,0) )
	{
		if (GetDebug())
		{
			std::cout << "ParaToSPHARMMesh surfSPHARM file computed correctly : " << filename_found << std::endl;
		}

		vtkPolyDataReader *VTKreader = vtkPolyDataReader::New();
		vtkPolyData *surfaceMesh = vtkPolyData::New();
		VTKreader->SetFileName(filename_found.c_str() );
		VTKreader->Update();
		surfaceMesh = VTKreader->GetOutput();

		if (surfaceMesh->GetNumberOfPoints() == 0 )
		{
			std::cout << filename_found << " has been computed, but has no points... deleting!" << std::endl;
			std::vector<const char *> args;
			args.push_back("rm");
			args.push_back(filename_found.c_str());
			args.push_back(0);

			// Run the application
			itksysProcess* gp = itksysProcess_New();
			itksysProcess_SetCommand(gp, &*args.begin() );
			itksysProcess_SetOption(gp, itksysProcess_Option_HideWindow, 1);
			itksysProcess_Execute(gp);
			data_flags[DataNumber-1]=1;
		} 
	}
	else
	{
		std::cout << "One or more surfSPHARM file(s) missing : " << surfSPHARM_file << std::endl;
		data_flags[DataNumber-1]=1;
	}
    }

    for( int DataNumber = 2; DataNumber <= GetDataNumber(); DataNumber++ )
    {
    	//FOR ALL FILES IN THE INPUT FILE, DO THE PERTINENT CHECKINGS 
	char *aux2, *aux;
	string input_file = GetNthDataListValue(DataNumber,GetColumnVolumeFile());
	aux = new char [input_file.size()+1];
	strcpy (aux, input_file.c_str());

	aux2 = strtok(aux,".");
	string GenParaMesh_path = GetOutputDirectory();
	GenParaMesh_path.append("Mesh/SPHARM/");

	// Find surfSPHARM files
	string surfSPHARM_file = aux2;
	surfSPHARM_file.append("_pp_surfSPHARM_ellalign.vtk");
	string aux3 = GenParaMesh_path;
	aux3.append(itksys::SystemTools::GetFilenameName(surfSPHARM_file));
	surfSPHARM_file = aux3;
	string filename_found;

	if ( itksys::SystemTools::LocateFileInDir(surfSPHARM_file.c_str(),GenParaMesh_path.c_str(),filename_found,0) )
	{
		if (GetDebug())
		{
			std::cout << "ParaToSPHARMMesh surfSPHARM file computed correctly : " << filename_found << std::endl;
		}

		vtkPolyDataReader *VTKreader = vtkPolyDataReader::New();
		vtkPolyData *surfaceMesh = vtkPolyData::New();
		VTKreader->SetFileName(filename_found.c_str() );
		VTKreader->Update();
		surfaceMesh = VTKreader->GetOutput();

		if (surfaceMesh->GetNumberOfPoints() == 0 )
		{
			std::cout << filename_found << " has been computed, but has no points... deleting!" << std::endl;
			std::vector<const char *> args;
			args.push_back("rm");
			args.push_back(filename_found.c_str());
			args.push_back(0);

			// Run the application
			itksysProcess* gp = itksysProcess_New();
			itksysProcess_SetCommand(gp, &*args.begin() );
			itksysProcess_SetOption(gp, itksysProcess_Option_HideWindow, 1);
			itksysProcess_Execute(gp);
			data_flags[DataNumber-1]=1;
		} 
	}
	else
	{
		std::cout << "One or more surfSPHARM file(s) missing : " << surfSPHARM_file << std::endl;
		data_flags[DataNumber-1]=1;
	}
    }

    std::cout << std::endl; 
    std::cout << "------------------------------------------" << std::endl;
    std::cout << "CHECK INPUT DATA" << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    std::cout << std::endl;


    for ( unsigned long i = 0 ; i < data_flags.size() ; i ++ )
    {
	string input_file = GetNthDataListValue(i+1,GetColumnVolumeFile());
	if ( data_flags[i]==1 )
	{
		char *aux2, *aux;
		aux = new char [input_file.size()+1];
		strcpy (aux, input_file.c_str());
	
		aux2 = strtok(aux,".");
		string Euler_path = GetOutputDirectory();
		Euler_path.append("EulerFiles/");
	
		// Find Euler files
		string Euler_file = aux2;
		Euler_file.append("_euler.txt");
		string aux3 = Euler_path;
		aux3.append(itksys::SystemTools::GetFilenameName(Euler_file));
		Euler_file = aux3;

		std::ifstream input;
		char line[100];
    		//open input
		input.open ( Euler_file.c_str() ); 
		input.seekg(0,std::ios::beg);
		input.getline(line, 1000);
		float EulerValue;
		sscanf(line, "%f", &EulerValue);
		//close input	
		input.close () ;
	
		if (EulerValue == 0)
			std::cout << input_file  << " does not have spherical parameterization" << std::endl;
		else
			std::cout << input_file  << std::endl;
		
	}
    }

    std::cout << std::endl;
}
