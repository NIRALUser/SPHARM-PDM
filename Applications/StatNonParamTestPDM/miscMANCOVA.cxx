#include "miscMANCOVA.h"
#include <sstream>
#include <boost/math/distributions/students_t.hpp>
#include <itksys/SystemTools.hxx>

using namespace std;
vnl_vector<double> TiedRankCalc(vnl_vector<double>& data)
{
  // computes ranking of values in data vector, but resolves ties by
  // returning the average rank

  vnl_vector<double> output(data.size() );

  typedef double T;
  std::vector<std::pair<T, int> > pairVector;
  for( unsigned int i = 0; i < data.size(); ++i )
    {
    pairVector.push_back(std::pair<T, int>(data(i), i) );
    }

  std::sort(pairVector.begin(), pairVector.end() );

  // need to go through values and average if there are ties

  unsigned int lastIndex = 0;
  double       lastValue = (pairVector[lastIndex]).first;
  unsigned int currentIndex = 0;

  if( data.size() <= 1 ) // there is no data or only one data point, so just return
    {
    for( unsigned int i = 0; i < data.size(); ++i )
      {
      output(i) = i;
      }
    return output;
    }

  double nextValue = (pairVector[currentIndex + 1]).first;

  bool keepGoing = true;

  do
    {

    while( lastValue == nextValue ) // check for ties
      {
      currentIndex++;
      if( currentIndex + 1 >= data.size() ) // we are at the end
        {
        break;
        }
      else
        {
        nextValue = (pairVector[currentIndex + 1]).first; // check the
        // next value
        }
      }

    unsigned int nrOfTies = currentIndex - lastIndex + 1;
    double       averageRankSum = 0;
    for( unsigned int i = lastIndex; i <= currentIndex; ++i )
      {
      // averageRankSum += (pairVector[i]).second;
      averageRankSum += i;
      }

    averageRankSum /= nrOfTies;
    for( unsigned int i = lastIndex; i <= currentIndex; ++i )
      {
      output( (pairVector[i]).second) = averageRankSum;
      }

    if( currentIndex >= data.size() - 1 ) // we are at the end, so stop processing
      {
      keepGoing = false;
      }
    else
      {
      currentIndex++; // let's look at the next values
      lastIndex = currentIndex;
      lastValue = (pairVector[lastIndex]).first;

      if( currentIndex + 1 >= data.size() ) // we are now at the end, so
      // simply return the result
        {
        output( currentIndex ) = (pairVector[currentIndex]).second;
        keepGoing = false;
        }
      else
        { // there is a next value and it is the following ...
        nextValue = (pairVector[currentIndex + 1]).first;
        }
      }

    }
  while( keepGoing );    // that's all folks ...

  return output;

}

double computePearsonCorrelation( vnl_vector<double>& x, vnl_vector<double>& y )
{
  unsigned int numSubjects = x.size();

  // first compute the means

  double dXM = 0;
  double dYM = 0;

  for( unsigned int sub = 0; sub < numSubjects; sub++ )
    {
    dXM += x[sub];
    dYM += y[sub];
    }

  dXM /= numSubjects;
  dYM /= numSubjects;

  // now compute the variances

  double dXV = 0;
  double dYV = 0;
  for( unsigned int sub = 0; sub < numSubjects; sub++ )
    {
    dXV += (x[sub] - dXM) * (x[sub] - dXM);
    dYV += (y[sub] - dYM) * (y[sub] - dYM);
    }

  dXV /= (numSubjects - 1);
  dYV /= (numSubjects - 1);

  // which gives the standard deviations

  double dXS = sqrt( dXV );
  double dYS = sqrt( dYV );

  // from which we can compute the correlation coefficient

  double dR = 0;
  for( unsigned int sub = 0; sub < numSubjects; sub++ )
    {
    dR += (x[sub] - dXM) * (y[sub] - dYM);
    }
  dR /= (numSubjects - 1) * dXS * dYS;

  return dR;

}

double computePearsonCorrelationWithP( vnl_vector<double>& x, vnl_vector<double>& y, double & dPDiff,
                                       double & dPGreater,
                                       double & dPSmaller )
{

  // Equations from the book
  // Handbook of Parametric and Nonparametric Statistical Procedures
  // David J. Sheskin, Chapman & Hall/CRC
  //
  // as well as from the BOOST documentation

  unsigned int numSubjects = x.size();

  double dR = computePearsonCorrelation( x, y );

  // Degrees of freedom:
  unsigned v = numSubjects - 1;
  // t-statistic:
  double t_stat = dR * sqrt(numSubjects - 2) / sqrt(1 - dR * dR);

  boost::math::students_t dist(v);
  double                  q = boost::math::cdf(boost::math::complement(dist, fabs(t_stat) ) );

  dPDiff = 2 * q; // need to multiply by 2, because this is a two-sided test
  // just want to reject the possibility that there is no correlation

  dPGreater = boost::math::cdf(boost::math::complement(dist, t_stat) );

  dPSmaller = boost::math::cdf(dist, t_stat);

  return dR;

}

void output_vector(vnl_vector<double> data, std::string outbase, std::string toAppend)
{
  // Create and write out difference vector file:
  std::string outFile = outbase + toAppend;
  // std::cout << "outputFile: " << outFile << std::endl	;
  std::ofstream outputFile;

  outputFile.open(outFile.c_str(), std::ios::out);

  // Set up header:
  outputFile << "NUMBER_OF_POINTS = " << data.size() << std::endl;
  outputFile << "DIMENSION = 1" << std::endl;
  outputFile << "TYPE = Scalar" << std::endl;
  for( unsigned int feat = 0; feat < data.size(); ++feat )
    {
    outputFile << data(feat) << std::endl;
    }
  outputFile.close();
}

void output_matrix(vnl_matrix<double> data, std::string outbase, std::string toAppend)
{
  // Create and write out difference vector file:
  std::string outFile = outbase + toAppend;

  std::ofstream outputFile;

  outputFile.open(outFile.c_str(), std::ios::out);

  // Set up header:
  outputFile << "NUMBER_OF_POINTS = " << data.rows() << std::endl;
  outputFile << "DIMENSION = " << data.columns() << std::endl;
  outputFile << "TYPE = Vector" << std::endl;
  for( unsigned int feat = 0; feat < data.rows(); ++feat )
    {
    for( unsigned int dim = 0; dim < data.columns(); ++dim )
      {
      outputFile << data(feat, dim) << " ";
      }
    outputFile << std::endl;
    }
  outputFile.close();
}

void write_SubjectPoints( std::string fileName, unsigned int sub,
                          vnl_matrix<double>  * & featureValue,
                          unsigned int numFeatures,
                          MeshType::Pointer & surfaceMesh,
                          MeshSpatialObjectType::Pointer & SOMesh )
{
  // get the points
  PointsContainerPointer currentPoints;

  currentPoints = surfaceMesh->GetPoints();

  PointType currentPoint;
  for( unsigned int iI = 0; iI < numFeatures; ++iI )
    {
    for( unsigned int iJ = 0; iJ < dimension; iJ++ )
      {
      currentPoint[iJ] = (*featureValue)[sub][iI * dimension + iJ];
      }
    currentPoints->SetElement(iI, currentPoint);
    }

  // Write out total mean:
  surfaceMesh->SetPoints(currentPoints);
  SOMesh->SetMesh(surfaceMesh);
  MeshWriterType::Pointer writer = MeshWriterType::New();
  writer->SetInput(SOMesh);
  writer->SetFileName(fileName.c_str() );
  writer->Update();

}

void
load_MeshList_file( std::string filename, int surfaceColumn, unsigned int numIndependent,
                    unsigned int numGroupTypes, std::vector<int>  groupTypeColumns, std::vector<int> independentColumns,
                    bool scaleOn, int scaleColumn, unsigned int & numSubjects, unsigned int & numA, unsigned int & numB,
                    vnl_matrix<int> * & groupLabel, vnl_matrix<double> * & scaleFactor,
                    std::string * & meshFileNames,  vnl_matrix<double> * & indValue,
                    bool computeScaleFactorFromVolumes )
{

  cout << "filename: " << filename << " numIndependent: " << numIndependent << endl;
  // bool interactionTest =false; // Replace
  // unsigned int testColumn=1;
  const int   MAXLINE  = 5000;
  static char line[MAXLINE];

  char *        extension = const_cast<char *>(strrchr(filename.c_str(), '.') );
  std::ifstream datafile(filename.c_str(), std::ios::in);

  if( !datafile.is_open() )
    {
    std::cerr << "ERROR: Mesh list file does not exist" << filename << std::endl;
    exit(-1);
    }

  numSubjects = 0;

  datafile.clear();
  datafile.getline(line, MAXLINE);
  while( !datafile.eof() )
    {
    if( line[0] != '#' ) // skip over comments
      {
      numSubjects++;
      }
    datafile.getline(line, MAXLINE);

    }

  // if the input file is a csv file
  if( !strcmp(extension, ".csv") )
    {
    // because of the header
    numSubjects = numSubjects - 1;

    if( debug )
      {
      std::cout << "Num Subjects: " << numSubjects << std::endl;
      }

    scaleFactor =  new vnl_matrix<double>; scaleFactor->set_size(numSubjects, 1);
    meshFileNames = new std::string[numSubjects];

    // Create the data matrix:
    indValue = new vnl_matrix<double>;
    indValue->set_size(numSubjects, numIndependent); // Each subjects data is flattened, this data will go into X

    groupLabel = new vnl_matrix<int>;
    groupLabel->set_size(numSubjects, numGroupTypes);

    // read the list
    unsigned int curLine = 0;
    // char * cur_token;
    string Line;

    datafile.clear();
    datafile.seekg(ios::beg); // go back to the beginning of the
    // file to read the actual data

    getline(datafile, Line); // skip the csv file header

    int    temp_int = 0;
    double temp_double = 0;
    int    indexColumn = 0;

    while( !datafile.eof() )
      {
      getline(datafile, Line);

      if( Line[0] != '#' && !(Line.empty() ) )
        {
        char filename[MAXLINE];
        // int retval;

        istringstream buffer(Line);
        string        cur_token;

        while( getline(buffer, cur_token, ',') )
          {
          if( !cur_token.empty() )
            {
            for( unsigned int gr = 0; gr < numGroupTypes; ++gr ) // Read all of the group types
              {
              if( indexColumn == groupTypeColumns[gr] ) // if the curent column contains a group type
                {
                if( cur_token.empty() )
                  {
                  std::cerr << "ERROR: Could not read group types properly." << std::endl;
                  exit(-1);
                  }
                temp_int = atoi(cur_token.c_str() );
                (*groupLabel)[curLine][gr] = temp_int;
                }
              }
            if( scaleOn )
              {
              if( indexColumn == scaleColumn ) // if the current column contains a scale factor
                {
                if( cur_token.empty() )
                  {
                  std::cerr << "ERROR: Could not read scaling factor properly." << std::endl;
                  exit(-1);
                  }
                temp_double = atof(cur_token.c_str() );
                (*scaleFactor)[curLine][0] = temp_double;
                }
              }

            if( indexColumn == surfaceColumn ) // if the current column contains the surface file
              {
              if( cur_token.empty() )
                {
                std::cerr << "ERROR: Could not read filename properly." << std::endl;
                exit(-1);
                }
              strcpy(filename, cur_token.c_str() );
              meshFileNames[curLine] = std::string(filename);
              }
            for( unsigned int ind_vars = 0; ind_vars < numIndependent; ++ind_vars ) // Read all of the independent
                                                                                    // variables
              {

              if( indexColumn == independentColumns[ind_vars] ) // if the current column contain a independent variable
                {
                if( cur_token.empty() )
                  {
                  std::cerr << "ERROR: Could not read independent variable properly." << std::endl;
                  exit(-1);
                  }
                temp_double = atof(cur_token.c_str() );
                (*indValue)[curLine][ind_vars] = temp_double;
                }
              }
            indexColumn++;
            }

          else
            {
            std::cerr << "WARNING: Could not read group token. Skipping line" << std::endl;
            }

          }

        curLine++;
        indexColumn = 0;

        }

      }

    datafile.close();
    }

  else
    {
    if( debug )
      {
      std::cout << "Num Subjects: " << numSubjects << std::endl;
      }

    scaleFactor =  new vnl_matrix<double>; scaleFactor->set_size(numSubjects, 1);
    meshFileNames = new std::string[numSubjects];

    // Create the data matrix:
    indValue = new vnl_matrix<double>;
    indValue->set_size(numSubjects, numIndependent); // Each subjects data is flattened, this data will go into X

    groupLabel = new vnl_matrix<int>;
    groupLabel->set_size(numSubjects, numGroupTypes);

    // read the list
    unsigned int curLine = 0;
    char *       cur_token;

    datafile.clear();
    datafile.seekg(0, std::ios::beg); // go back to the beginning of the
    // file to read the actual data
    datafile.getline(line, MAXLINE);

    int    temp_int = 0;
    double temp_double = 0;

    while( !datafile.eof() )
      {
      if( line[0] != '#' && !(line[0] == 0) )
        {
        char filename[MAXLINE];
        int  retval;

        cur_token = strtok(line, " \t");

        if( cur_token != NULL )
          {
          for( unsigned int gr = 0; gr < numGroupTypes; ++gr ) // Read all of the group types
            {
            if( cur_token == NULL )
              {
              std::cerr << "ERROR: Could not read group types properly." << std::endl;
              exit(-1);
              }
            retval = sscanf(cur_token, "%d", &temp_int);
            (*groupLabel)[curLine][gr] = temp_int;
            cur_token = strtok(NULL, " \t");
            }

          if( cur_token == NULL )
            {
            std::cerr << "ERROR: Could not read scaling factor properly." << std::endl;
            exit(-1);
            }
          retval = sscanf(cur_token, "%lf", &temp_double); cur_token = strtok(NULL, " \t");
          (*scaleFactor)[curLine][0] = temp_double;

          if( cur_token == NULL )
            {
            std::cerr << "ERROR: Could not read filename properly." << std::endl;
            exit(-1);
            }
          retval = sscanf(cur_token, "%s", filename); cur_token = strtok(NULL, " \t");

          meshFileNames[curLine] = std::string(filename);
          for( unsigned int ind_vars = 0; ind_vars < numIndependent; ++ind_vars ) // Read all of the independent
                                                                                  // variables
            {
            if( cur_token == NULL )
              {
              std::cerr << "ERROR: Could not read independent variable properly." << std::endl;
              exit(-1);
              }
            retval = sscanf(cur_token, "%lf", &temp_double); cur_token = strtok(NULL, " \t");
            (*indValue)[curLine][ind_vars] = temp_double;
            }

          curLine++;

          }
        else
          {
          std::cerr << "WARNING: Could not read group token. Skipping line" << std::endl;
          }
        }
      datafile.getline(line, MAXLINE);

      }

    datafile.close();
    }

  if( computeScaleFactorFromVolumes )
    {
    std::cout << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "WARNING: Reinterpreting scaling column as volumes." << std::endl;
    std::cout << "... computing the scaling factor from the volumes." << std::endl;
    std::cout << "This is different from the old file format." << std::endl;
    std::cout << "MAKE SURE YOU UNDERSTAND WHAT YOU ENABLED WITH --computeScaleFactorFromVolumes" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    // reinterpret the scale column in the file as
    // volumes and compute the scale factor from them

    // first compute the average volume

    double dAverageVolume = 0;
    for( unsigned int iI = 0; iI < numSubjects; iI++ )
      {
      dAverageVolume += (*scaleFactor)[iI][0];
      }
    dAverageVolume /= numSubjects;

    std::cout << "Average volume is = " << dAverageVolume << std::endl;
    std::cout << "Computing scale factor by (volume/averageVolume)^(1/3)" << std::endl;

    double one_third = 1.0 / 3;
    for( unsigned int iI = 0; iI < numSubjects; iI++ )
      {
      (*scaleFactor)[iI][0] = pow( ( (*scaleFactor)[iI][0]) / dAverageVolume, one_third);
      }

    }

  std::cout << std::endl << std::endl;

  // possibly relabel the group associations
  // TODO: Check if this really should be done for the MANCOVA case
  // TODO: It may only be a historic remnant from the old
  // TODO: StatNonParam program

  int          preLabelA, preLabelB;
  unsigned int i;
  numA = 0, numB = 0;
  for( unsigned int group_type = 0; group_type < numGroupTypes; ++group_type )
    {
    // change labels to the predefined ones, this will fail if there are more than 2 labels in the file
    preLabelA = (*groupLabel)[0][group_type];
    preLabelB = (*groupLabel)[0][group_type];
    for( i = 0; i < numSubjects; i++ )
      {
      if( preLabelA != (*groupLabel)[i][group_type] )
        {
        if( preLabelB != preLabelA && preLabelB  != (*groupLabel)[i][group_type] )
          {
          std::cout << "Error: more than 2 labels in file" << std::endl;
          }
        else
          {
          preLabelB = (*groupLabel)[i][group_type];
          }
        }
      }
    for( i = 0; i < numSubjects; i++ )
      {
      if( preLabelA == (*groupLabel)[i][group_type] )
        {
        (*groupLabel)[i][group_type] = GROUP_A_LABEL;
        ++numA;
        }
      else if( preLabelB == (*groupLabel)[i][group_type] )
        {
        (*groupLabel)[i][group_type] = GROUP_B_LABEL;
        ++numB;
        }
      }
    if( debug )
      {
      std::cout << "data in group_type " << group_type << " has been relabeled: " <<  preLabelA << " --> group A = "
                << GROUP_A_LABEL
                << " ; " << preLabelB << " --> group B = " << GROUP_B_LABEL << std::endl;
      std::cout << "#(A)= " << numA << "; #(B)= " << numB << std::endl;
      }
    } // END of For-loop

  // TODO Remove leading and ending " from meshFileNames if present

  if( meshFileNames[1].compare(1, 1, "/") == 0 )
    {
    std::vector<string> meshFileNamesCopy;
    for( unsigned int i = 0; i < numSubjects; i++ )
      {
      meshFileNamesCopy.push_back(".");
      }

    if( debug )
      {
      for( unsigned int i = 0; i < numSubjects; i++ )
        {
        meshFileNamesCopy.at(i).assign(meshFileNames[i], 1, meshFileNames[i].size() - 2);
        std::cout << "reading Mesh.Copy " << meshFileNamesCopy.at(i) << std::endl;
        }
      }
    for( unsigned int i = 0; i < numSubjects; i++ )
      {
      meshFileNames[i] = meshFileNamesCopy[i];
      }
    }
  // debug info
  if( debug )
    {
    for( unsigned int i = 0; i < numSubjects; i++ )
      {
      std::cout << meshFileNames[i] << " " << (*groupLabel)[i][0] << " " << (*scaleFactor)[i][0] << std::endl;
      }
    }

}

void load_Meshes( bool scaleOn,  unsigned int numSubjects,
                  unsigned int numIndependent,
                  vnl_matrix<double> * & scaleFactor,  vnl_matrix<double> * & indValue,
                  std::string * & meshFileNames, unsigned int & numFeatures, vnl_matrix<double>  * & featureValue,
                  MeshType::Pointer & surfaceMesh,
                  MeshSpatialObjectType::Pointer & SOMesh )
{
  // Read the meshes
  featureValue = new vnl_matrix<double>;

  MeshReaderType::Pointer reader = MeshReaderType::New();
  PointsContainerPointer  points;
  for( unsigned int index = 0; index < numSubjects; index++ )
    {

    // if the surface file is a vtk file
    char *extension = const_cast<char *>(strrchr(meshFileNames[index].c_str(), '.') );
    std::cout << "reading Mesh " << meshFileNames[index].c_str() << " - " << extension << std::endl;

    if( !strcmp(extension, ".vtk") ) // TODO
      {
      // read vtk file
      vtkPolyDataReader *mesh = vtkPolyDataReader::New();
      mesh->SetFileName(meshFileNames[index].c_str() );

      // std::cout << "converting vtk to itk" << std::endl;

      try
        {
        mesh->Update();
        }
      catch( itk::ExceptionObject ex )
        {
        std::cout << "Error reading VTK meshfile:  " << meshFileNames[index] << std::endl << "ITK error: "
                  << ex.GetDescription() << std::endl;
        exit(-3);
        }

      vtkPolyData *polyData = mesh->GetOutput();

      // convert polydata to itk mesh
      int length = strlen(meshFileNames[index].c_str() ) - 4;

      char *meshFile = new char[length + 6];
      strcpy(meshFile, meshFileNames[index].c_str() );
      meshFile[length] = '\0';
      strcat(meshFile, ".meta");

      vtkPolyDataToitkMesh *vtkItkConverter = new vtkPolyDataToitkMesh();
      vtkItkConverter->SetInput( polyData );

      // write out the itk meta mesh file
      MeshConverterType::Pointer itkConverter = MeshConverterType::New();
      itkMeshSOType::Pointer meshSO = itkMeshSOType::New();
      meshSO->SetMesh( vtkItkConverter->GetOutput() );
      itkConverter->WriteMeta( meshSO, meshFile );

      mesh->Delete();

      // read mesh file
      try
        {
        reader->SetFileName(meshFile);
        reader->Update();
        }
      catch( itk::ExceptionObject ex )
        {
        std::cout << "Error reading META meshfile:  " << meshFileNames[index] << std::endl << "ITK error: "
                  << ex.GetDescription() << std::endl;
        }

      }

    else
      {
      try
        {
        reader->SetFileName(meshFileNames[index].c_str() );
        reader->Update();
        }
      catch( itk::ExceptionObject ex )
        {
        std::cout << "Error reading META meshfile:  " << meshFileNames[index] << std::endl << "ITK error: "
                  << ex.GetDescription() << std::endl;
        }
      }

    MeshReaderType::SceneType::Pointer          scene = reader->GetScene();
    MeshReaderType::SceneType::ObjectListType * objList =  scene->GetObjects(1, NULL);

    MeshReaderType::SceneType::ObjectListType::iterator it = objList->begin();

    itk::SpatialObject<3> * curObj = *it;
    SOMesh = dynamic_cast<MeshSpatialObjectType *>(curObj);
    surfaceMesh = SOMesh->GetMesh();
    points = surfaceMesh->GetPoints();

    if( index == 0 )
      {
      numFeatures = points->Size();
      featureValue->set_size(numSubjects, numFeatures * 3 + numIndependent); // Each subjects data is flattened, this
                                                                             // data will go into X
      }

    float xCOG = 0; float yCOG = 0; float zCOG = 0;
    for( unsigned int pointID = 0; pointID < numFeatures; pointID++ )
      {
      PointType curPoint =  points->GetElement(pointID);
      xCOG = xCOG + curPoint[0];
      yCOG = yCOG + curPoint[1];
      zCOG = zCOG + curPoint[2];
      }
    xCOG = xCOG / numFeatures;
    yCOG = yCOG / numFeatures;
    zCOG = zCOG / numFeatures;
    float COG[3]; COG[0] = xCOG; COG[1] = yCOG; COG[2] = zCOG;
    for( unsigned int pointID = 0; pointID < numFeatures; pointID++ )
      {
      PointType curPoint =  points->GetElement(pointID);
      for( unsigned int dim = 0; dim < 3; dim++ )
        {
        if( scaleOn )
          {
          (*featureValue)[index][pointID * 3 + dim] = (curPoint[dim] - COG[dim]) / (*scaleFactor)[index][0] + COG[dim];
          }
        else
          {
          (*featureValue)[index][pointID * 3 + dim] = curPoint[dim];
          }
        }
      }
    // Store Independent variables after tupel data:
    for( unsigned int ind_vars = 0; ind_vars < numIndependent; ++ind_vars )
      {
      featureValue->set_column(numFeatures * 3 + ind_vars, indValue->get_column(ind_vars) );
      }
    }

}

//
// Function added by bp2009
// Added to include the longitudinal analysis in the study...
//
void load_KWMreadableInputFile( bool scaleOn,  unsigned int numSubjects,
                                unsigned int numIndependent,
                                vnl_matrix<double> * & scaleFactor,  vnl_matrix<double> * & indValue,
                                std::string * & meshFileNames, unsigned int & numFeatures,
                                vnl_matrix<double>  * & featureValue, MeshType::Pointer & surfaceMesh,
                                MeshSpatialObjectType::Pointer & SOMesh )
{
  // Read the KWMfile
  // By now, and only for testing purposes, I am going to ignore a lot of the input parameters to the function... I
  // copied directly from load_Meshes, and there are a lot of things I dont need

  featureValue = new vnl_matrix<double>;
  for( unsigned int index = 0; index < numSubjects; index++ )
    {
    try
      {
      // std::cout<< "Working on file...  "<< meshFileNames[index] << std::endl << std::endl;

      char               line[70];
      ifstream           inputFile;
      char *             aux;
      vnl_vector<double> curPoint(3, 0.0);
      int                NPoints, pointID, pointCont;

      inputFile.open(meshFileNames[index].c_str(), std::ios::in);
      inputFile.getline(line, 70, '\n');
      aux = strtok(line, " = ");
      aux = strtok(NULL, " = ");
      NPoints = atoi(aux);
      inputFile.getline(line, 70, '\n');
      inputFile.getline(line, 70, '\n');

      pointID = 0; pointCont = 0;
      // Start reading the Input File
      while( !inputFile.getline(line, 70, '\n').eof() )
        {
        aux = strtok(line, " ");
        while( aux != NULL )
          {
          curPoint[pointCont] = atof(aux);
          // std::cout<< aux << std::endl;
          aux = strtok(NULL, " ");
          pointCont++;
          }

        pointCont = 0;
        // std::cout<< "" "Point " << pointID << " values " << curPoint[0] << " " << curPoint[1] << " " << curPoint[2]
        // << std::endl;

        if( index == 0 )
          {
          numFeatures = NPoints;
          featureValue->set_size(numSubjects, numFeatures * 3 + numIndependent); // Each subjects data is flattened,
                                                                                 // this data will go into X
          }
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          if( scaleOn )
            {
            // std::cout<< "Aqui..." << std::endl;
            (*featureValue)[index][pointID * 3 + dim] = curPoint[dim] / (*scaleFactor)[index][0];
            }
          else
            {
            // std::cout<< "Aqui no scale..." << std::endl;
            (*featureValue)[index][pointID * 3 + dim] = curPoint[dim];
            }
          }

        // std::cout<< "File id " << index << " Point " << pointID << " values " << curPoint[0] << " " << curPoint[1] <<
        // " " << curPoint[2] << std::endl;
        pointID++;
        }

      inputFile.close();
      // End reading the Input File

      } // end FOR index
    catch( itk::ExceptionObject ex )
      {
      std::cout << "Error reading KWMfile:  " << meshFileNames[index] << std::endl << "ITK error: "
                << ex.GetDescription() << std::endl;
      exit(-3);
      }
    // std::cout<< "Finished reading file id " << index << std::endl;
    // Store Independent variables after tupel data:
    for( unsigned int ind_vars = 0; ind_vars < numIndependent; ++ind_vars )
      {
      featureValue->set_column(numFeatures * 3 + ind_vars, indValue->get_column(ind_vars) );
      }
    }

  std::cout << "KWM Files Readed..." << std::endl;

}

//
// Function modified by bp2009
// Included new input argument KWMreadableInputFile...
// Controlling unallowed operations wrt to Mesh Operations (no Meshes appear in a featural analysis)
//
void write_SurfaceProperties(  std::string outbase,  PointsContainerPointer & meanPoints,
                               PointsContainerPointer & meanPointsA,  PointsContainerPointer & meanPointsB,
                               vnl_matrix<double>& diffVectors,
                               vnl_vector<double>& normProjections, vnl_vector<double>& DiffMagnitude,
                               vnl_matrix<double>& zScoresProjected,
                               vnl_matrix<double>& zScores, bool writeOutZScores, MeshType::Pointer & surfaceMesh,
                               MeshSpatialObjectType::Pointer & SOMesh, std::string * & meshFileNames,
                               int KWMreadableInputFile)
{

  // std::cout<< "Starting to write surface properties..." << std::endl;

  if( KWMreadableInputFile == 0 )
    {

    // if the input file is a csv file
//  char *extension=strchr(outbase.c_str(),'.');
//  if(!strcmp(extension,".csv"))
//  {
//    int index=outbase.find(".csv",0);
//    outbase=outbase.erase(index,outbase.size());
//  }
// Write out total mean:
    surfaceMesh->SetPoints(meanPoints);
    SOMesh->SetMesh(surfaceMesh);
    MeshWriterType::Pointer writer = MeshWriterType::New();
    writer->SetInput(SOMesh);
    std::string FilenameAverage(outbase);

    FilenameAverage = FilenameAverage + std::string("_meanAll_uncorrected.meta");
    writer->SetFileName(FilenameAverage.c_str() );
    writer->Update();

    // Write out Group A mean:
    surfaceMesh->SetPoints(meanPointsA);
    SOMesh->SetMesh(surfaceMesh);
    writer->SetInput(SOMesh);
    std::string FilenameAverageA(outbase);
    FilenameAverageA = FilenameAverageA + std::string("_meanA.meta");
    writer->SetFileName(FilenameAverageA.c_str() );
    writer->Update();

    // Write out group B mean:
    surfaceMesh->SetPoints(meanPointsB);
    SOMesh->SetMesh(surfaceMesh);
    writer->SetInput(SOMesh);
    std::string FilenameAverageB(outbase);
    FilenameAverageB = FilenameAverageB + std::string("_meanB.meta");
    writer->SetFileName(FilenameAverageB.c_str() );
    writer->Update();
    } // end if (KWMreadableInputFile == 0)

  // Create and write out difference vector file:

  output_matrix(diffVectors, outbase, std::string("_diffMesh.txt") );

  // normal projections

  output_vector(normProjections, outbase, std::string("_normProjections.txt") );
  output_vector(DiffMagnitude, outbase, std::string("_DiffMagnitude.txt") );

  // write out z-scores

  if( writeOutZScores )
    {

    unsigned int numSubjects = zScores.columns();
    for( unsigned sub = 0; sub < numSubjects; sub++ )
      {
      boost::filesystem::path my_path( meshFileNames[sub] );

      std::ostringstream fileNameSuffixP;
      fileNameSuffixP << "-" << my_path.stem() << "-zScore-projected" << ".txt";
      output_vector(zScoresProjected.get_column(sub), outbase, fileNameSuffixP.str() );
      std::ostringstream fileNameSuffix;
      fileNameSuffix << "-" << my_path.stem() << "-zScore-mahalanobis" << ".txt";
      output_vector(zScores.get_column(sub), outbase, fileNameSuffix.str() );

      }

    }

}

//
// Function modified by bp2009
// Included new input argument KWMreadableInputFile...
// Controlling unallowed operations wrt to Mesh Operations (no Meshes appear in a featural analysis)
//
void compute_SurfaceProperties( PointsContainerPointer & meanPoints, PointsContainerPointer & meanPointsA,
                                PointsContainerPointer & meanPointsB, vnl_matrix<double> & diffVectors,
                                vnl_vector<double>& normProjections,
                                vnl_vector<double>& DiffMagnitude, vnl_matrix<double>& zScores,
                                vnl_matrix<double>& zScoresProjected,
                                unsigned int numSubjects, unsigned int numA, unsigned int numB,
                                unsigned int numFeatures, vnl_matrix<int> * & groupLabel,
                                vnl_matrix<double> * & featureValue, vtkPolyDataNormals *& MeshNormals,
                                MeshType::Pointer & surfaceMesh,
                                std::string outbase, std::string * & meshFileNames,
                                int KWMreadableInputFile)
{
  double       tmp_pnt[dimension], tmpA[dimension], tmpB[dimension];
  PointType    A_pnt, B_pnt;
  unsigned int dim;

  // need to make sure that the overall mean
  // is computed correctly
  // let's do this by averaging the mean of A and B
  // rather than by computing the overall mean
  // which would bias towards the group with more subjects
  for( unsigned int pointID = 0; pointID < numFeatures; pointID++ )
    {
    for( dim = 0; dim < dimension; dim++ )
      {
      tmpA[dim] = 0; tmpB[dim] = 0;
      }
    for( unsigned int sub = 0; sub < numSubjects; ++sub )
      {
      if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
        {
        for( dim = 0; dim < dimension; ++dim )
          {
          tmpA[dim] += (*featureValue)[sub][pointID * dimension + dim];
          }
        }
      else
        {
        for( dim = 0; dim < dimension; ++dim )
          {
          tmpB[dim] += (*featureValue)[sub][pointID * dimension + dim];
          }
        }
      }

    // normalize to compute the mean shapes
    // for the individual groups as well as
    // an overall mean shape

    if( numA > 0 ) // only if there are elements in group A
      {
      for( dim = 0; dim < dimension; dim++ ) // Normalize
        {
        tmpA[dim] /= numA;
        }
      }

    if( numB > 0 ) // only if there are elements in group B
      {
      for( dim = 0; dim < dimension; dim++ ) // Normalize
        {
        tmpB[dim] /= numB;
        }
      }

    if( numA > 0 && numB > 0 )
      {
      for( dim = 0; dim < dimension; dim++ ) // Normalize
        {
        tmp_pnt[dim] = 0.5 * (tmpA[dim] + tmpB[dim]);
        }
      }
    else if( numA > 0 )
      {
      for( dim = 0; dim < dimension; dim++ ) // Normalize
        {
        tmp_pnt[dim] = tmpA[dim];
        }
      }
    else if( numB > 0 )
      {
      for( dim = 0; dim < dimension; dim++ ) // Normalize
        {
        tmp_pnt[dim] = tmpB[dim];
        }
      }
    else
      {
      for( dim = 0; dim < dimension; dim++ ) // Normalize
        {
        tmp_pnt[dim] = 0;
        }
      }

    meanPoints->InsertElement(pointID, PointType(tmp_pnt) );
    meanPointsA->InsertElement(pointID, PointType(tmpA) );
    meanPointsB->InsertElement(pointID, PointType(tmpB) );

    if( numA > 0 && numB > 0 )
      {
      for( dim = 0; dim < dimension; ++dim ) // Initialize
        {
        diffVectors[pointID][dim] = tmpA[dim] - tmpB[dim];
        }
      }
    else // one of the groups is empty, so we cannot compute a
    // difference and therefore set everything to zero
      {
      for( dim = 0; dim < dimension; ++dim ) // Initialize
        {
        diffVectors[pointID][dim] = 0;
        }
      }

    }

  // compute the mean surface normals
  // as well as the signed distances by projecting the surface
  // difference onto the surface normals

  surfaceMesh->SetPoints(meanPoints);

  itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
  convertMeshToVTK->SetInput(surfaceMesh);
  vtkPolyData * vtkMesh = convertMeshToVTK->GetOutput();

  MeshNormals->SetComputePointNormals(1);
  MeshNormals->SetComputeCellNormals(0);
  MeshNormals->SetSplitting(0);
  MeshNormals->AutoOrientNormalsOn();
  MeshNormals->ConsistencyOn();
  MeshNormals->FlipNormalsOff();    // all normals are outward pointing
  MeshNormals->SetInput(vtkMesh);
  MeshNormals->Update();

  vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
  vtkMeshNormals->Update();

  vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
  vtkDataArray * ArrayNormal = NormalPoints->GetNormals();

  // now compute the signed distances

  double *tempNorm;
  for( unsigned int feat = 0; feat < numFeatures; ++feat )
    {
    if( KWMreadableInputFile == 0 )
      {
      tempNorm = ArrayNormal->GetTuple(feat);  // here is crashing... why

      // all of this should automatically become zero if there is only
      // one group, because in this case diffVectors will be all zero

      normProjections(feat) = diffVectors(feat, 0) * tempNorm[0] + diffVectors(feat, 1) * tempNorm[1] + diffVectors(
          feat, 2) * tempNorm[2];
      } // end if (KWMreadableInputFile == 0)
    DiffMagnitude(feat) = sqrt(diffVectors(feat, 0) * diffVectors(feat, 0) + diffVectors(feat, 1) * diffVectors(feat,
                                                                                                                1)
                               + diffVectors(feat, 2) * diffVectors(feat, 2) );
    if( normProjections(feat) < 0 ) // Negative, means pointing in
      {
      DiffMagnitude(feat) = -DiffMagnitude(feat);
      }
    }

  // Computing the z-scores ...

  // first the projected case
  // let's compute the mean projected distance of A wrt. B and the
  // other way around with respect to the mean surface

  if( KWMreadableInputFile == 0 )
    {
    for( unsigned int feat = 0; feat < numFeatures; feat++ )
      {
      double *           tempNorm = ArrayNormal->GetTuple(feat);
      PointType          currentMeanPoint =  meanPoints->GetElement(feat);
      vnl_vector<double> dProj(numSubjects, 0);
      double             mu_aP = 0;
      double             mu_bP = 0;
      for( unsigned int sub = 0; sub < numSubjects; sub++ )
        {
        for( unsigned int tup = 0; tup < dimension; ++tup )
          {
          // compute the projection along the mean normal
          dProj(sub) += ( (*featureValue)[sub][feat * dimension + tup] - currentMeanPoint[tup]) * tempNorm[tup];
          }
        // now switch depending on type

        if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
          {
          mu_aP += dProj(sub);
          }
        else
          {
          mu_bP += dProj(sub);
          }
        }
      // now compute the mean
      if( numA > 0 )
        {
        mu_aP /= numA;
        }
      if( numB > 0 )
        {
        mu_bP /= numB;
        }

      // now compute the variance

      double sigma_aP = 0;
      double sigma_bP = 0;
      for( unsigned int sub = 0; sub < numSubjects; sub++ )
        {
        // now switch depending on type

        if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
          {
          sigma_aP += (dProj(sub) - mu_aP) * (dProj(sub) - mu_aP);
          }
        else
          {
          sigma_bP += (dProj(sub) - mu_bP) * (dProj(sub) - mu_bP);
          }
        }
      // now compute the standard deviation
      if( numA > 0 )
        {
        sigma_aP /= numA;
        sigma_aP = sqrt(sigma_aP);
        }
      if( numB > 0 )
        {
        sigma_bP /= numB;
        sigma_bP = sqrt(sigma_bP);
        }
      // and finally the z-scores
      for( unsigned int sub = 0; sub < numSubjects; sub++ )
        {

        // now switch depending on type

        if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
          {
          zScoresProjected(feat, sub) = (dProj(sub) - mu_bP) / sigma_bP;
          }
        else
          {
          zScoresProjected(feat, sub) = (dProj(sub) - mu_aP) / sigma_aP;
          }
        }
      }

    } // end of if KWMreadableInputFile
  // now do all of the z-score business using the Mahalanobis distance
  for( unsigned int feat = 0; feat < numFeatures; feat++ )
    {
    vnl_vector<double> mu_b(dimension, 0);
    vnl_vector<double> mu_a(dimension, 0);

    vnl_vector<double> currentVec(dimension, 0);
    for( unsigned int sub = 0; sub < numSubjects; sub++ )
      {
      for( unsigned int tup = 0; tup < dimension; ++tup )
        {
        currentVec(tup) = (*featureValue)[sub][feat * dimension + tup];
        }
      // now switch depending on type

      if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
        {
        mu_a += currentVec;
        }
      else
        {
        mu_b += currentVec;
        }
      }
    // now compute the mean
    if( numA > 0 )
      {
      mu_a /= numA;
      }
    if( numB > 0 )
      {
      mu_b /= numB;
      }

    // now compute the variances

    vnl_matrix<double> var_b(dimension, dimension, 0);
    vnl_matrix<double> var_a(dimension, dimension, 0);
    for( unsigned int sub = 0; sub < numSubjects; sub++ )
      {
      for( unsigned int tup = 0; tup < dimension; ++tup )
        {
        currentVec(tup) = (*featureValue)[sub][feat * dimension + tup];
        }

      // now switch depending on type

      if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
        {
        var_a += outer_product<double>(currentVec - mu_a, currentVec - mu_a);
        }
      else
        {
        var_b += outer_product<double>(currentVec - mu_b, currentVec - mu_b);
        }
      }
    // now compute the standard deviation
    if( numA > 0 )
      {
      var_a /= numA;
      }
    if( numB > 0 )
      {
      var_b /= numB;
      }

    // and finally compute the z-scores by computing
    // the Mahalanobis distance

    // compute the inverse of the variance

    vnl_matrix<double> var_bInv(dimension, dimension);
    vnl_matrix<double> var_aInv(dimension, dimension);

    var_aInv = vnl_matrix_inverse<double>(var_a);
    var_bInv = vnl_matrix_inverse<double>(var_b);
    for( unsigned int sub = 0; sub < numSubjects; sub++ )
      {
      for( unsigned int tup = 0; tup < dimension; ++tup )
        {
        currentVec(tup) = (*featureValue)[sub][feat * dimension + tup];
        }

      // now switch depending on type

      if( (*groupLabel)[sub][0] == GROUP_A_LABEL )
        {
        zScores(feat, sub) = sqrt( dot_product(currentVec - mu_b, var_bInv * (currentVec - mu_b) ) );
        }
      else
        {
        zScores(feat, sub) = sqrt( dot_product(currentVec - mu_a, var_aInv * (currentVec - mu_a) ) );
        }
      }

    }
}

vnl_vector<double> fdrCorrection( vnl_vector<double>& rawP, double fdrLevel, double & fdrThresh )
{
  // to determine the adjusted level, we first need to sort the raw p
  // values

  vnl_vector<double> pSorted( rawP );
  std::sort( pSorted.begin(), pSorted.end() );

  fdrThresh = 0;

  unsigned int numFeatures = rawP.size();
  for( int pointID = (int)(numFeatures - 1); pointID >= 0; --pointID )
    {
    if( pSorted(pointID) <= (pointID + 1) * fdrLevel / numFeatures )
      {
      fdrThresh = pSorted(pointID);
      break;
      }
    }

  vnl_vector<double> fdrP( rawP );
  double             fdrFactor = 0;
  if( fdrThresh > 0 )
    {
    fdrFactor = fdrLevel / fdrThresh;
    }
  if( fdrFactor < 1.0 )
    {
    fdrFactor = 1.0;
    }
  for( unsigned int pointID = 0; pointID < numFeatures; ++pointID )
    {
    if( fdrThresh > 0 )
      {
      fdrP(pointID) = rawP[pointID] * fdrFactor;
      }
    else
      {
      fdrP(pointID) = 1.0;  // there is nothing signicant (sorry), so let's set the p-value to 1
      }
    if( fdrP(pointID) > 1.0 )
      {
      fdrP(pointID) = 1.0; // bound the p-value at 1.0, p-values over 1.0 do not make sense
      }
    }

  return fdrP;

}

vnl_vector<double> bonferroniCorrection( vnl_vector<double>& rawP )
{
  vnl_vector<double> bP( rawP );
  unsigned int       sz = rawP.size();
  for( unsigned int iI = 0; iI < sz; iI++ )
    {
    bP(iI) = rawP(iI) * sz;
    if( bP(iI) > 1 )
      {
      bP(iI) = 1;
      }
    }
  return bP;
}

void Meta2VTK(char * infile, char* outfile)
{
  typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
  typedef itk::Mesh<double, 3, MeshTraitsType>                        itkMeshType;
  typedef itk::MeshSpatialObject<itkMeshType>                         itkMeshSOType;
  typedef itk::MetaMeshConverter<3, double, MeshTraitsType>           MeshConverterType;

  // read the data in meta format
  MeshConverterType::Pointer itkConverter = MeshConverterType::New();
  itkMeshSOType::Pointer meshSO =
    dynamic_cast<itkMeshSOType *>(itkConverter->ReadMeta(infile).GetPointer() );
  itkMeshType::Pointer   mesh = meshSO->GetMesh();
  // delete (itkConverter);

  // convert to vtk format
  itkMeshTovtkPolyData * ITKVTKConverter = new itkMeshTovtkPolyData;
  ITKVTKConverter->SetInput( mesh );

  // write out the vtk mesh
  vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
  writer->SetInput( ITKVTKConverter->GetOutput() );
  writer->SetFileName( outfile );
  writer->Update();

  writer->Delete();
  // delete (ITKVTKConverter);

}

void write_ColorMap(std::string outbase, bool interactionTest, double significanceLevel, double FDRdiscoveryLevel)
{
  int end;

  if( interactionTest )
    {
    end = 8;
    }                       // nb of file have to be created
  else
    {
    end = 3;
    }          // nb of file have to be created
  char InputMetaFile[512];
  char InputVTKFile[512];
  char TextFile[512];
  char OutputFile[512];
  strcpy(OutputFile, outbase.c_str() );
  strcat(OutputFile, "_meanAll_uncorrected.vtk");
  strcpy(InputMetaFile, outbase.c_str() );
  strcpy(InputVTKFile, outbase.c_str() );
  strcat(InputVTKFile, "_meanAll_uncorrected.vtk");
  strcat(InputMetaFile, "_meanAll_uncorrected.meta");
  Meta2VTK(InputMetaFile, InputVTKFile);
  for( int i = 0; i < end; i++ )
    {

    std::vector<const char *> args;
    char*                     data = NULL;
    int                       length;
    double                    timeout = 0.05;
    int                       result;
    std::ostringstream        tmp1;
    tmp1 << significanceLevel;
    std::string        RawPvalueColorMapString(tmp1.str() );
    std::ostringstream tmp2;
    tmp2 << FDRdiscoveryLevel;
    std::string FDRPvalueColorMapString(tmp2.str() );

    switch( i )
      {
      case 0:

        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_mancovaRawP.txt");
        args.push_back("MeshMath");
        args.push_back(OutputFile);
        args.push_back(OutputFile);
        args.push_back("-significanceLevel");
        args.push_back(RawPvalueColorMapString.c_str() );
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("RawP");
        args.push_back(0);
        break;

      case 1:

        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_mancovaFDRP.txt");
        args.push_back("MeshMath");
        args.push_back(OutputFile);
        args.push_back(OutputFile);
        args.push_back("-significanceLevel");
        args.push_back(FDRPvalueColorMapString.c_str() );
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("FDRP");
        args.push_back(0);
        break;

      case 2:
        if( interactionTest )
          {
          strcpy(TextFile, outbase.c_str() );
          strcat(TextFile, "_normProjectionsPearsonPval.txt");
          args.push_back("MeshMath");
          args.push_back(OutputFile);
          args.push_back(OutputFile);
          args.push_back("-significanceLevel");
          args.push_back(RawPvalueColorMapString.c_str() );
          args.push_back("-KWMtoPolyData");
          args.push_back(TextFile);
          args.push_back("normProjectionsPearsonPval");
          }

        else
          {
          strcpy(TextFile, outbase.c_str() );
          strcat(TextFile, "_DiffMagnitude.txt");
          args.push_back("MeshMath");
          args.push_back(OutputFile);
          args.push_back(OutputFile);
          args.push_back("-KWMtoPolyData");
          args.push_back(TextFile);
          args.push_back("DiffMagnitude");
          }
        args.push_back(0);

        break;

      case 3:

        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_normProjectionsSpearmanPvalFDR.txt");
        args.push_back("MeshMath");
        std::cout << OutputFile << std::endl;
        args.push_back(OutputFile);
        args.push_back(OutputFile);
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("normProjectionsSpearmanPvalFDR");
        args.push_back("-significanceLevel");
        args.push_back(FDRPvalueColorMapString.c_str() );
        args.push_back(0);

        break;

      // case 6:
      case 4:

        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_normProjectionsPearsonPvalFDR.txt");
        args.push_back("MeshMath");
        args.push_back(OutputFile);
        args.push_back(OutputFile);
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("normProjectionsPearsonPvalFDR");
        args.push_back("-significanceLevel");
        args.push_back(FDRPvalueColorMapString.c_str() );
        args.push_back(0);

        break;

      // case 7:
      case 5:
        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_normProjectionsPearson.txt");
        args.push_back("MeshMath");
        args.push_back(OutputFile);
        args.push_back(OutputFile);
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("normProjectionsPearson");
        args.push_back(0);

        break;

      // case 8:
      case 6:
        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_normProjectionsSpearman.txt"); // TODO fichier txt avec les info pour les overlay
                                                          // (activescalar)
        args.push_back("MeshMath");                       // PROG appele
        args.push_back(OutputFile);                       // origine
        args.push_back(OutputFile);                       // sortie (ou serotn stoquees les info de l overlay
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("normProjectionsSpearman");   // non de ton active scalar
        args.push_back(0);

        break;

      // case 9:
      case 7:

        strcpy(TextFile, outbase.c_str() );
        strcat(TextFile, "_normProjectionsSpearmanPval.txt");
        args.push_back("MeshMath");
        args.push_back(OutputFile);
        args.push_back(OutputFile);
        args.push_back("-KWMtoPolyData");
        args.push_back(TextFile);
        args.push_back("normProjectionsSpearmanPval");
        args.push_back("-significanceLevel");
        args.push_back(RawPvalueColorMapString.c_str() );
        args.push_back(0);

        break;

      }

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
        for( int i = 0; i < length; i++ )
          {
          st << data[i];
          }
        string dim = st.str();
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

int minmax(const std::string & filename, double *min, double* max)
{
  std::ifstream fichier(filename.c_str(), std::ios::in);
  std::string   s;

  if( fichier )
    {
    *min = 0;
    *max = 0;
    double i;
    while( std::getline(fichier, s) )
      {
      const char *x = s.c_str();
      i = strtod(x, NULL);

      if( i < *min )
        {
        *min = i;
        }
      else if( i > *max )
        {
        *max = i;
        }
      }

    }
  else
    {
    std::cout << "can not open the file " << filename << std::endl;
    }

  fichier.close();
  return 0;
}

void write_MRMLScene(std::string outbase, bool interactionTest)
{

  string m_directionToDisplay;

  std::vector<const char *> args;

  char*  data = NULL;
  int    length;
  double timeout = 0.05;
  int    result;

  char nameVTK[512];
  std::strcpy(nameVTK, outbase.c_str() );
  std::strcat(nameVTK, "_meanAll_uncorrected.vtk");
  stringstream ss;
  string       s_nameVTK;
  ss << nameVTK;
  ss >> s_nameVTK;
  size_t found;
  found = s_nameVTK.find_last_of("/\\");

  // create a directory to save the transfoms files
  std::string s_outputdirectory;
  s_outputdirectory.append(s_nameVTK);
  s_outputdirectory.erase(found, s_outputdirectory.size() - 1);
  s_outputdirectory.append("/transformFiles");
  itksys::SystemTools::MakeDirectory(s_outputdirectory.c_str() );

  // the name of the vtk
  s_nameVTK.erase(0, found + 1);
  std::strcpy(nameVTK, s_nameVTK.c_str() );

  vector<double> Dim;
  char           nameVTK2[512];
  std::strcpy(nameVTK2, outbase.c_str() );
  std::strcat(nameVTK2, "_meanAll_uncorrected.vtk");
  Dim = SetImageDimensions( (char *)(nameVTK2) );

  std::string              pos;
  std::vector<std::string> pos_all;
  std::string              fidupos;
  std::vector<std::string> fidupos_all;
  std::string              meandiff;

  // remove transforms
  s_outputdirectory.append("/");
  bool transformDirectoryEmpty = DirectoryIsEmpty(s_outputdirectory.c_str() );
  if( !transformDirectoryEmpty )
    {
    int length = strlen(s_outputdirectory.c_str() );
    length = length + 10;

    DIR *          pdir = NULL;
    struct dirent *pent;
    pdir = opendir(s_outputdirectory.c_str() );
    while( (pent = readdir(pdir) ) )
      {
      char *file = NULL;
      file = new char[512];
      strcpy(file, s_outputdirectory.c_str() );
      strcat(file, pent->d_name);
      if( file[length] != '.' )
        {
        remove(file);
        }
      }
    }

  if( interactionTest )
    {
    char mrmlfile[512];
    std::strcpy(mrmlfile, outbase.c_str() );
    std::strcat(mrmlfile, "_InteractionTest_MRMLscene.mrml");

    args.push_back("CreateMRML");
    args.push_back(mrmlfile);

    double x, y, z, x2, fidux, fiduy, fiduz, fidux2;
    x = -165;
    y = 27;
    z = -62;
    /*fidux= x +Dim[0]/2+5;
    fiduy= y+Dim[1]/2;
    fiduz=90;*/
    fidux = x + Dim[4] - 10;
    fiduy = y + Dim[5];
    fiduz = (Dim[7] + Dim[8]) / 2 - z;

    // shape and tranfoms

// RawP
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransRawP.tfm"); args.push_back("-n");
    args.push_back("transRawP"); args.push_back("-l");

    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back("RawP");
    args.push_back("-p"); args.push_back("transRawP"); args.push_back("-as"); args.push_back("RawP"); args.push_back(
      "-cc"); args.push_back("customLUT_RawP.txt");

    fidupos.append(Convert_Double_To_CharArray(fidux) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("MANCOVA_RawP"); args.push_back("-lbl"); args.push_back(
      "MANCOVA_RawP"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

    pos.clear();
    fidupos.clear();

// FDRP
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransFDRP.tfm"); args.push_back("-n");
    args.push_back("transFDRP"); args.push_back("-l");

    z = z - Dim[2] - 10;
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back("FDRP");
    args.push_back("-p"); args.push_back("transFDRP"); args.push_back("-as"); args.push_back("FDRP"); args.push_back(
      "-cc"); args.push_back("customLUT_FDRP.txt");

    fiduz = fiduz + Dim[2] + 5;
    fidupos.append(Convert_Double_To_CharArray(fidux) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("MANCOVA_FDRP"); args.push_back("-lbl"); args.push_back(
      "MANCOVA_FDRP"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

    fidupos.clear();
    pos.clear();

// normProjectionsPearsonPvalFDR
    args.push_back("-t"); args.push_back("-f"); args.push_back(
      "./transformFiles/TransnormProjectionsPearsonPvalFDR.tfm"); args.push_back("-n"); args.push_back(
      "transnormProjectionsPearsonPvalFDR"); args.push_back("-l");

    z = z - Dim[2] - 10;
    x2 = x + Dim[0];
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x2) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "normProjectionsPearsonPvalFDR"); args.push_back("-p"); args.push_back("transnormProjectionsPearsonPvalFDR");
    args.push_back("-as"); args.push_back("normProjectionsPearsonPvalFDR"); args.push_back("-cc"); args.push_back(
      "customLUT_normProjectionsPearsonPvalFDR.txt");

    fiduz = fiduz + 2 * Dim[2] + 10;
    fidux2 = fidux + 2 * Dim[0] + 10;
    fidupos.append(Convert_Double_To_CharArray(fidux2) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);
// -q
    args.push_back("-q"); args.push_back("-id"); args.push_back("normProjectionsPearsonPvalFDR"); args.push_back("-lbl");
    args.push_back("normProjectionsPearsonPvalFDR"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

    fidupos.clear();
    pos.clear();

// ProjectionsSpearmanPvalFDR
    args.push_back("-t"); args.push_back("-f"); args.push_back(
      "./transformFiles/TransnormProjectionsSpearmanPvalFDR.tfm"); args.push_back("-n"); args.push_back(
      "transnormProjectionsSpearmanPvalFDR"); args.push_back("-l");

    //	z=z-Dim[2]-5;
    x = x - Dim[0];
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "normProjectionsSpearmanPvalFDR"); args.push_back("-p"); args.push_back("transnormProjectionsSpearmanPvalFDR");
    args.push_back("-as"); args.push_back("normProjectionsSpearmanPvalFDR"); args.push_back("-cc"); args.push_back(
      "customLUT_normProjectionsSpearmanPvalFDR.txt");

    fidux = fidux - 2 * Dim[0] - 10;
    fidupos.append(Convert_Double_To_CharArray(fidux) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);
// -q
    args.push_back("-q"); args.push_back("-id"); args.push_back("normProjectionsSpearmanPvalFDR"); args.push_back(
      "-lbl"); args.push_back("normProjectionsSpearmanPvalFDR"); args.push_back("-pos"); args.push_back(
      fidupos_all.back().c_str() );

    fidupos.clear();
    pos.clear();

// normProjectionsPearsonPval
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransnormProjectionsPearsonPval.tfm");
    args.push_back("-n"); args.push_back("transnormProjectionsPearsonPval"); args.push_back("-l");

    z = z - Dim[2] - 10;
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x2) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "normProjectionsPearsonPval"); args.push_back("-p"); args.push_back("transnormProjectionsPearsonPval");
    args.push_back("-as"); args.push_back("normProjectionsPearsonPval"); args.push_back("-cc"); args.push_back(
      "customLUT_normProjectionsPearsonPval.txt");

    fiduz = fiduz + Dim[2] + 5;
    fidupos.append(Convert_Double_To_CharArray(fidux2) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("transnormProjectionsPearsonPval"); args.push_back(
      "-lbl"); args.push_back("transnormProjectionsPearsonPval"); args.push_back("-pos"); args.push_back(
      fidupos_all.back().c_str() );

    fidupos.clear();
    pos.clear();

// normProjectionsSpearmanPval
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransnormProjectionsSpearmanPval.tfm");
    args.push_back("-n"); args.push_back("transnormProjectionsSpearmanPval"); args.push_back("-l");

    // z=z-Dim[2]-5;
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "normProjectionsSpearmanPval"); args.push_back("-p"); args.push_back("transnormProjectionsSpearmanPval");
    args.push_back("-as"); args.push_back("normProjectionsSpearmanPval"); args.push_back("-cc"); args.push_back(
      "customLUT_normProjectionsSpearmanPval.txt");

    fidupos.append(Convert_Double_To_CharArray(fidux) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("transnormProjectionsSpearmanPval"); args.push_back(
      "-lbl"); args.push_back("transnormProjectionsSpearmanPval"); args.push_back("-pos"); args.push_back(
      fidupos_all.back(
        ).c_str() );

    fidupos.clear();
    pos.clear();

// normProjectionsPearson
    args.push_back("-t"); args.push_back("-f"); args.push_back("transformFiles/TransnormProjectionsPearson.tfm");
    args.push_back("-n"); args.push_back("transnormProjectionsPearson"); args.push_back("-l");

    z = z - Dim[2] - 10;
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x2) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "normProjectionsPearson"); args.push_back("-p"); args.push_back("transnormProjectionsPearson"); args.push_back(
      "-as");
    args.push_back("normProjectionsPearson"); args.push_back("-cc"); args.push_back(
      "customLUT_normProjectionsPearson.txt");

    fiduz = fiduz + Dim[2] + 5;
    fidupos.append(Convert_Double_To_CharArray(fidux2) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("transnormProjectionsPearson"); args.push_back("-lbl");
    args.push_back("transnormProjectionsPearson"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

    fidupos.clear();
    pos.clear();

// normProjectionsSpearman
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransnormProjectionsSpearman.tfm");
    args.push_back("-n"); args.push_back("transnormProjectionsSpearman"); args.push_back("-l");

    // z=z-Dim[2]-5;
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) ); pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "normProjectionsSpearman"); args.push_back("-p"); args.push_back("transnormProjectionsSpearman"); args.push_back(
      "-as"); args.push_back("normProjectionsSpearman"); args.push_back("-cc"); args.push_back(
      "customLUT_normProjectionsSpearman.txt");

    fidupos.append(Convert_Double_To_CharArray(fidux) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) ); fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("transnormProjectionsSpearman"); args.push_back("-lbl");
    args.push_back("transnormProjectionsSpearman"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

    fidupos.clear();
    pos.clear();

    }
  else
    {

    char mrmlfile[512];
    std::strcpy(mrmlfile, outbase.c_str() );
    std::strcat(mrmlfile, "_GroupTest_MRMLscene.mrml");

    char nameMeanA[512];
    std::strcpy(nameMeanA, outbase.c_str() );
    std::strcat(nameMeanA, "_meanA.meta");
    stringstream ssa;
    string       s_nameMeanA;
    ssa << nameMeanA;
    ssa >> s_nameMeanA;
    size_t found1;
    found1 = s_nameMeanA.find_last_of("/\\");

    char nameMeanB[512];
    std::strcpy(nameMeanB, outbase.c_str() );
    std::strcat(nameMeanB, "_meanB.meta");
    stringstream ssb;
    string       s_nameMeanB;
    ssb << nameMeanB;
    ssb >> s_nameMeanB;
    s_nameMeanA.erase(0, found1 + 1);
    std::strcpy(nameMeanA, s_nameMeanA.c_str() );
    s_nameMeanB.erase(0, found1 + 1);
    std::strcpy(nameMeanB, s_nameMeanB.c_str() );

    double x, y, z, x1, x2, fidux, fiduy, fiduz, fidux2;
    x = 10;
    y = 0;
    z = -60;
    fidux = x + Dim[4] - 10;
    fiduy = y + Dim[5];
    fiduz = (Dim[7] + Dim[8]) / 2 - z;

    args.push_back("CreateMRML");
    args.push_back(mrmlfile);

    // shapes and transforms
// RawP

    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransRawP.tfm"); args.push_back("-n");
    args.push_back("transRawP"); args.push_back("-l");

    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back("RawP");
    args.push_back("-p"); args.push_back("transRawP"); args.push_back("-as"); args.push_back("RawP"); args.push_back(
      "-cc"); args.push_back("customLUT_RawP.txt");

    fidupos.append(Convert_Double_To_CharArray(fidux) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("RawP"); args.push_back("-lbl"); args.push_back("RawP");
    args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

// FDRP
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransFDRP.tfm"); args.push_back("-n");
    args.push_back("transFDRP"); args.push_back("-l");

    z = z + Dim[2] + 5; x1 = x;
    pos.clear();
    fidupos.clear();
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back("FDRP");
    args.push_back("-p"); args.push_back("transFDRP"); args.push_back("-as"); args.push_back("FDRP"); args.push_back(
      "-cc"); args.push_back("customLUT_FDRP.txt");

    fiduz = fiduz - Dim[2] - 5;
    fidupos.append(Convert_Double_To_CharArray(fidux) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("FDRP"); args.push_back("-lbl"); args.push_back("FDRP");
    args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

// overlays Right
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransMeanOverlayRight.tfm");
    args.push_back("-n"); args.push_back("transMeanOverlayRight"); args.push_back("-l");

    z = z + Dim[2] + 5;
    x2 = x + 2 * Dim[0] + 10;

    pos.clear();
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x2) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameMeanA); args.push_back("-n"); args.push_back(
      "MeanAOverlayRight"); args.push_back("-p"); args.push_back("transMeanOverlayRight"); args.push_back("-as");
    args.push_back("MeanAOverlayRight"); args.push_back("-dc"); args.push_back("1,0,0"); args.push_back("-op");
    args.push_back("0.6");
    args.push_back("-m"); args.push_back("-f"); args.push_back(nameMeanB); args.push_back("-n"); args.push_back(
      "MeanBOverlayRight"); args.push_back("-p"); args.push_back("transMeanOverlayRight"); args.push_back("-as");
    args.push_back("MeanBOverlayRight"); args.push_back("-dc"); args.push_back("0,0,1"); args.push_back("-op");
    args.push_back("0.4");

// overlays Left
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransMeanOverlayLeft.tfm");
    args.push_back("-n"); args.push_back("transMeanOverlayLeft"); args.push_back("-l");

    x = x - 2 * Dim[0] - 10;
    pos.clear();
    fidupos.clear();
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameMeanA); args.push_back("-n"); args.push_back(
      "MeanAOverlayLeft"); args.push_back("-p"); args.push_back("transMeanOverlayLeft"); args.push_back("-as");
    args.push_back("MeanAOverlayLeft"); args.push_back("-dc"); args.push_back("1,0,0"); args.push_back("-op");
    args.push_back("0.4");
    args.push_back("-m"); args.push_back("-f"); args.push_back(nameMeanB); args.push_back("-n"); args.push_back(
      "MeanBOverlayLeft"); args.push_back("-p"); args.push_back("transMeanOverlayLeft"); args.push_back("-as");
    args.push_back("MeanBOverlayLeft"); args.push_back("-dc"); args.push_back("0,0,1"); args.push_back("-op");
    args.push_back("0.6");

    fiduz = fiduz - 1.5 * Dim[2] - 5;
    fidux = fidux + Dim[4] / 2;
    fidupos.append(Convert_Double_To_CharArray(fidux) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("Overlays"); args.push_back("-lbl"); args.push_back(
      "Overlays"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

// meanA
    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransMeanRight.tfm"); args.push_back(
      "-n"); args.push_back("transMeanRight"); args.push_back("-l");

    z = z + Dim[2] + 5;

    fidupos.clear();
    pos.clear();
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x2) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameMeanA); args.push_back("-n"); args.push_back(
      "MeanRight"); args.push_back("-p"); args.push_back("transMeanRight"); args.push_back("-as"); args.push_back(
      "MeanRight"); args.push_back("-dc"); args.push_back("1,0,0");

    fiduz = fiduz - Dim[2] - 5;
    fidupos.append(Convert_Double_To_CharArray(fidux) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    args.push_back("-q"); args.push_back("-id"); args.push_back("Means"); args.push_back("-lbl"); args.push_back(
      "Means"); args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );

// meanB

    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransMeanLeft.tfm"); args.push_back(
      "-n"); args.push_back("transMeanLeft"); args.push_back("-l");

    pos.clear();
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameMeanB); args.push_back("-n"); args.push_back(
      "MeanLeft"); args.push_back("-p"); args.push_back("transMeanLeft"); args.push_back("-as"); args.push_back(
      "MeanLeft");
    args.push_back("-dc"); args.push_back("0,0,1");

// Diff

    args.push_back("-t"); args.push_back("-f"); args.push_back("./transformFiles/TransDiffMagnitude.tfm");
    args.push_back("-n"); args.push_back("transDiffMagnitude"); args.push_back("-l");

    z = z + Dim[2] + 5;
    pos.clear();
    fidupos.clear();
    pos.append("1,0,0,0,1,0,0,0,1,");
    pos.append(Convert_Double_To_CharArray(x1) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(y) );
    pos.append(",");
    pos.append(Convert_Double_To_CharArray(z) );
    pos_all.push_back(pos);
    args.push_back( (pos_all.back() ).c_str() );

    args.push_back("-m"); args.push_back("-f"); args.push_back(nameVTK); args.push_back("-n"); args.push_back(
      "DiffMagnitude"); args.push_back("-p"); args.push_back("transDiffMagnitude"); args.push_back("-as");
    args.push_back(
      "DiffMagnitude"); args.push_back("-cc"); args.push_back("customLUT_DiffMagnitude.txt");

    fiduz = fiduz - 1.5 * Dim[2] - 5;
    fidux = fidux + Dim[4] / 2;
    fidupos.append(Convert_Double_To_CharArray(fidux) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduy) );
    fidupos.append(",");
    fidupos.append(Convert_Double_To_CharArray(fiduz) );
    fidupos_all.push_back(fidupos);

    char minmaxf[512];
    std::strcpy(minmaxf, outbase.c_str() );
    std::strcat(minmaxf, "_DiffMagnitude.txt");
    double min, max;
    minmax(minmaxf, &min, &max);

    std::vector<std::string> vec_min_max;
    meandiff.append("Mean_Difference_Min:");
    meandiff.append(Convert_Double_To_CharArray(min) );
    meandiff.append("_Max:");
    meandiff.append(Convert_Double_To_CharArray(max) );
    vec_min_max.push_back(meandiff);
    std::cout << meandiff << std::endl;
    std::cout << vec_min_max.back() << std::endl;

    args.push_back("-q"); args.push_back("-id"); args.push_back(meandiff.c_str() ); args.push_back("-lbl");
    args.push_back(meandiff.c_str() );
    args.push_back("-pos"); args.push_back(fidupos_all.back().c_str() );
    }
  for( int i = 0; i < args.size(); i++ )
    {
    std::cout << args.at(i) << " ";
    }
  // end
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
      for( int i = 0; i < length; i++ )
        {
        st << data[i];
        }
      string dim = st.str();
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

std::vector<double> SetImageDimensions(char *filename)
{

  vector<double> vectDims;
  // read vtk file
  vtkPolyDataReader *meshin = vtkPolyDataReader::New();
  meshin->SetFileName(filename);

  meshin->Update();

  vtkPolyData *poly = vtkPolyData::New();
  poly = meshin->GetOutput();
  vtkIdType idNumPointsInFile = poly->GetNumberOfPoints();

  vtkPoints * pts;
  double      minCoord[3];
  double      maxCoord[3];

  double *firstCoord;

  // find the max and min coordinates
  pts = poly->GetPoints();
  firstCoord = pts->GetPoint(0);

  minCoord[0] = firstCoord[0];
  minCoord[1] = firstCoord[1];
  minCoord[2] = firstCoord[2];

  maxCoord[0] = firstCoord[0];
  maxCoord[1] = firstCoord[1];
  maxCoord[2] = firstCoord[2];
  for( unsigned int i = 1; i < idNumPointsInFile; i++ )
    {
    double *p;
    p = pts->GetPoint(i);

    if( p[0] <= minCoord[0] )
      {
      minCoord[0] = p[0];
      }

    if( p[1] <= minCoord[1] )
      {
      minCoord[1] = p[1];
      }

    if( p[2] <= minCoord[2] )
      {
      minCoord[2] = p[2];
      }

    if( p[0] >= maxCoord[0] )
      {
      maxCoord[0] = p[0];
      }

    if( p[1] >= maxCoord[1] )
      {
      maxCoord[1] = p[1];
      }

    if( p[2] >= maxCoord[2] )
      {
      maxCoord[2] = p[2];
      }
    }

  // vectDims.clear();
  // find the biggest dimension
  vectDims.push_back(maxCoord[0] - minCoord[0]);
  vectDims.push_back(maxCoord[1] - minCoord[1]);
  vectDims.push_back(maxCoord[2] - minCoord[2]);
  vectDims.push_back(minCoord[0]);
  vectDims.push_back(maxCoord[0]); // 4
  vectDims.push_back(minCoord[1]);
  vectDims.push_back(maxCoord[1]);
  vectDims.push_back(minCoord[2]);
  vectDims.push_back(maxCoord[2]);

  return vectDims;

}

/*
vector <double> GetImageDimensions()
{
  return vectDims;
}*/

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

// TODO use this function to reduce the size of the creation MRML function

//
// add_Shape_MRMLScene(args,"RawP",nameVTK,1,"transRawP","./transformFiles/TransRawP.tfm",pos_all.back(),1,"customLUT_RawP.txt",1,fidupos_all.back());
void add_Shape_MRMLScene(std::vector<const char *>& args, std::string name, std::string mesh, bool transform,
                         std::string nametransf, std::string pathtransf, std::string trans_value, bool colorMap,
                         std::string colorMap_value,
                         bool fidu,
                         std::string fidu_value)
// void add_Shape_MRMLScene(std::string name,std::string mesh,bool transform,std::string trans_value,bool
// colorMap,std::string colorMap_value,bool fidu,std::string fidu_value)
{

  // transf
  if( transform )
    {
    /*pathtransf.append("./transformFiles/");
    pathtransf.append(name);
    pathtransf.append(".tfm");
    nametransf.append("trans");
    nametransf.append(name);
std::cout<<pathtransf<<std::endl;*/
    /*std::cout<<pathtransf<<std::endl;
    std::cout<<nametransf.c_str()<<std::endl;
    std::cout<<trans_value<<std::endl;*/

    args.push_back("-t"); args.push_back("-f"); args.push_back(pathtransf.c_str() ); args.push_back("-n");
    args.push_back(nametransf.c_str() ); args.push_back("-l"); args.push_back(trans_value.c_str() );

    // args.push_back("-t"); args.push_back("-f"); args.push_back(pathtransf.c_str()); args.push_back("-n");
    // args.push_back(nametransf.c_str());args.push_back("-l");args.push_back(trans_value.c_str());
    }

  // shape
  args.push_back("-m"); args.push_back("-f"); args.push_back(mesh.c_str() ); args.push_back("-n"); args.push_back(
    name.c_str() );
  if( transform )
    {
    args.push_back("-p"); args.push_back(nametransf.c_str() );
    }

  if( colorMap )
    {
    args.push_back("-as"); args.push_back(name.c_str() ); args.push_back("-cc"); args.push_back(colorMap_value.c_str() );
    }
  else
    {
    args.push_back("-as"); args.push_back(name.c_str() ); args.push_back("-dc"); args.push_back(colorMap_value.c_str() );
    }

  // fidu
  if( fidu )
    {
    args.push_back("-q"); args.push_back("-id"); args.push_back(name.c_str() ); args.push_back("-lbl"); args.push_back(
      name.c_str() ); args.push_back("-pos"); args.push_back(fidu_value.c_str() );
    }

}

bool DirectoryIsEmpty(const char * path)
{
  DIR *          dir = opendir(path);
  struct dirent *mydir;
  int            nbFiles = -2;

  if( dir != NULL )
    {
    while( (mydir = readdir(dir) ) != NULL )
      {
      nbFiles++;

      }
    }
  else
    {
    cerr << " Cannot read the directory " << path << endl;
    }
  if( nbFiles == 0 )
    {
    return true;
    }
  else
    {
    return false;
    }
}
