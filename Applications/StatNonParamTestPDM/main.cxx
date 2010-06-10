#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>

#include <libgen.h>

#include <itkMeshSpatialObject.h>
#include <itkMesh.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTriangleCell.h>
#include <itkDefaultDynamicMeshTraits.h>
#include <itkProcessObject.h>

#include <vnl/algo/vnl_symmetric_eigensystem.h>

#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataToitkMesh.h"

#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>

#include "computeStatPDM.h"
#include "argio.h"

using namespace std;

#define dimension 3

typedef itk::DefaultDynamicMeshTraits< double , dimension, dimension, double, double, double > MeshTrait;
typedef itk::Mesh<double,dimension, MeshTrait> MeshType;

/** Hold on to the type information specified by the template parameters. */
typedef  MeshType::Pointer              MeshPointer;
typedef  MeshTrait::PointType           MeshPointType;
typedef  MeshTrait::PixelType           MeshPixelType; 
typedef  MeshType::Pointer              MeshPointer;
typedef  MeshType::CellTraits           CellTraits;
typedef  MeshType::PointsContainerPointer PointsContainerPointer;
typedef  MeshType::PointsContainer      PointsContainer;
typedef  MeshType::CellsContainerPointer CellsContainerPointer;
typedef  MeshType::CellsContainer       CellsContainer;
typedef  MeshType::PointType            PointType;
typedef  MeshType::CellType             CellType;
typedef  itk::TriangleCell<CellType>   TriangleType;

typedef itk::MeshSpatialObject<MeshType> MeshSpatialObjectType;
typedef itk::SpatialObjectWriter<3,double,MeshTrait> MeshWriterType;
typedef itk::SpatialObjectReader<3,double,MeshTrait> MeshReaderType;

static int debug = 0;

void
load_MeshList_file(char * filename, char * subjInfo, bool scaleOn, bool signDistOn, 
		   int &numSubjects, int &numFeatures,
		   int * &groupLabel, double * &featureValue, char * outbase) 
{

  const int MAXLINE  = 5000; 
  static char line [ MAXLINE ];

  ifstream datafile(filename,ios::in); 
  if (!datafile.is_open()) {
    cerr << "file does not exist" << filename << endl;
    exit(0);
  }

  numSubjects = 0;

  datafile.clear(); 
  datafile.getline(line,MAXLINE);
  while (!datafile.eof()){
    if (line[0] != '#') {
	numSubjects++;
    }
    datafile.getline(line,MAXLINE);
  }
  if (debug) std::cout << "Num Subjects: " << numSubjects << std::endl;
  groupLabel = new int [ numSubjects ];
  double * scaleFactor = new double [numSubjects];
  std::string * meshFileName = new std::string [numSubjects];
    
  // read the list
  int curLine = 0;
  datafile.clear();
  datafile.seekg(0, std::ios::beg);
  datafile.getline(line,MAXLINE);
  while (!datafile.eof()){
    if (line[0] != '#') {
      char filename[MAXLINE];
      int retval;
      if (subjInfo) {
	retval = sscanf(line," %s ",  filename);
      } else {
	retval = sscanf(line," %d %lf %s ", &(groupLabel[curLine]) , &(scaleFactor[curLine]), filename);
      }
      if (retval > 0 ) {
	meshFileName[curLine] = std::string(filename);
	curLine++;
      }
    }
    datafile.getline(line,MAXLINE);
  }
  datafile.close();

  if (subjInfo){
    ifstream infofile(subjInfo,ios::in); 
    if (!infofile.is_open()) {
      std::cerr << "file does not exist" << subjInfo << std::endl;
      exit(0);
    }
    infofile.clear(); 
    infofile.getline(line,MAXLINE);
    int numSubjectsInfo = 0;
    while (!infofile.eof()){
      if (line[0] != '#') {
	numSubjectsInfo++;
      }
      infofile.getline(line,MAXLINE);
    }    
    int * groupLabelInfo = new int [ numSubjectsInfo ];
    double * scaleFactorInfo = new double [numSubjectsInfo];
    std::string * subjectIDInfo = new std::string [numSubjectsInfo];

    // Read all the subject info
    int curLine = 0;
    infofile.clear();
    infofile.seekg(0, std::ios::beg);
    infofile.getline(line,MAXLINE);
    while (!infofile.eof()){
      if (line[0] != '#') {
	char subjectID[MAXLINE];
        int retval = sscanf(line," %s %d %lf ", subjectID, &(groupLabelInfo[curLine]) , &(scaleFactorInfo[curLine]));
	if (retval > 0) {
	  subjectIDInfo[curLine] = std::string(subjectID);
	  curLine++;
	}
      }
      infofile.getline(line,MAXLINE);

    }
    infofile.close();

    // match the subject info to the mesh files
    for (int subjIndex = 0; subjIndex < numSubjects; subjIndex++) {
      bool foundMatch = false;
      for (int infoIndex = 0; infoIndex < numSubjectsInfo; infoIndex++) {
	if (strstr(meshFileName[subjIndex].c_str(), subjectIDInfo[infoIndex].c_str())) {
	  if (foundMatch) {
	    std::cout << "Double match for mesh file name with multiple subject id's : " << meshFileName[subjIndex] 
		      << "-" << subjectIDInfo[infoIndex] << std::endl;
	    exit(-1);
	  } else {
	    //std::cout << " match  " << meshFileName[subjIndex] << "-" << subjectIDInfo[infoIndex] << std::endl;
	    foundMatch = true;
	    groupLabel[subjIndex] = groupLabelInfo[infoIndex];
	    scaleFactor[subjIndex] = scaleFactorInfo[infoIndex];
	  }
	}
      }
      if (!foundMatch) {
	std::cout << "Found no match for mesh file" << meshFileName[subjIndex] << std::endl;
	exit(-2);
      }
    }
  }
  // debug info
  if (debug) {
    for (int i = 0; i < numSubjects; i++) {
      std::cout << meshFileName[i] << " " << groupLabel[i] << " " << scaleFactor[i] << std::endl;
    }
  }


  // Read the meshes
  MeshType::Pointer surfaceMesh = MeshType::New();
  MeshSpatialObjectType::Pointer  SOMesh;
  MeshReaderType::Pointer reader = MeshReaderType::New();
  PointsContainerPointer points;
  for (int index = 0; index < numSubjects; index++) {
    try
      {
       reader->SetFileName(meshFileName[index].c_str());
       reader->Update();
       MeshReaderType::SceneType::Pointer scene = reader->GetScene();  
       MeshReaderType::SceneType::ObjectListType * objList =  scene->GetObjects(1,NULL);
    
       MeshReaderType::SceneType::ObjectListType::iterator it = objList->begin();
       itk::SpatialObject<3> * curObj = *it;
       SOMesh = dynamic_cast<MeshSpatialObjectType*> (curObj);
       surfaceMesh = SOMesh->GetMesh();
       points = surfaceMesh->GetPoints();

       if (index == 0) {
         numFeatures = points->Size() * 3;
         featureValue = new double [ numSubjects * numFeatures ];
       }
    
       for (unsigned int pointID = 0; pointID < points->Size();pointID++) {
         PointType curPoint =  points->GetElement(pointID);
         for (unsigned int dim = 0; dim < 3; dim++) {
	   if (scaleOn) {
	     featureValue[index*numFeatures + (pointID * 3) + dim] = curPoint[dim] / scaleFactor[index];
	   } else {
	     featureValue[index*numFeatures + (pointID * 3) + dim] = curPoint[dim];
	   }
         }
       }
      }
    catch(itk::ExceptionObject ex)
      {
	std::cout<< "Error reading meshfile:  "<< meshFileName[index] << std::endl << "ITK error: " << ex.GetDescription()<< std::endl;
           exit(-3);
    }
  }

  // change labels to the predefined ones, this will fail if there are more than 2 labels in the file
  int preLabelA = groupLabel[0];
  int preLabelB = groupLabel[0];
  int i;
  for (i = 0; i < numSubjects; i++) {
    if (preLabelA != groupLabel[i]) {
      if (preLabelB != preLabelA && preLabelB  != groupLabel[i]) {
	std::cout << "Error: more than 2 labels in file" << std::endl;
      } else {
          preLabelB = groupLabel[i];
      }
    }
  }
  for (i = 0; i < numSubjects; i++) {
    if (preLabelA == groupLabel[i]) {
      groupLabel[i] = GROUP_A_LABEL ;
    } else if (preLabelB == groupLabel[i]) {
      groupLabel[i] = GROUP_B_LABEL ;
    }
  }
  if (debug) {
    std::cout << "data has been relabeled: " <<  preLabelA << " --> group A = " << GROUP_A_LABEL
	      << " ; " << preLabelB << " --> group B = " << GROUP_B_LABEL << std::endl;
  }

  // compute and save averages
  static double * meanA = new double [numFeatures];
  static double * meanB = new double [numFeatures];
  int numSubjA = 0;
  int numSubjB = 0;
  for (int feat = 0; feat < numFeatures; feat++) {
    meanA[feat] = 0;
    meanB[feat] = 0;
  }
  for (int subj = 0; subj < numSubjects; subj++) {
    int subjIndex = subj * numFeatures;
    if (groupLabel[subj] == GROUP_A_LABEL) {
      numSubjA++;
      for (int feat = 0; feat < numFeatures; feat++) {
        meanA[feat] = meanA[feat] + featureValue[subjIndex + feat];
      }
    } else if (groupLabel[subj] == GROUP_B_LABEL) {
      numSubjB++;
      for (int feat = 0; feat < numFeatures; feat++) {
        meanB[feat] = meanB[feat] + featureValue[subjIndex + feat];
      }
    } else {
      std::cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
    }
  }
  for (int feat = 0; feat < numFeatures; feat++) {
    meanA[feat] = meanA[feat] / numSubjA;
    meanB[feat] = meanB[feat] / numSubjB;
  }
  if (debug) std::cout << "# A: " << numSubjA << ", #B: " << numSubjB << std::endl;
  // adapt the last loaded mesh
  surfaceMesh = SOMesh->GetMesh();
  PointsContainerPointer pointsA = PointsContainer::New();
  PointsContainerPointer pointsB = PointsContainer::New();
    
  for (unsigned int i = 0; i < points->Size(); i++) {
    double vert[3];
    vert[0] = meanA[i*3 + 0]; vert[1] = meanA[i*3 + 1]; vert[2] = meanA[i*3 + 2];
    pointsA->InsertElement(i, PointType(vert));
    vert[0] = meanB[i*3 + 0]; vert[1] = meanB[i*3 + 1]; vert[2] = meanB[i*3 + 2];
    pointsB->InsertElement(i, PointType(vert));
  }
  surfaceMesh->SetPoints(pointsA); 
  SOMesh->SetMesh(surfaceMesh);
  MeshWriterType::Pointer writer = MeshWriterType::New();
  writer->SetInput(SOMesh);
  string FilenameA(outbase);
  FilenameA = FilenameA + string("_meanA.meta");
  writer->SetFileName(FilenameA.c_str());
  writer->Update();
  surfaceMesh->SetPoints(pointsB); 
  SOMesh->SetMesh(surfaceMesh);
  writer->SetInput(SOMesh);
  string FilenameB(outbase);
  FilenameB = FilenameB + string("_meanB.meta");
  writer->SetFileName(FilenameB.c_str());
  writer->Update();
  
  ofstream diffFile ;

  // mesh Difference as Vector map
  string FilenameDiff(outbase);
  FilenameDiff = FilenameDiff + string("_meanDiffVec.txt");
  diffFile.open ( FilenameDiff.c_str()) ;
  diffFile << "NUMBER_OF_POINTS = " << points->Size() << endl ;
  diffFile << "DIMENSION = 3" << endl ;
  diffFile << "TYPE = Vector" << endl ;
  for (unsigned int i = 0 ; i < points->Size(); i++ )
  {
    diffFile << meanB[i*3 + 0] - meanA[i*3 + 0] << " " << meanB[i*3 + 1] - meanA[i*3 + 1] << " " 
	     << meanB[i*3 + 2] - meanA[i*3 + 2] << endl ;
  }
  diffFile.close () ;

  // compute Covariance ellipsoids of groupA
  string FilenameCovarA(outbase);
  FilenameCovarA = FilenameCovarA + string("_CovarA.txt");
  diffFile.open ( FilenameCovarA.c_str()) ;
  diffFile << "NUMBER_OF_POINTS = " << points->Size() << endl ;
  diffFile << "DIMENSION = 9" << endl ;
  diffFile << "TYPE = Ellipsoid" << endl ;
  
  for (int feat = 0 ; feat < numFeatures / MeasurementVectorSize; feat++ )
  {
    SampleType::Pointer sample = SampleType::New();  
    sample->SetMeasurementVectorSize( MeasurementVectorSize );

    MeasurementVectorType mv;
    for (int subj = 0 ; subj < numSubjects; subj++) {
      if (groupLabel[subj] == GROUP_A_LABEL) {
	for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
	  mv[dim] = featureValue[subj*numFeatures + (feat * MeasurementVectorSize) + dim];
	}      
	sample->PushBack(mv) ;
      }
    }

    CovarianceCalculatorType::Pointer covarianceCalculator = CovarianceCalculatorType::New() ;
    covarianceCalculator->SetInputSample(sample.GetPointer()) ;
    covarianceCalculator->Update();
    
    const CovarianceCalculatorType::OutputType*  outMatrix = covarianceCalculator->GetOutput();

    // compute ellipsoid by solving the eigensystem of the covariance matrix
    EigenSystemType EigenSystem;
    MeasurementMatrixType EigenVectorMatrix, resMatrix;
    MeasurementVectorType EigenValues;
    for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
      for (int dim2 = 0 ; dim2 < MeasurementVectorSize; dim2++) {
	resMatrix[dim][dim2] = (float) (*outMatrix)[dim][dim2];
      }
    }
    EigenSystem.SetDimension(MeasurementVectorSize);
    EigenSystem.ComputeEigenValuesAndVectors(resMatrix, EigenValues, EigenVectorMatrix);
    
    // the Eigenvalues are sorted smallest to largest! Who came up with this sorting?
    // write eigenvectors * eigenvalue
    
    for (int dim = MeasurementVectorSize - 1  ; dim >= 0; dim--) {
      for (int dim2 = 0 ; dim2 < MeasurementVectorSize; dim2++) {
	diffFile << EigenVectorMatrix[dim][dim2] * sqrt(EigenValues[dim]) << " ";
      }
    }
    diffFile << std::endl;
    
  }
  diffFile.close () ;
  
  // compute Covariance ellipsoids of groupB
  string FilenameCovarB(outbase);
  FilenameCovarB = FilenameCovarB + string("_CovarB.txt");
  diffFile.open ( FilenameCovarB.c_str()) ;
  diffFile << "NUMBER_OF_POINTS = " << points->Size() << endl ;
  diffFile << "DIMENSION = 9" << endl ;
  diffFile << "TYPE = Ellipsoid" << endl ;
  
  for (int feat = 0 ; feat < numFeatures / MeasurementVectorSize; feat++ )
  {
    SampleType::Pointer sample = SampleType::New();  
    sample->SetMeasurementVectorSize( MeasurementVectorSize );

    MeasurementVectorType mv;
    for (int subj = 0 ; subj < numSubjects; subj++) {
      if (groupLabel[subj] == GROUP_B_LABEL) {
	for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
	  mv[dim] = featureValue[subj*numFeatures + (feat * MeasurementVectorSize) + dim];
	}      
	sample->PushBack(mv) ;
      }
    }
    
    CovarianceCalculatorType::Pointer covarianceCalculator = CovarianceCalculatorType::New() ;
    covarianceCalculator->SetInputSample(sample.GetPointer()) ;
    covarianceCalculator->Update();
    
    const CovarianceCalculatorType::OutputType*  outMatrix = covarianceCalculator->GetOutput();
    
    // compute ellipsoid by solving the eigensystem of the covariance matrix
    EigenSystemType EigenSystem;
    MeasurementMatrixType EigenVectorMatrix, resMatrix;
    MeasurementVectorType EigenValues;
    for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
      for (int dim2 = 0 ; dim2 < MeasurementVectorSize; dim2++) {
	resMatrix[dim][dim2] = (float) (*outMatrix)[dim][dim2];
      }
    }
    
    EigenSystem.SetDimension(MeasurementVectorSize);
    EigenSystem.ComputeEigenValuesAndVectors(resMatrix, EigenValues, EigenVectorMatrix);
    
    // the Eigenvalues are sorted smallest to largest! Who came up with this sorting? 
    // write eigenvectors * eigenvalue
    
    for (int dim = MeasurementVectorSize - 1  ; dim >= 0; dim--) {
      for (int dim2 = 0 ; dim2 < MeasurementVectorSize; dim2++) {
	diffFile << EigenVectorMatrix[dim][dim2] * sqrt(EigenValues[dim]) << " ";
      }
    }
    diffFile << std::endl;
    
  }
  diffFile.close () ;


  // compute Covariance ellipsoids and Mean of joined population
  string FilenameCovar(outbase);
  FilenameCovar = FilenameCovar + string("_CovarAll.txt");
  diffFile.open ( FilenameCovar.c_str()) ;
  diffFile << "NUMBER_OF_POINTS = " << points->Size() << endl ;
  diffFile << "DIMENSION = 9" << endl ;
  diffFile << "TYPE = Ellipsoid" << endl ;

  PointsContainerPointer pointsAB = PointsContainer::New();
  for (int feat = 0 ; feat < numFeatures / MeasurementVectorSize; feat++ )
  {
    SampleType::Pointer sample = SampleType::New();  
    sample->SetMeasurementVectorSize( MeasurementVectorSize );

    MeasurementVectorType mv;
    for (int subj = 0 ; subj < numSubjects; subj++) {
      for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
	mv[dim] = featureValue[subj*numFeatures + (feat * MeasurementVectorSize) + dim];
      }      
      sample->PushBack(mv) ;
    }
    MeanCalculatorType::Pointer meanCalculator = MeanCalculatorType::New() ;
    meanCalculator->SetInputSample(sample.GetPointer());
    meanCalculator->Update() ;

    double vert[MeasurementVectorSize];
    MeanCalculatorType::OutputType* meanOutput = meanCalculator->GetOutput();
    for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
      vert[dim] = (*meanOutput)[dim];
    }
    pointsAB->InsertElement(feat, PointType(vert));
    
    CovarianceCalculatorType::Pointer covarianceCalculator = CovarianceCalculatorType::New() ;
    covarianceCalculator->SetInputSample(sample.GetPointer()) ;
    covarianceCalculator->SetMean(meanCalculator->GetOutput()) ;
    covarianceCalculator->Update();
    
    const CovarianceCalculatorType::OutputType*  outMatrix = covarianceCalculator->GetOutput();

    // compute ellipsoid by solving the eigensystem of the covariance matrix
    EigenSystemType EigenSystem;
    MeasurementMatrixType EigenVectorMatrix, resMatrix;
    MeasurementVectorType EigenValues;
    for (int dim = 0 ; dim < MeasurementVectorSize; dim++) {
      for (int dim2 = 0 ; dim2 < MeasurementVectorSize; dim2++) {
	resMatrix[dim][dim2] = (float) (*outMatrix)[dim][dim2];
      }
    }
    EigenSystem.SetDimension(MeasurementVectorSize);
    EigenSystem.ComputeEigenValuesAndVectors(resMatrix, EigenValues, EigenVectorMatrix);
    
    // the Eigenvalues are sorted smallest to largest! Who came up with this sorting? It should be reverse...
    // write eigenvectors * eigenvalue
    
    for (int dim = MeasurementVectorSize - 1  ; dim >= 0; dim--) {
      for (int dim2 = 0 ; dim2 < MeasurementVectorSize; dim2++) {
	diffFile << EigenVectorMatrix[dim][dim2] * sqrt(EigenValues[dim]) << " ";
      }
    }
    diffFile << std::endl;
    
  }
  diffFile.close () ;
  
  // write mesh
  surfaceMesh->SetPoints(pointsAB); 
  SOMesh->SetMesh(surfaceMesh);
  writer->SetInput(SOMesh);
  string FilenameAB(outbase);
  FilenameAB = FilenameAB + string("_meanAll.meta");
  writer->SetFileName(FilenameAB.c_str());
  writer->Update();

  if (signDistOn) {
    // compute features to be the signed distance to the overall mean 
    // the sign is determined by the surface normal of the mean surface
    int numFeaturesDist = pointsAB->Size();
    double * featureValueDist = new double [ numSubjects * numFeaturesDist ];
    
    // compute mean surface normals
    if (debug)   std::cout << "Computing normals" << std::endl;

    itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
    convertMeshToVTK->SetInput(surfaceMesh);
    vtkPolyData * vtkMesh = convertMeshToVTK->GetOutput();
    vtkPolyDataNormals *MeshNormals = vtkPolyDataNormals::New();

    MeshNormals->SetComputePointNormals(1);
    MeshNormals->SetComputeCellNormals(0);
    MeshNormals->SetSplitting(0);
    MeshNormals->SetInput(vtkMesh);
    MeshNormals->Update();
    vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
    vtkMeshNormals->Update();

    vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
    vtkDataArray * ArrayNormal = NormalPoints->GetNormals();

    if (debug)   std::cout << "Computing normal projections for all meshes" << std::endl;
    for (int subj = 0; subj < numSubjects; subj++) {
      int subjIndex = subj * numFeatures;
      int subjIndexDist = subj * numFeaturesDist;
      for (int featDist = 0; featDist < numFeaturesDist; featDist++) {
	//val = (point_i - meanAB_i) . normalMeanAB_i
	double diffVec[3], * Normal;
	PointType curPoint =  pointsAB->GetElement(featDist);
	Normal = ArrayNormal->GetTuple3(featDist);
	for (int dim = 0; dim < 3; dim++) {
	  diffVec[dim] = featureValue[subjIndex + featDist * 3 + dim] - curPoint[dim];
	}
	featureValueDist[subjIndexDist + featDist] = 
	  diffVec[0] * Normal[0] + diffVec[1]  * Normal[1] + diffVec[2] * Normal[2];
      }
      // write feature maps as KWMeshVisu files
      string FilenameMag(outbase);
      char * basestring = strdup(basename((char *) meshFileName[subj].c_str()));
      char * ext = strrchr(basestring , '.');
      if (ext) ext[0] = '\0';
      FilenameMag = FilenameMag + string(basestring) + string("_meanSignDiff.txt");
      diffFile.open ( FilenameMag.c_str()) ;
      diffFile << "NUMBER_OF_POINTS = " << pointsAB->Size() << endl ;
      diffFile << "DIMENSION = 1" << endl ;
      diffFile << "TYPE = Scalar" << endl ;
      for (int featDist = 0; featDist < numFeaturesDist; featDist++) 
	{
	  diffFile << featureValueDist[subjIndexDist + featDist] << endl ;
	}
      diffFile.close () ;
    }
    
    numFeatures = numFeaturesDist;
    featureValue = featureValueDist;

  }

  return;
}

void
load_simpleStat_file(char * filename, int featSelStart, 
           int featSelLen, int featgroupID,
           int &numSubjects, int &numFeatures,
		     int * &groupLabel, double * &featureValue, char * outbase) 
{
  const int MAXLINE  = 10 * 50000;  // maximal 50'000 entries each 10 digits long
  static char line [ MAXLINE ];
  //int subjIDindex = 0;
  //int groupNameindex = 1;
  int groupIDindex = 2;
  int firstDataindex = 3;
  int endDataindex = -1;
  int featureLength = -1;

  if (featgroupID >= 0 ) groupIDindex = featgroupID;
  if (featSelStart >= 0 ) firstDataindex = featSelStart;
  if (featSelLen >= 0 ) featureLength = featSelLen;

  ifstream datafile(filename,ios::in); 
  if (!datafile.is_open()) {
    cerr << "file does not exist" << filename << endl;
    exit(0);
  }

  datafile.getline(line,MAXLINE);
  while (line[0] == '#' && ! datafile.eof()) {
    datafile.getline(line,MAXLINE);
  }
  if (datafile.eof()){
    cerr << "not a single line in file " << filename << endl;
    exit(0);
  }
  
  int cnt = 0;
  numFeatures = 0;
  while (line[cnt] != '\0') {
    // skip delimiters
    while  (line[cnt] != '\0' && 
         ( line[cnt] == ' ' || line[cnt] == ',' || line[cnt] == ';' || 
           line[cnt] == '\t')) {
      cnt++;
    }
    if (line[cnt] != '\0') numFeatures++;
    // skip non-delimiters
    while  (line[cnt] != '\0' && 
         line[cnt] != ' ' && line[cnt] != ',' && line[cnt] != ';' && 
         line[cnt] != '\t') {
      cnt++;
    }
  }
  
  if (featureLength > 0 ) numFeatures = featureLength;
  else numFeatures = numFeatures - firstDataindex;

  endDataindex = firstDataindex + numFeatures - 1;
  if (debug) std::cout << "Num Features: " << numFeatures << " from  " << firstDataindex << " to " << endDataindex << std::endl;

  numSubjects = 0;

  datafile.clear(); 
  datafile.seekg(0, std::ios::beg);
  datafile.getline(line,MAXLINE);
  while (!datafile.eof()){
    if (line[0] != '#') {
      numSubjects++;
    }
    datafile.getline(line,MAXLINE);
  }

  if (debug) std::cout << "Num Subjects: " << numSubjects << std::endl;

  groupLabel = new int [ numSubjects ];
  featureValue = new double [ numSubjects * numFeatures ];

  int curLine = 0;
  datafile.clear();
  datafile.seekg(0, std::ios::beg);
  datafile.getline(line,MAXLINE);
  while (!datafile.eof()){
    // scan features per line
    if (line[0] != '#') {
      int cnt = 0;
      int curFeature = 0;
      while (line[cnt] != '\0' && curFeature <= endDataindex) {
	// skip delimiters
	while  (line[cnt] != '\0' && 
		( line[cnt] == ' ' || line[cnt] == ',' || line[cnt] == ';' || 
		  line[cnt] == '\t')) {
	  cnt++;
	}
	// read in feature
      bool numberRead = false;
      int numChar = 0;
      if (line[cnt] != '\0') {
        if (curFeature == groupIDindex) {
          //scan in groupID
          int label = -1;
          int numRead = sscanf(&(line[cnt]),"%d%n",&label, &numChar);
          if (!numRead || !numChar) { std::cout << "Error while reading " << std::endl;}
          groupLabel[curLine] = label;
	  numberRead = true;
        } else if (curFeature >= firstDataindex && curFeature <= endDataindex) {
          //scan in feature
          double feature = -1;
          int numRead = sscanf(&(line[cnt]),"%lf%n",&feature, &numChar);
          if (!numRead || !numChar) { std::cout << "Error while reading " << std::endl;}
          int featIndex = curFeature - firstDataindex;
          featureValue[curLine*numFeatures + featIndex] = feature;
	  numberRead = true;
        }
	cnt = cnt + numChar;
	curFeature++;
      } 
      // skip non-delimiters
      while  (line[cnt] != '\0' && 
	      line[cnt] != ' ' && line[cnt] != ',' && line[cnt] != ';' && 
	      line[cnt] != '\t') {
	if (numberRead) { 
	  // we were reading a number and now there are just adjacent to this number other characters -> error
	  // when reading number, then the next character needs to be a delimiter
	  std::cout << "Error while reading number, encountering non-number on line " << curLine 
	       << " feature: " << curFeature << " numChar : " << numChar << std::endl;
	}
        cnt++;
      }
      }
      if (curFeature <= endDataindex) { std::cout << "Not enough features" << std::endl;}
      // next line
      curLine++;
    }
    datafile.getline(line,MAXLINE);  
  }
  datafile.close();

  // change labels to the predefined ones, this will fail if there are more than 2 labels in the file
  int preLabelA = groupLabel[0];
  int preLabelB = groupLabel[0];
  int i;
  for (i = 0; i < numSubjects; i++) {
    if (preLabelA != groupLabel[i]) {
      if (preLabelB != preLabelA && preLabelB  != groupLabel[i]) {
          std::cout << "Error more than 2 labels in file" << std::endl;
      } else {
                      preLabelB = groupLabel[i];
      }
    }
  }
  for (i = 0; i < numSubjects; i++) {
    if (preLabelA == groupLabel[i]) {
      groupLabel[i] = GROUP_A_LABEL ;
    } else if (preLabelB == groupLabel[i]) {
      groupLabel[i] = GROUP_B_LABEL ;
    }
  }
  if (debug){
    std::cout << "data has been relabeled: " <<  preLabelA << " --> " << GROUP_A_LABEL
	      << " ; " << preLabelB << " --> " << GROUP_B_LABEL << std::endl;
  }

  // compute and save averages
  static double * meanA = new double [numFeatures];
  static double * meanB = new double [numFeatures];
  static double * meanAB = new double [numFeatures];
  int numSubjA = 0;
  int numSubjB = 0;
  int numSubjAB = 0;
  for (int feat = 0; feat < numFeatures; feat++) {
    meanA[feat] = 0;
    meanB[feat] = 0;
    meanAB[feat] = 0;
  }
  for (int subj = 0; subj < numSubjects; subj++) {
    int subjIndex = subj * numFeatures;
    numSubjAB++;
    for (int feat = 0; feat < numFeatures; feat++) {
      meanAB[feat] = meanAB[feat] + featureValue[subjIndex + feat];
    }
    if (groupLabel[subj] == GROUP_A_LABEL) {
      numSubjA++;
      for (int feat = 0; feat < numFeatures; feat++) {
        meanA[feat] = meanA[feat] + featureValue[subjIndex + feat];
      }
    } else if (groupLabel[subj] == GROUP_B_LABEL) {
      numSubjB++;
      for (int feat = 0; feat < numFeatures; feat++) {
        meanB[feat] = meanB[feat] + featureValue[subjIndex + feat];
      }
    } else {
      std::cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
    }
  }
  for (int feat = 0; feat < numFeatures; feat++) {
    meanA[feat] = meanA[feat] / numSubjA;
    meanB[feat] = meanB[feat] / numSubjB;
    meanAB[feat] = meanAB[feat] / numSubjAB;
  }
  if (debug) std::cout << "# A: " << numSubjA << ", #B: " << numSubjB << std::endl;

  ofstream meanFile ;
  string FilenameMeanA(outbase);
  FilenameMeanA = FilenameMeanA + string("_meanA.txt");
  meanFile.open ( FilenameMeanA.c_str()) ;
  meanFile << "NUMBER_OF_POINTS = " << numFeatures << endl ;
  meanFile << "DIMENSION = 1" << endl ;
  meanFile << "TYPE = Scalar" << endl ;
  for (unsigned int i = 0 ; i < numFeatures; i++ )
  {
    meanFile <<  meanA[i]  << endl ;
  }
  meanFile.close () ;
  string FilenameMeanB(outbase);
  FilenameMeanB = FilenameMeanB + string("_meanB.txt");
  meanFile.open ( FilenameMeanB.c_str()) ;
  meanFile << "NUMBER_OF_POINTS = " << numFeatures << endl ;
  meanFile << "DIMENSION = 1" << endl ;
  meanFile << "TYPE = Scalar" << endl ;
  for (unsigned int i = 0 ; i < numFeatures; i++ )
  {
    meanFile <<  meanB[i] << endl ;
  }
  meanFile.close () ;
  string FilenameMeanAB(outbase);
  FilenameMeanAB = FilenameMeanAB + string("_meanAB.txt");
  meanFile.open ( FilenameMeanAB.c_str()) ;
  meanFile << "NUMBER_OF_POINTS = " << numFeatures << endl ;
  meanFile << "DIMENSION = 1" << endl ;
  meanFile << "TYPE = Scalar" << endl ;
  for (unsigned int i = 0 ; i < numFeatures; i++ )
  {
    meanFile <<  meanAB[i] << endl ;
  }
  meanFile.close () ;
}


int main (int argc, const char ** argv) 
{
  if (argc <=1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
    std::cout << argv[0] << std::endl;
    std::cout << " computing  local statistics corrected for multiple comparisons via non-parametric permutations tests" << endl;
    std::cout << "  the metric for univariate data is mean difference and for multivariate data Hotelling T^2" << std::endl;
    std::cout << "usage:  infile [options] "<< std::endl;
    std::cout << std::endl;
    std::cout << "infile        input file " << std::endl;
    std::cout << "-surfList     if present, then the input file does not contain the data, but rather is a list of" << std::endl;
    std::cout << "              lines each with an itk readable mesh file" << std::endl;
    std::cout << "-scale        enable scaling of surfaces, this option is only used if -surfList is specified, scaling is applied by division with the individual scaling information per mesh" << std::endl;
    std::cout << "-signDist     compute statistics using univariate signed distance to the overall mean" << std::endl;
    std::cout << "-subjInfoFile <file> if present then the information about group ID and scaling factor is stored in this file" << std::endl;
    std::cout << "              each subject have a single line with information separated by white spaces." << std::endl;
    std::cout << "              the first column has to be the subjectID which is also encoded in the meshfilename" << std::endl;
    std::cout << "              the second column is the group ID and the third the scale information/divisor, any other column is disregarded" << std::endl;
    std::cout << "              if this option is not present and -surfList is used then each meshfile has to be preceded by the" << std::endl;
    std::cout << "              group id (number) and the scaling factor" << std::endl; 
    std::cout << "-out <base>   output base name, if omitted then infile is used " << std::endl << std::endl;
    std::cout << "-featgroupID  column number of group id" << std::endl;
    std::cout << "-featSelStart column number for first feature" << std::endl;
    std::cout << "-featSelLen   number of consecutive columns with features" << std::endl;
    std::cout << "              indexes for featSelStart and -featgroupID begn at 0" << std::endl << std::endl;
    std::cout << "-log          perform log(1+p) transform necessary for magnitude data (all features are transformed)" << std::endl;
    std::cout << "-numPerms <i> number of permutations" << std::endl;
    std::cout << "-signLevel <f> maximum level of significance,choose 1 if all levels should be computed" << std::endl;
    std::cout << "-signSteps <i> number of steps from 0..signLevel" << std::endl;
    std::cout << "-3DvectorData  3D statistics with Hotelling T^2 "<< std::endl;
    std::cout << "-2DvectorData  2D statistics with Hotelling T^2 "<< std::endl;
    std::cout << "-v             verbose mode " << std::endl;
    std::cout << std::endl << std::endl;
    std::cout << " the uncorrrected p-values will be written out as well as the corrected p-vals" << std::endl;
    std::cout << std::endl << std::endl;
    exit(0);
  }

  char *infile     = strdup(argv[1]);
  char *outbase    = ipGetStringArgument(argv, "-out", NULL);
  if (!infile) {
    std::cout << " no input file specified " << std::endl;exit(1);
    }
  if (!outbase) {
    outbase = strdup(infile);
  } 

  const int numPerms = ipGetIntArgument(argv, "-numPerms", 10000);
  const double significanceLevel = ipGetDoubleArgument(argv, "-signLevel", 0.05);
  const int significanceSteps = ipGetIntArgument(argv, "-signSteps", 100);

  const int featSelStart  = ipGetIntArgument(argv, "-featSelStart", -1);
  const int featSelLen  = ipGetIntArgument(argv, "-featSelLen", -1);
  const int featgroupID  = ipGetIntArgument(argv, "-featgroupID", -1);

  int logOn = ipExistsArgument(argv, "-log");

  int numSubjects, numFeatures;
  int * groupLabel = NULL;
  double * featureValue = NULL; 

  bool surfListInput = ipExistsArgument(argv,"-surfList");
  bool surfListScale = ipExistsArgument(argv,"-scale");
  bool surfSignDist = ipExistsArgument(argv,"-signDist");
  char * surfListSubjInfo = ipGetStringArgument(argv, "-subjInfoFile", NULL); 

  bool vector2DDataFlag = ipExistsArgument(argv, "-2DvectorData");
  bool vector3DDataFlag = ipExistsArgument(argv, "-3DvectorData");

  debug = ipExistsArgument(argv, "-v");

  if (!surfListInput) {
    load_simpleStat_file(infile,featSelStart,featSelLen,featgroupID,numSubjects, 
			 numFeatures,groupLabel,featureValue, outbase);
  } else {
    // load the surfListInput
    if (! surfSignDist) {
      vector3DDataFlag = true;
    } 
    logOn = false;
    load_MeshList_file(infile, surfListSubjInfo, surfListScale, surfSignDist, numSubjects, numFeatures,groupLabel, featureValue, outbase);
  }

  if (logOn) {
    // log transform the Feature data
    for (int subj = 0; subj < numSubjects; subj++) {
      for (int feat = 0; feat < numFeatures; feat++) {
        int subjIndex = subj * numFeatures;
        if (featureValue[subjIndex + feat] < 0) {
          std::cerr << "Features are negative, and thus do not represent magnitude data" << endl;
          std::cerr << "Cowardly refusing to log transform data... " << std::endl;
          exit(-4);
        }
        featureValue[subjIndex + feat] = log (1 + featureValue[subjIndex + feat]);
      }
    }
  }

  int tupel_size = 1;
  if (vector3DDataFlag) {
    tupel_size = 3;
  } else  if (vector2DDataFlag) {
    tupel_size = 2;
  } else {
    tupel_size = 1;
  }

  if (debug) std::cout << "tupel: " << tupel_size << " , numFeatures: "  << numFeatures << std::endl;

  double * pValue = new double [numFeatures];
  
  doTesting(numSubjects, numFeatures, numPerms, tupel_size, groupLabel, featureValue,
        significanceLevel, significanceSteps,pValue, outbase);

  return 0; 
}
