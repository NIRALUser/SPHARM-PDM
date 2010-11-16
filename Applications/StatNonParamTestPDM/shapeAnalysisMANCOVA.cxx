/***********************************
*** Written by Marc Macenko      ***
*** 09/25/2008                   ***
*** Modified by Marc Niethammer  ***
*** 05/2009                      ***
*** Modified by Beatriz Paniagua ***
*** 08/2009 and 09/2009          ***
************************************/

#include "MANCOVA.h"
#include "miscMANCOVA.h"
#include "nonMANCOVA.h"

#include <iostream>
#include <fstream>

#include <string>
#include <time.h>
#include "shapeAnalysisMANCOVACLP.h"

#include "constants.h"
#include "newTypes.h"

// TODO: It may be useful to introduce some structures
// to pass around all the important information
// this should make the function interfaces much less
// cumbersome (MN)

void doMANCOVATesting( unsigned int numSubjects, unsigned int numFeatures, unsigned int numGroupTypes, unsigned int numIndependent, unsigned int testColumn, unsigned int numPerms, vnl_matrix<int> * groupLabel, vnl_matrix<double>  * &featureValue, bool interactionTest, unsigned int numA, unsigned int numB, double significanceLevel, int computeStatisticsType, vnl_vector<double>& rawP, PointsContainerPointer &meanPoints  );

int main(int argc, char *argv[])
{std::cout << " quel" << std::endl;
  /*
    %%% 
    % Y = X*B+U  This is a standard Linear model:
    %     Y is the dependent variable, such as caudate displacement (nXd)
    %     X is the independent variables (nXp)
    %     B is the parameter matrix  (pXd)
    %     
    % A*B = C     This is the null hypothesis to test against
    %     A describes the variables that should be investigated with the rest
    %       'factored out' (qXp)
    %     B is the parameter matrix  (pXd)
    %     C completes the equation. It is often filled with zeros if the
    %       null-Hypothesis is that there is no correlation (qXd)
  */
  
  // computes the local pValues for each feature according to non-parametric permutation tests
//
// numSubjects = Number of Subjects Group A + Group B
// numFeatures = Number of Scalar features per Subject
// numPerms    = Number of Permutations used in computation
// tupelSize   = Numbers of features grouped in tupels 
// groupLabel = int array [numSubjects] with group Labels, Group A = 0, Group B = 1
// featureValue = double array [numSubjects*numFeatures] with scalar feature values for
//   all Subjects e.g. n-th feature of m-th subjects is at featureValue[m * numFeatures + n]
// significanceLevel = level of significance (usually 0.05 or 0.01)
// significanceSteps = number of significance steps

  PARSE_ARGS;
  std::cout << std::endl;
  if (infile.empty()) {
    std::cout << " no input file specified " << std::endl;exit(1);
  }
  else {
    std::cout << " input file specified: " << infile<< std::endl;
  }
  //if (outbase.empty()) {
    outbase = infile;

  //if the input file is a csv file
  char *extension=strchr(outbase.c_str(),'.');
  if(!strcmp(extension,".csv"))
  {
     int index=outbase.find(".csv",0);
     //erase the extension
     outbase=outbase.erase(index,outbase.size());
  }

  //}   ..inde column (in interaction test) or list of group columns (in group test), find match
  //testColumn = is index of match in independentColumns or groupTypeColumns
//TODO
if(interactionTest)
	{for( int i =0; i<numIndependent;i++)
		{ 
std::cout << " vect indecol " <<i <<" "<<independentColumns[i]<< std::endl;
if(independentColumns[i]==testColumn) { testColumn = i; }} }
else{
	for( int i =0; i<numGroupTypes;i++)
		{ if(groupTypeColumns[i]==testColumn) { testColumn = i; }} }
  
  unsigned int numSubjects, numFeatures;
  vnl_matrix<int> * groupLabel = NULL;
  vnl_matrix<double> * featureValue = NULL; 
  vnl_matrix<double> * scaleFactor = NULL;
  vnl_matrix<double> * indValue = NULL;
  std::string* meshFileNames = NULL;
  unsigned int numA;  // number of group A members
  unsigned int numB;  // number of group B members

  // load the mesh list file that defines all the inputs
  
  load_MeshList_file( infile,infileColumn, numIndependent, numGroupTypes,groupTypeColumns,independentColumns, surfListScale, scaleColumn, numSubjects, numA, numB, groupLabel, scaleFactor, meshFileNames, indValue, computeScaleFactorFromVolumes );


//TODO was here




  
  // load the actual mesh information
  MeshType::Pointer surfaceMesh = MeshType::New();
  MeshSpatialObjectType::Pointer  SOMesh;

  if (KWMreadableInputFile==0)
    load_Meshes( surfListScale, numSubjects, numIndependent, scaleFactor, indValue, meshFileNames, numFeatures, featureValue, surfaceMesh, SOMesh );
  //This part added by bp2009 to include the longitudinal analysis in the study...
  else
    load_KWMreadableInputFile( surfListScale, numSubjects, numIndependent, scaleFactor, indValue, meshFileNames, numFeatures, featureValue, surfaceMesh, SOMesh );
  //bp2009
  // compute surface properties

  PointsContainerPointer meanPoints  = PointsContainer::New();
  PointsContainerPointer meanPointsA = PointsContainer::New();
  PointsContainerPointer meanPointsB = PointsContainer::New();
  PointsContainerPointer meanPointsCorrected = PointsContainer::New();

  vnl_matrix<double> diffVectors(numFeatures,dimension);
  vnl_vector<double> normProjections(numFeatures);
  vnl_vector<double> normDistProjections(numFeatures);

  vnl_matrix<double> zScores(numFeatures,numSubjects);
  vnl_matrix<double> zScoresProjected(numFeatures,numSubjects);

  vtkPolyDataNormals *meanSurfaceNormals = vtkPolyDataNormals::New();
  vtkPolyDataNormals *meanSurfaceNormalsCorrected = vtkPolyDataNormals::New();


  // keep track of time

  time_t start,end;

  time (&start);

  compute_SurfaceProperties( meanPoints, meanPointsA, meanPointsB, diffVectors, normProjections, normDistProjections, zScores, zScoresProjected, numSubjects, numA, numB, numFeatures, groupLabel, featureValue, meanSurfaceNormals, surfaceMesh, outbase, meshFileNames , KWMreadableInputFile);

  // write out additional mesh information (means, differences, etc.)

  write_SurfaceProperties( outbase, meanPoints, meanPointsA, meanPointsB, diffVectors, normProjections, normDistProjections, zScoresProjected, zScores, writeZScores, surfaceMesh, SOMesh, meshFileNames, KWMreadableInputFile );


  // Now do the actual statistical testing based on Dimitrio's MANCOVA theory
  
  // select statistic

  int computeStatisticsType = -1;

  if ( useRoy )
    {
    computeStatisticsType = ROY;
    std::cout << "Using Roy statistics." << std::endl;
    }
  else if ( useHotelling )
    {
    computeStatisticsType = HOTELLING;
    std::cout << "Using Hotelling statistics." << std::endl;
    }
  else if ( usePillai )
    {
    computeStatisticsType = PILLAI;
    std::cout << "Using Pillai statistics." << std::endl;
    }
  else if ( useWilks )
    {
    computeStatisticsType = WILKS;
    std::cout << "Using Wilks statistics." << std::endl;
    }

  if ( computeStatisticsType<0 )
    {
    std::cerr << "No statistic type selected, defaulting to Hotelling." << std::endl;
    computeStatisticsType = HOTELLING;
    }

  std::cout << "Now running nonparametric tests:\n";

  vnl_vector<double> mancovaRawP;
  doMANCOVATesting( numSubjects, numFeatures, numGroupTypes, numIndependent, testColumn, numPerms, groupLabel, featureValue, interactionTest, numA, numB, significanceLevel, computeStatisticsType, mancovaRawP, meanPointsCorrected );

  // now compute the FDR corrected versions

  double fdrThresh;
  vnl_vector<double> fdrP;
  vnl_vector<double> bonferroniP;

  fdrP = fdrCorrection( mancovaRawP, significanceLevel, fdrThresh );
  std::cout << "fdr thresh (MANCOVA) is = " << fdrThresh << std::endl;

  bonferroniP = bonferroniCorrection( mancovaRawP );

  output_vector( mancovaRawP, outbase + std::string("_mancovaRawP"), ".txt" );
  output_vector( fdrP, outbase + std::string("_mancovaFDRP"), ".txt" );
  output_vector( bonferroniP, outbase + std::string("_mancovaBonferroniP"), ".txt" );



  // Write out the corrected mean mesh:
if (KWMreadableInputFile==0)
{
  surfaceMesh->SetPoints(meanPointsCorrected); 
  SOMesh->SetMesh(surfaceMesh);
  MeshWriterType::Pointer writer = MeshWriterType::New();
  writer->SetInput(SOMesh);
  
  std::string FilenameCorrectedAverage(outbase);
  FilenameCorrectedAverage += std::string("_GLM_correctedMean.meta");

  writer->SetFileName( FilenameCorrectedAverage.c_str() );
  writer->Update();
} //end if (KWMreadableInputFile==0)
  // compute simple interaction tests if desired

  if ( simpleCorrs )
    {

    if ( interactionTest ) 
      {

      int correlationType = -1;

      if ( negativeCorrelation )
	{
	correlationType = NEGCORR;
	}
      else if ( positiveCorrelation )
	{
	correlationType = POSCORR;
	}
      else if ( trendCorrelation )
	{
	correlationType = TCORR;
	}

      if ( correlationType<0 )
	{
	std::cerr << "No correlation type selected, defaulting to trend (=two-tail)." << std::endl;
	correlationType = TCORR;
	}
      
      vnl_vector<double> spearmanRhoDist(numFeatures);
      vnl_vector<double> spearmanRhoPro(numFeatures);
      
      vnl_vector<double> spearmanRhoProPval(numFeatures,0);// Calculated through permutations
      vnl_vector<double> spearmanRhoDistPval(numFeatures,0);
      
      vnl_vector<double> pearsonRhoDist(numFeatures);
      vnl_vector<double> pearsonRhoPro(numFeatures);
      
      vnl_vector<double> pearsonRhoProPval(numFeatures,0);// Calculated through permutations
      vnl_vector<double> pearsonRhoDistPval(numFeatures,0);
      
      /*if ( simpleCorrsCorrectedMean ) // needs to be fixed, if we do
	testing that we would need to correct every individual
	measurement as well, otherwise there will be some form of inconsistency
	{
	std::cout << "Using GLM corrected mean for the simple correlation test." << std::endl;

	// need to compute a new surface normal

	surfaceMesh->SetPoints(meanPointsCorrected); 

	itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
	convertMeshToVTK->SetInput(surfaceMesh);
	vtkPolyData * vtkMesh = convertMeshToVTK->GetOutput();
	
	meanSurfaceNormalsCorrected->SetComputePointNormals(1);
	meanSurfaceNormalsCorrected->SetComputeCellNormals(0);
	meanSurfaceNormalsCorrected->SetSplitting(0);
	meanSurfaceNormalsCorrected->AutoOrientNormalsOn();
	meanSurfaceNormalsCorrected->ConsistencyOn();
	meanSurfaceNormalsCorrected->FlipNormalsOff();  // all normals are outward pointing
	meanSurfaceNormalsCorrected->SetInput(vtkMesh);
	meanSurfaceNormalsCorrected->Update();

	do_ScalarInteractionTest( numSubjects, numFeatures, testColumn, featureValue, meanPointsCorrected, meanSurfaceNormalsCorrected, spearmanRhoDist, spearmanRhoPro, spearmanRhoProPval, spearmanRhoDistPval, pearsonRhoDist, pearsonRhoPro, pearsonRhoProPval, pearsonRhoDistPval, correlationType, computeParametricP, numPerms );

	}
      else
	{*/
	std::cout << "Using uncorrected mean for the simple correlation test." << std::endl;

	int writeOutVariabilityScores = 1; 
	do_ScalarInteractionTest( numSubjects, numFeatures, testColumn, featureValue, meanPoints, meanSurfaceNormals, spearmanRhoDist, spearmanRhoPro, spearmanRhoProPval, spearmanRhoDistPval, pearsonRhoDist, pearsonRhoPro, pearsonRhoProPval, pearsonRhoDistPval, correlationType, computeParametricP, numPerms, outbase );
	//}      

      // now compute the FDR and the Bonferroni corrected versions
      
      double fdrThresh;
      vnl_vector<double> fdrP;
      vnl_vector<double> bonferroniP;

      bonferroniP = bonferroniCorrection( pearsonRhoDistPval );
      fdrP = fdrCorrection( pearsonRhoDistPval, significanceLevel, fdrThresh );
      std::cout << "fdr thresh (Pearson distance) is = " << fdrThresh << std::endl;
      
      output_vector(bonferroniP, outbase, std::string("_normDistProjectionsPearsonPvalBonferroni.txt"));
      output_vector(fdrP, outbase, std::string("_normDistProjectionsPearsonPvalFDR.txt"));
      
      bonferroniP = bonferroniCorrection( pearsonRhoProPval );
      fdrP = fdrCorrection( pearsonRhoProPval, significanceLevel, fdrThresh );
      std::cout << "fdr thresh (Pearson projection) is = " << fdrThresh << std::endl;
      
      output_vector(bonferroniP, outbase, std::string("_normProjectionsPearsonPvalBonferroni.txt"));
      output_vector(fdrP, outbase, std::string("_normProjectionsPearsonPvalFDR.txt"));

      bonferroniP = bonferroniCorrection( spearmanRhoDistPval );
      fdrP = fdrCorrection( spearmanRhoDistPval, significanceLevel, fdrThresh );
      std::cout << "fdr thresh (Spearman distance) is = " << fdrThresh << std::endl;
      
      output_vector(bonferroniP, outbase, std::string("_normDistProjectionsSpearmanPvalBonferroni.txt"));
      output_vector(fdrP, outbase, std::string("_normDistProjectionsSpearmanPvalFDR.txt"));
      
      bonferroniP = bonferroniCorrection( spearmanRhoProPval );
      fdrP = fdrCorrection( spearmanRhoProPval, significanceLevel, fdrThresh );
      std::cout << "fdr thresh (Spearman projection) is = " << fdrThresh << std::endl;
      
      output_vector(bonferroniP, outbase, std::string("_normProjectionsSpearmanPvalBonferroni.txt"));
      output_vector(fdrP, outbase, std::string("_normProjectionsSpearmanPvalFDR.txt"));
      
      // write the simple correlation test results
      
      output_vector(spearmanRhoPro, outbase, std::string("_normProjectionsSpearman.txt"));
      output_vector(spearmanRhoDist, outbase, std::string("_normDistProjectionsSpearman.txt"));
      
      output_vector(spearmanRhoProPval, outbase, std::string("_normProjectionsSpearmanPval.txt"));
      output_vector(spearmanRhoDistPval, outbase, std::string("_normDistProjectionsSpearmanPval.txt"));
      
      output_vector(pearsonRhoPro, outbase, std::string("_normProjectionsPearson.txt"));
      output_vector(pearsonRhoDist, outbase, std::string("_normDistProjectionsPearson.txt"));
      
      output_vector(pearsonRhoProPval, outbase, std::string("_normProjectionsPearsonPval.txt"));
      output_vector(pearsonRhoDistPval, outbase, std::string("_normDistProjectionsPearsonPval.txt"));
      
      }

    }

  time (&end);
  double dif = round(difftime (end,start));
  unsigned int iif = (unsigned int) dif;

  std::cout << "Runtime was "<< iif/60 << " minute(s) and " << iif%60 << " second(s)." << std::endl;

  // clean up everything
if (KWMreadableInputFile==0)
{

  meanSurfaceNormals->Delete();
  meanSurfaceNormalsCorrected->Delete();
  surfaceMesh->Delete();
  SOMesh->Delete();
}

  write_ColorMap(outbase,interactionTest);

  write_MRMLScene(outbase,interactionTest);

  if ( groupLabel!=NULL ) delete groupLabel;
  if ( featureValue!=NULL ) delete featureValue;
  if ( scaleFactor!=NULL ) delete scaleFactor;
  if ( indValue!=NULL ) delete indValue;
  
  if ( meshFileNames!=NULL ) delete [] meshFileNames;

  // done cleaning up

  return 0;
}



