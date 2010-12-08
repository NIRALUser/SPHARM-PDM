#include "nonMANCOVA.h"

void do_ScalarInteractionTest( unsigned int numSubjects, unsigned int numFeatures, unsigned int testColumn, vnl_matrix<double>  * &featureValue,  PointsContainerPointer &meanPoints, vtkPolyDataNormals *& MeshNormals, vnl_vector<double>& spearmanRhoDist, vnl_vector<double>& spearmanRhoPro, vnl_vector<double>& spearmanRhoProPval, vnl_vector<double>& spearmanRhoDistPval, vnl_vector<double>& pearsonRhoDist, vnl_vector<double>& pearsonRhoPro, vnl_vector<double>& pearsonRhoProPval, vnl_vector<double>& pearsonRhoDistPval, int correlationType, bool useParametricP, unsigned int numPerms , std::string outbase)
{

  vtkDataArray *ArrayNormal = MeshNormals->GetOutput()->GetPointData()->GetNormals();

  vnl_vector<double> dataToCor = featureValue->get_column(numFeatures*tupelSize +testColumn);
  vnl_vector<double> dataRanked = TiedRankCalc(dataToCor);
  vnl_vector<double> currentDataRanked = dataRanked;
  vnl_vector<double> currentDataToCor = dataToCor;
  
  itk::RandomPermutation permutation(numSubjects);
  double temp;
  double *tempNorm; 

  vnl_vector<double> d,projsRanked,distsRanked;
  vnl_vector<double> projs(numSubjects,0.0),dists(numSubjects,0.0);

  std::cout << "Performing simple correlation testing ... ";
  
  if (true)// Correlate with normProjections
    {
    for (unsigned int feat=0; feat<numFeatures; ++feat)
      {

      tempNorm=ArrayNormal->GetTuple(feat);
      PointType currentMeanPoint =  meanPoints->GetElement(feat);

      for (unsigned int sub=0; sub<numSubjects; ++sub)
	{ 
	projs(sub) = 0.0;
	dists(sub) = 0.0;
	for (unsigned int tup=0; tup<tupelSize; ++tup)
	  {
	  // compute the projection along the mean normal
	  projs(sub)+= ((*featureValue)[sub][feat*tupelSize+tup]-currentMeanPoint[tup])*tempNorm[tup];
	  // compute the squared distance to the mean shape
	  dists(sub)+= ((*featureValue)[sub][feat*tupelSize+tup]-currentMeanPoint[tup])*((*featureValue)[sub][feat*tupelSize+tup]-currentMeanPoint[tup]);	
	  }
	if (projs(sub)<0)// Negative, means pointing in
	  dists(sub)=-sqrt(dists(sub));
	else
	  dists(sub)=sqrt(dists(sub));
	}

      projsRanked=TiedRankCalc(projs);
      distsRanked=TiedRankCalc(dists);


      output_vector(projs, outbase, std::string("_normProjectionsVariabilityScore.txt"));
      output_vector(dists, outbase, std::string("_normDistProjectionsVariabilityScore.txt"));
      
      // compute the Pearson correlation coefficient

      // parametric p-values

      double dPSP_diff, dPSP_positive, dPSP_negative;
      double dPSD_diff, dPSD_positive, dPSD_negative;
      double dPPP_diff, dPPP_positive, dPPP_negative;
      double dPPD_diff, dPPD_positive, dPPD_negative;

      spearmanRhoPro(feat) = computePearsonCorrelationWithP( dataRanked, projsRanked, dPSP_diff, dPSP_positive, dPSP_negative );
      spearmanRhoDist(feat) = computePearsonCorrelationWithP( dataRanked, distsRanked, dPSD_diff, dPSD_positive, dPSD_negative );
      pearsonRhoPro(feat) = computePearsonCorrelationWithP( dataToCor, projs, dPPP_diff, dPPP_positive, dPPP_negative );
      pearsonRhoDist(feat) = computePearsonCorrelationWithP( dataToCor, dists, dPPD_diff, dPPD_positive, dPPD_negative );

      spearmanRhoProPval(feat) = 0.0;
      spearmanRhoDistPval(feat) = 0.0;
      pearsonRhoProPval(feat) = 0.0;
      pearsonRhoDistPval(feat) = 0.0;

      if ( useParametricP )
	{
	if ( correlationType==NEGCORR )
	  {
	  spearmanRhoProPval(feat) = dPSP_negative;
	  spearmanRhoDistPval(feat) = dPSD_negative;
	  pearsonRhoProPval(feat) = dPPP_negative;
	  pearsonRhoDistPval(feat) = dPPD_negative;
	  }
	else if ( correlationType==POSCORR )
	  {
	  spearmanRhoProPval(feat) = dPSP_positive;
	  spearmanRhoDistPval(feat) = dPSD_positive;
	  pearsonRhoProPval(feat) = dPPP_positive;
	  pearsonRhoDistPval(feat) = dPPD_positive;
	  }
	else if ( correlationType==TCORR )
	  {
	  spearmanRhoProPval(feat) = dPSP_diff;
	  spearmanRhoDistPval(feat) = dPSD_diff;
	  pearsonRhoProPval(feat) = dPPP_diff;
	  pearsonRhoDistPval(feat) = dPPD_diff;
	  }
	else
	  {
	  std::cerr << "ERROR: Unknown correlation type " << correlationType << ". ABORT." << std::endl;
	  exit(-1);
	  }

	}
      else
	{

	// do permutation testing to compute the p-values

	for (unsigned int perm=0;perm<numPerms;++perm) 
	  { 
	  
	  permutation.Shuffle();
	  for(unsigned int sub=0;sub<numSubjects;++sub)
	    {
	    currentDataRanked(sub) = dataRanked(permutation[sub]);
	    currentDataToCor(sub) = dataToCor(permutation[sub]);
	    }
	  
	  // now compute the test statistics
	  
	  temp = computePearsonCorrelation( currentDataRanked, projsRanked );
	  
	  if (temp>spearmanRhoPro(feat))
	    spearmanRhoProPval(feat)=spearmanRhoProPval(feat)+1;
	  
	  temp= computePearsonCorrelation( currentDataRanked, distsRanked );
	  
	  if (temp>spearmanRhoDist(feat))
	    spearmanRhoDistPval(feat)=spearmanRhoDistPval(feat)+1;
	  
	  temp = computePearsonCorrelation( currentDataToCor, projs );
	  
	  if (temp>pearsonRhoPro(feat))
	    pearsonRhoProPval(feat)=pearsonRhoProPval(feat)+1;
	  
	  temp= computePearsonCorrelation( currentDataToCor, dists );
	  
	  if (temp>pearsonRhoDist(feat))
	    pearsonRhoDistPval(feat)=pearsonRhoDistPval(feat)+1;
	  
	  
	  } // end of permutation loop
	
	spearmanRhoProPval(feat)=spearmanRhoProPval(feat)/numPerms;
	spearmanRhoDistPval(feat)=spearmanRhoDistPval(feat)/numPerms;
	pearsonRhoProPval(feat)=pearsonRhoProPval(feat)/numPerms;
	pearsonRhoDistPval(feat)=pearsonRhoDistPval(feat)/numPerms;
	
	}
      }
      
    }
  //output_vector(spearmanRhoPro, std::string("./"), std::string("debug.txt"));
  std::cout << "done." << std::endl;

}


