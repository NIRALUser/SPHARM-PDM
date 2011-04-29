
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>

#include "computeStatPDM.h"

#include "itkSampleFalseDiscoveryRateCorrectionFilter.h"  




using namespace std;
 
const int debug = 1;
const int writePvalDEBUG = 1;

void computeGroupMeanDiffStat(int numSubjects, int numFeatures,
                  int * groupLabel, double * featureValue,
			      double * statistic, double &avgStat, double &quant95Stat)
// computes the group difference statistic for the given group Labels 
//
// the group difference statistic is computed as the absolute difference of the means for each
// feature
//
// numSubjects = Number of Subjects Group A + Group B
// numFeatures = Number of scalar features per Subject
// groupLabel = int array [numSubj2ects] with group Labels, Group A = 0, Group B = 1
// featureValue = double array [numSubjects*numFeatures] with scalar feature values for
//   all Subjects e.g. n-th feature of m-th subjects is at featureValue[m * numFeatures + n]
// statistic = double array [numFeatures] output with group difference statistics
//
//
{
    // compute mean for group A and group B for each feature
    double * meanA = new double [numFeatures];
    double * meanB = new double [numFeatures];
    int feat;
    for (feat = 0; feat < numFeatures; feat++) {
     meanA[feat] = 0;
     meanB[feat] = 0;
    }
    
    int numSubjA = 0;
    int numSubjB = 0;
    
    for (int subj = 0; subj < numSubjects; subj++) {
     int subjIndex = subj * numFeatures;
     if (groupLabel[subj] == GROUP_A_LABEL) {
         numSubjA++;
         for (feat = 0; feat < numFeatures; feat++) {
          meanA[feat] = meanA[feat] + featureValue[subjIndex + feat];
         }
     } else if (groupLabel[subj] == GROUP_B_LABEL) {
         numSubjB++;
         for (feat = 0; feat < numFeatures; feat++) {
          meanB[feat] = meanB[feat] + featureValue[subjIndex + feat];         
         }
     } else {
         cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
     }
    }
    for (feat = 0; feat < numFeatures; feat++) {
      meanA[feat] = meanA[feat] / numSubjA;
      meanB[feat] = meanB[feat] / numSubjB;
    }

    // absolute difference is statistic    
    for (feat = 0; feat < numFeatures; feat++) {
      double val = fabs(meanA[feat] - meanB[feat]);
      statistic[feat] = val;
    }
    
    // compute 95th quantile
    quant95Stat = computeQuantile(numFeatures, statistic, 0.95);

    // compute mean 
    avgStat = 0;
    for (feat = 0; feat < numFeatures; feat++) {
      avgStat += statistic[feat];
    }
    avgStat /= numFeatures;


    delete meanA;
    delete meanB;

}

void computeGroupStudentTDiffStat(int numSubjects, int numFeatures,
                 int * groupLabel, double * featureValue,
                 double * statistic, double &avgStat, double &quant95Stat)
// computes the group difference statistic for the given group Labels using the Student t statistic
//
// the group difference statistic is computed as the absolute difference of the means for each
// feature weighted with their variances (Student t statistic)
//
// numSubjects = Number of Subjects Group A + Group B
// numFeatures = Number of scalar features per Subject
// groupLabel = int array [numSubj2ects] with group Labels, Group A = 0, Group B = 1
// featureValue = double array [numSubjects*numFeatures] with scalar feature values for
//   all Subjects e.g. n-th feature of m-th subjects is at featureValue[m * numFeatures + n]
// statistic = double array [numFeatures] output with group difference statistics
//
//
{
    // compute mean for group A and group B for each feature
    double * meanA = new double [numFeatures];
    double * meanB = new double [numFeatures];
    double * varA = new double [numFeatures];
    double * varB = new double [numFeatures];

    int feat;
    for (feat = 0; feat < numFeatures; feat++) {
      meanA[feat] = 0;
      meanB[feat] = 0;
      varA[feat] = 0;
      varB[feat] = 0;
    }
    
    int numSubjA = 0;
    int numSubjB = 0;
    
    // compute means
    int subj;
    for (subj = 0; subj < numSubjects; subj++) {
     int subjIndex = subj * numFeatures;
     if (groupLabel[subj] == GROUP_A_LABEL) {
         numSubjA++;
         for (feat = 0; feat < numFeatures; feat++) {
          meanA[feat] = meanA[feat] + featureValue[subjIndex + feat];
         }
     } else if (groupLabel[subj] == GROUP_B_LABEL) {
         numSubjB++;
         for (feat = 0; feat < numFeatures; feat++) {
          meanB[feat] = meanB[feat] + featureValue[subjIndex + feat];         
         }
     } else {
         cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
     }
    }
    for (feat = 0; feat < numFeatures; feat++) {
     meanA[feat] = meanA[feat] / numSubjA;
     meanB[feat] = meanB[feat] / numSubjB;
    }
    // compute variances
    for (subj = 0; subj < numSubjects; subj++) {
     int subjIndex = subj * numFeatures;
     double diff;
     if (groupLabel[subj] == GROUP_A_LABEL) {
         for (feat = 0; feat < numFeatures; feat++) {
                         diff = featureValue[subjIndex + feat] - meanA[feat];
                         varA[feat] = varA[feat] + diff * diff;
         }
     } else if (groupLabel[subj] == GROUP_B_LABEL) {
         for (feat = 0; feat < numFeatures; feat++) {
                         diff = featureValue[subjIndex + feat] - meanB[feat];
                         varB[feat] = varB[feat] + diff*diff;         
         }
     } else {
         cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
     }
    }
    // varA = sigma_a^2 * (n_a - 1)
    // varB = sigma_b^2 * (n_b - 1)

    // absolute difference is statistic
    for (feat = 0; feat < numFeatures; feat++) {
      double val = fabs(meanA[feat] - meanB[feat]) / sqrt( (varA[feat] + varB[feat]) * ( numSubjA  + numSubjB) 
                        / numSubjA / numSubjB / (numSubjA + numSubjB - 2)) ;
      statistic[feat] = val;
    }

    // compute 95th quantile
    quant95Stat = computeQuantile(numFeatures, statistic, 0.95);

    // compute mean
    avgStat = 0;
    for (feat = 0; feat < numFeatures; feat++) {
      avgStat += statistic[feat];
    }
    avgStat /= numFeatures;

    delete meanA;
    delete meanB;
    delete varA;
    delete varB;

}

int doTesting(int numSubjects, int numFeatures, int numPerms, int tupelSize,
           int * groupLabel, double * featureValue,
           double significanceLevel, double FDRdiscoveryLevel,int significanceSteps,
           double * pValue, char *outbase)
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
// pValue = double array [numFeatures] output with p-values at each feature location
//          the pValue is == 1 for all features with pValue greater than significanceLevel
//          for features with pValue less than significanceLevel, the actual p-values are 
//          computed in step sizes of significanceLevel/significanceSteps; 
//  e.g. significanceValue = 0.05, significanceSteps = 10; so, the possible p-values are:
//    1 (non-significant), 0.05, 0.045, 0.040, 0.035, 0.030, 0.025, 0.020, 0.015, 0.010, 0.005
//
{
    int numSubjA, numSubjB;
    numSubjA = 0; numSubjB = 0;
    for (int subj = 0; subj < numSubjects; subj++) {
      if (groupLabel[subj] == GROUP_A_LABEL) {
        numSubjA++;
      } else if (groupLabel[subj] == GROUP_B_LABEL) {
        numSubjB++;         
      } else {
        cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
      }
    }
    std::cout << "group sizes: " << numSubjA << "/" << numSubjB << " = " << numSubjects << std::endl;
    std::cout << "Features : " << numFeatures << " , real: " << numFeatures / tupelSize << std::endl;

    double * permStat = new double [numPerms * numFeatures / tupelSize];
    int * newGroupLabel = new int [numSubjects]; 
    double * realStudTStat = new double [numFeatures / tupelSize];
    double * realMeanDiffStat = new double [numFeatures / tupelSize];
    double * avgStat = new double [numPerms];
    double * quant95Stat = new double [numPerms];

    double * realStat = permStat; // the first permutation is the real one
    int updateInc = numPerms/40+1;

    int perm;
    // generate the permutations and the statistics of the permutations
    for (perm = 0; perm < numPerms; perm++) {
      if (perm % updateInc == 0) {
       std::cout << "."; std::cout.flush();
      }

      if (perm == 0) {  // the first permutation is the real one 
	for (int i=0; i<numSubjects; i++) {newGroupLabel[i]=groupLabel[i];}
      } else {
	generatePermGroup(groupLabel, numSubjA, numSubjB, newGroupLabel );
      }

      if (tupelSize == 1) {
        computeGroupMeanDiffStat(numSubjects, numFeatures,  
                                 newGroupLabel, featureValue, 
				 &(permStat[perm * numFeatures / tupelSize]), avgStat[perm], quant95Stat[perm]);
        if (perm == 0) {
          double avgTDiff,quant95TDiff;
          computeGroupStudentTDiffStat(numSubjects, numFeatures,
                                       newGroupLabel, featureValue,realStudTStat, avgTDiff,quant95TDiff);
        }
      } else  {
        computeGroupVecHotelDiffStat(numSubjects, numFeatures,  tupelSize,
                                     newGroupLabel, featureValue, 
                                     &(permStat[perm * numFeatures / tupelSize]), avgStat[perm], quant95Stat[perm]); 
        if (perm == 0) {
          double avgMeanDiff,quant95Diff;
          computeGroupVecMeanDiffStat(numSubjects, numFeatures,  tupelSize,
                                      newGroupLabel, featureValue, realMeanDiffStat, avgMeanDiff, quant95Diff); 
        }
      }
      //std::cout << avgStat[perm] << std::endl;
    }
    std::cout << std::endl;

    // compute p-values
    double * permStatPval = new double [numPerms * numFeatures / tupelSize];

    computePermStatPval(numFeatures / tupelSize, numPerms, permStat, permStatPval);

    double avgPval =  computePValue(numPerms, avgStat, avgStat[0]);
    double quant95Pval =  computePValue(numPerms, quant95Stat, quant95Stat[0]);
    
    std::cout << "writing avg shape diff pvalue " << avgPval << std::endl;
    // write uncorrected p-values
    string avgPvalFile(outbase);
    avgPvalFile = avgPvalFile + string("_avgPval.txt");
    std::ofstream efile;
    efile.open(avgPvalFile.c_str(), ios::out); 
    efile << "Pval of average statistic across surface:" << avgPval << std::endl;
    efile << "Pval of 95 percentile statistic across surface:" << quant95Pval << std::endl;
    efile.close();
    
    std::cout << "writing raw pvalues " << std::endl;
    // write uncorrected p-values
    string rawPvalFile(outbase);
    rawPvalFile = rawPvalFile + string("_rawPval.txt");
    std::ofstream e3file;
    e3file.open(rawPvalFile.c_str(), ios::out); 
    e3file << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
    e3file << "DIMENSION = 1" << std::endl;
    e3file << "TYPE = Scalar" << std::endl;
    int feat;
    for (feat = 0; feat < numFeatures / tupelSize; feat++)
      e3file << permStatPval[feat] << " ";
    e3file << std::endl;
    e3file.close();
    
    std::cout << "writing statistic (T^2)" << std::endl;
    // write uncorrected p-values
    string realStatFile(outbase);
    realStatFile = realStatFile + string("_realStat.txt");
    std::ofstream e4file;
    e4file.open(realStatFile.c_str(), ios::out);
    e4file << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
    e4file << "DIMENSION = 1" << std::endl;
    e4file << "TYPE = Scalar" << std::endl; 

    for (feat = 0; feat < numFeatures / tupelSize; feat++)
      e4file << realStat[feat] << " ";
    e4file << std::endl;
    e4file.close();

    if (tupelSize == 1) {
      std::cout << "writing t statistic" << std::endl;
      // write uncorrected p-values
      string realTStatFile(outbase);
      realTStatFile = realTStatFile + string("_tstat.txt");
      std::ofstream e5file;
      e5file.open(realTStatFile.c_str(), ios::out); 
      e5file << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
      e5file << "DIMENSION = 1" << std::endl;
      e5file << "TYPE = Scalar" << std::endl;
      for (feat = 0; feat < numFeatures / tupelSize; feat++)
       e5file << realStudTStat[feat] << " ";
      e5file << std::endl;
      e5file.close();

      //       std::cout << "writing Cohen's d statistic" << std::endl;
      //       // write uncorrected p-values
      //       string realZStatFile(outbase);
      //       realZStatFile = realZStatFile + string("_Dstat.txt");
      //       efile.open(realZStatFile.c_str(), ios::out); 
      //       efile << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
      //       efile << "DIMENSION = 1" << std::endl;
      //       efile << "TYPE = Scalar" << std::endl;
      //       for (feat = 0; feat < numFeatures / tupelSize; feat++)
      //        efile << 2 * realStudTStat[feat]/sqrt((double) (numSubjA + numSubjB - 2)) << " ";
      //       efile << std::endl;
      //       efile.close();
      
      //       std::cout << "writing r-square statistic" << std::endl;
      //       // write uncorrected p-values
      //       string realRStatFile(outbase);
      //       realRStatFile = realRStatFile + string("_r2stat.txt");
      //       efile.open(realRStatFile.c_str(), ios::out); 
      //       efile << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
      //       efile << "DIMENSION = 1" << std::endl;
      //       efile << "TYPE = Scalar" << std::endl;
      //       for (feat = 0; feat < numFeatures / tupelSize; feat++) {
      //        double factor = realStudTStat[feat] * realStudTStat[feat] /(numSubjA + numSubjB - 2);
      //        efile << factor / (1 + factor) << " ";
      //       }
      //       efile << std::endl;
      //       efile.close();

    } else {
      // write T stats = sqrt (T^2)
      std::cout << "writing T statistic" << std::endl;
      string realTStatFile(outbase);
      realTStatFile = realTStatFile + string("_Tstat.txt");
      std::ofstream e5file;
      efile.open(realTStatFile.c_str(), ios::out); 
      e5file << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
      e5file << "DIMENSION = 1" << std::endl;
      e5file << "TYPE = Scalar" << std::endl;
      for (feat = 0; feat < numFeatures / tupelSize; feat++)
       e5file << sqrt(realStat[feat]) << " ";
      e5file << std::endl;
      e5file.close();

      //       // write r^2: it is not yet clear whether this is the correct formulation of r^2, as we are using 
      //       // a non-standard, robust T^2 formulation
      //       std::cout << "writing r-square statistic" << std::endl;
      //       string realRStatFile(outbase);
      //       realRStatFile = realRStatFile + string("_r2stat.txt");
      //       efile.open(realRStatFile.c_str(), ios::out); 
      //       efile << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
      //       efile << "DIMENSION = 1" << std::endl;
      //       efile << "TYPE = Scalar" << std::endl;
      //       for (feat = 0; feat < numFeatures / tupelSize; feat++) {
      //        double factor = realStat[feat]/(numSubjA + numSubjB - 2);
      //        efile << factor / (1 + factor) << " ";
      //       }
      //       efile << std::endl;
      //       efile.close();

      std::cout << "writing Mean Diff statistic" << endl;
      // write uncorrected effect size = difference of means
      string realMeanDiffStatFile(outbase);
      realMeanDiffStatFile = realMeanDiffStatFile + string("_MeanDiff.txt");
      std::ofstream e6file;
      e6file.open(realMeanDiffStatFile.c_str(), ios::out); 
      e6file << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
      e6file << "DIMENSION = 1" << std::endl;
      e6file << "TYPE = Scalar" << std::endl;
      for (feat = 0; feat < numFeatures / tupelSize; feat++)
       e6file << realMeanDiffStat[feat] << " ";
      e6file << std::endl;
      e6file.close();
    }

    // compute the correction via False Discovery Rate
    {

     typedef neurolib::Statistics::SampleFalseDiscoveryRateCorrectionFilter<PvalueSampleType> FDRFilterType;
      
      FDRFilterType::Pointer FDRFilter = FDRFilterType::New();
      FDRFilter->SetMaximumFalseDiscoveryRate(significanceLevel);
      FDRFilter->SetNumberOfSteps(significanceSteps);
  
      
      PvalueSampleType::Pointer pvalueSample = PvalueSampleType::New();  
      pvalueSample->SetMeasurementVectorSize( 1 );
      {
	FDRFilterType::MeasurementVectorType mv;
	for (int feat = 0; feat < numFeatures / tupelSize; feat++) 
	  {
	    mv[0] = permStatPval[feat];
	    pvalueSample->PushBack(mv);
	  }
      }
      std::cout << "computing FDR " << std::endl;
      //FDRFilter->SetInputSample(pvalueSample.GetPointer());
FDRFilter->SetInput(pvalueSample.GetPointer());
      FDRFilter->Update() ;
      
      pvalueSample = FDRFilter->GetOutput();         
      
      string outfileFDR(outbase);
      outfileFDR = outfileFDR + string("_FDRPval.txt");
      std::ofstream e2file;
      e2file.open(outfileFDR.c_str(), ios::out); 
      e2file << "NUMBER_OF_POINTS = " << pvalueSample->Size() << std::endl;
      e2file << "DIMENSION = 1" << std::endl;
      e2file << "TYPE = Scalar" << std::endl;
      PvalueSampleType::ConstIterator iterPval = pvalueSample->Begin();
      while (iterPval != pvalueSample->End())  
	{
	  e2file << (iterPval.GetMeasurementVector())[0] << " ";
	  ++iterPval;
	}
      e2file << std::endl;
      e2file.close();
    }

    // compute the correction via the distribution across all permutation of the minimum across the surface 
    // within all permutations

    std::cout << "computing Min distribution" << std::endl;
    
    // compute Minimum Statistic over all features for each permutation
    double * minStat = new double [numPerms]; 
    for (perm = 0; perm < numPerms; perm++) {
      double percLevel;
      percLevel = 0.0000001;
      double percValue = computeQuantile(numFeatures / tupelSize, 
					  &(permStatPval[perm * numFeatures  / tupelSize]), 
					  percLevel); 
      minStat[perm] = percValue;
    }

    

    // reset pValue
    for (feat = 0; feat < numFeatures / tupelSize; feat++) pValue[feat] = 1.0;
    
    // loop over all significance levels
    double significanceStepsize = FDRdiscoveryLevel / significanceSteps;
    for (double curSignLevel = FDRdiscoveryLevel; curSignLevel >= 0.00001;
      curSignLevel = curSignLevel - significanceStepsize) {

     // compute significance value in minimum statistic at significance level 
     double curSignValue = computeQuantile(numPerms, minStat, curSignLevel);

     // find statistic treshold for each feature f using 
     // the significance value of the minimum statistic
     for (feat = 0; feat < numFeatures  / tupelSize; feat++) {

         // if real value is equal or greater then 
         // feature Threshold then significant
         //if (realStat[feat] >= featThresh) {
         if (permStatPval[feat] <= curSignValue) {
          pValue[feat] = curSignLevel;
         }
      }
    }
     

    string outfile(outbase);
    outfile = outfile + string("_corrPval.txt");
    std::ofstream e7file;
    e7file.open(outfile.c_str(), ios::out); 
    e7file << "NUMBER_OF_POINTS = " << numFeatures / tupelSize << std::endl;
    e7file << "DIMENSION = 1" << std::endl;
    e7file << "TYPE = Scalar" << std::endl;
    for (feat = 0; feat < numFeatures / tupelSize; feat++)
      e7file << pValue[feat] << " ";
    e7file << std::endl;
    e7file.close();

    delete permStat;
    delete realStudTStat;
    delete permStatPval;
    delete minStat;
    delete avgStat;
    delete newGroupLabel;
    delete realMeanDiffStat;

    return 1;
}

double computePValue(int numPerms, double * stat, double statVal)
// computes the p-value at the given statistical Value using the permutations 
// computes single right tail significance level (appropriate for Maximum statistic)
//
// numPerms = Number of Permutations
// stat = double array[numPerms ] contains the test statisticss
{ 
  static int first = 0;
    if (!first) {
     first = 1;
     srand(time(NULL));
    }

    StatElement * sortPermStat= new StatElement[numPerms];

    for (int perm = 0; perm < numPerms; perm++) {
     sortPermStat[perm].statVal = stat[perm];
     sortPermStat[perm].index = perm;
    }
    
    // sort, smallest first
    qsort(sortPermStat, numPerms, sizeof(StatElement), 
       (int (*) (const void *, const void *)) smallerStatElem);

    // find smallest element
    bool found = false;
    int sortIndex = -1;
    for (int pIndex = numPerms - 1 ; pIndex >= 0 && !found ; pIndex--) {
      if (statVal >= sortPermStat[pIndex].statVal ) {
       found = true;
       sortIndex = pIndex;
       // std::cout << statVal << "," << sortPermStat[pIndex].statVal << "," << pIndex << std::endl;
      }
    }
    double retval = 1.0 - (double) sortIndex/numPerms;
    
    delete sortPermStat;

    return retval;

}

double computeQuantile(int numObs, double * stat, double quantile)
// computes the value in the provided statistic for the given quantile value 
//
// numObs = Number of Observations
// stat = double array[numObs ] contains the test statisticss
{ 
  static int first = 0;
    if (!first) {
     first = 1;
     srand(time(NULL));
    }

    StatElement * sortStat= new StatElement[numObs];

    for (int perm = 0; perm < numObs; perm++) {
     sortStat[perm].statVal = stat[perm];
     sortStat[perm].index = perm;
    }
    
    // sort, smallest first
    qsort(sortStat, numObs, sizeof(StatElement), 
       (int (*) (const void *, const void *)) smallerStatElem);

    // index at value
    double quantindex = (double) numObs * quantile;
    if (quantile == 1.0 ) quantindex = numObs;
    
    double retval = stat[sortStat[(int) quantindex].index];
    
    delete sortStat;

    return retval;

}

void computePermStatPval(int numFeatures, int numPerms, 
                double * permStat, double * permStatPval)
// computes the Pval for all permutation statistics
// the p-val is computed as the percentual ordered rank over all permutations
//
// numFeatures = Number of scalar features per Subject
// numPerms = Number of Permutations
// permStat = double array[numPerms * numFeatures] contains the test statistics
// permStatPval = double array [numPerms * numFeatures] returns the p-val of the statistics
{
    static int first = 0;
    if (!first) {
     first = 1;
     srand(time(NULL));
    }

    int feat;
    int perm;

    StatElement * sortPermStat= new StatElement[numPerms];

    for (feat = 0; feat < numFeatures; feat++) {

     for (perm = 0; perm < numPerms; perm++) {
         sortPermStat[perm].statVal = permStat[perm * numFeatures + feat];
         sortPermStat[perm].index = perm;
     }
     
     // sort, smallest first
     qsort(sortPermStat, numPerms, sizeof(StatElement), 
           (int (*) (const void *, const void *)) smallerStatElem);
     
     double prevPval = 0;
     double curPval = 0;
     for (perm = 0 ; perm < numPerms; perm++) {

       // percentual rank 0..1 -> cumulative probability -> p-val
       double nextPval = 1.0 - (double) (perm+1) / (double) numPerms; 
       
       int curIndex = sortPermStat[perm].index;
       
       if ((perm == 0) || (sortPermStat[perm].statVal != sortPermStat[perm-1].statVal)) {
         // current value is different from previous value (or first value), 
         // thus step up p-value
         prevPval = curPval;
         curPval = nextPval;
       }

       permStatPval[curIndex * numFeatures + feat] = curPval;

     }
    }
    delete sortPermStat;
}

void ALTcomputePermStatPval(int numFeatures, int numPerms, 
                double * permStat, double * permStatPval)
// computes the Pval for all permutation statistics
// the p-val is computed as the percentual ordered rank over all permutations
//
// numFeatures = Number of scalar features per Subject
// numPerms = Number of Permutations
// permStat = double array[numPerms * numFeatures] contains the test statistics
// permStatPval = double array [numPerms * numFeatures] returns the p-val of the statistics
{
    static int first = 0;
    if (!first) {
     first = 1;
     srand(time(NULL));
    }

    StatElement * sortPermStat= new StatElement[numPerms];
    int feat;
    int perm;

    for (feat = 0; feat < numFeatures; feat++) {

     for (perm = 0; perm < numPerms; perm++) {
         sortPermStat[perm].statVal = permStat[perm * numFeatures + feat];
         sortPermStat[perm].index = perm;
     }
     
     // sort, smallest first
     qsort(sortPermStat, numPerms, sizeof(StatElement), 
           (int (*) (const void *, const void *)) smallerStatElem);
     
     double CurPval = 1;
     double stepSize = fabs(sortPermStat[(int) ((double) 0.96 * numPerms)].statVal -
                      sortPermStat[(int) ((double) 0.95 * numPerms)].statVal);
     if (stepSize == 0) 
       stepSize = fabs(sortPermStat[numPerms - 1].statVal - sortPermStat[0].statVal);
     double CurBinMin = sortPermStat[0].statVal;

     for (perm = 0 ; perm < numPerms; perm++) {
       if (sortPermStat[perm].statVal > CurBinMin + stepSize) {
         CurBinMin = sortPermStat[perm].statVal;
         CurPval = 1.0 - (double) (perm+1) / (double) numPerms;
       }
       
       int curIndex = sortPermStat[perm].index;
       
       permStatPval[curIndex * numFeatures + feat] = CurPval;
     }
    }

    delete sortPermStat;

}


int smallerStatElem(StatElement * elem1, StatElement * elem2)
// comparison function for sorting
{
    if (elem1->statVal > elem2->statVal) return 1;
    else if (elem1->statVal < elem2->statVal) return -1;
    else return 0;
}
