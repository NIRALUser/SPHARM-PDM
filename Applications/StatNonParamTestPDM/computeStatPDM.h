#ifndef COMPUTE_PERM_P_VAL_H__MST
#define COMPUTE_PERM_P_VAL_H__MST

// #include <itkMeanCalculator.h>
#include <itkMeanSampleFilter.h>

// #include <itkCovarianceCalculator.h>
#include <itkCovarianceSampleFilter.h>

#include <itkArray.h>
#include <itkMatrix.h>
#include <itkListSample.h>
#include <itkSymmetricEigenAnalysis.h>

#include "createPerm.h"

// rewrite this to be all included in a single class

enum { MeasurementVectorSize = 3 };
typedef float                                                                      MeasurementType;
typedef itk::Vector<MeasurementType, MeasurementVectorSize>                        MeasurementVectorType;
typedef itk::Matrix<MeasurementType, MeasurementVectorSize, MeasurementVectorSize> MeasurementMatrixType;
typedef itk::Statistics::ListSample<MeasurementVectorType>                         SampleType;

// typedef itk::Statistics::MeanCalculator< SampleType > MeanCalculatorType;
typedef itk::Statistics::MeanSampleFilter<SampleType> MeanSampleFilterType;

// typedef itk::Statistics::CovarianceCalculator< SampleType > CovarianceCalculatorType;
typedef itk::Statistics::CovarianceSampleFilter<SampleType> CovarianceSampleFilterType;

typedef itk::SymmetricEigenAnalysis<MeasurementMatrixType, MeasurementVectorType> EigenSystemType;

typedef itk::Vector<MeasurementType, 1>                          PvalueMeasurementVectorType;
typedef itk::Statistics::ListSample<PvalueMeasurementVectorType> PvalueSampleType;

typedef struct
  {
  double statVal;
  int index;
  } StatElement;

int doTesting(int numSubjects, int numFeatures, int numPerms, int tupelSize, int * groupLabel, double * featureValue,
              double significanceLevel, double FDRdiscoveryLevel, int significanceSteps, double * pValue,
              char *outbase);

// computes the local pValues for each feature according to non-parametric permutation tests, considering p-values for
// equalizing sensitivity

void computeGroupVecHotelDiffStat(int numSubjects, int numFeatures, int tupelSize, int * groupLabel,
                                  double * featureValue, double * statistic, double & avgStat,
                                  double & quant95Stat);

// computes the group difference statistic for the given group Labels
// computes Hotelling metric using 3D coordinates

void computeGroupVecMeanDiffStat(int numSubjects, int numFeatures, int tupelSize, int * groupLabel,
                                 double * featureValue, double * statistic, double & avgStat,
                                 double & quant95Stat);

// computes the group difference statistic for the given group Labels
// computes Absolute Mean distance metric using 3D coordinates

void computeGroupDiffStat(int numSubjects, int numFeatures, int * groupLabel, double * featureValue, double * statistic,
                          double & avgStat,
                          double & quant95Stat);

// computes the group difference statistic for the given group Labels
// computes absolute mean difference

double computeQuantile(int numObs, double * stat, double quantile);

// computes the value in the provided statistic for the given quantile value
// numObs = Number of Observations
// stat = double array[numObs ] contains the test statisticss

double computePValue(int numPerms, double * stat, double statVal);

// computes the Pval for a single permutation statistics at the given Value
// single right tail (maximum) p-val
// same as given the value, what is the quantile %

void computePermStatPval(int numFeatures, int numPerms, double * permStat, double * permStatPval);

// computes the Pval for all permutation statistics at all values
// single right tail (maximum) p-val
// same as given the values, what are their respective quantile %

int smallerStatElem(StatElement * elem1, StatElement * elem2);

// comparison function for sorting

#endif
