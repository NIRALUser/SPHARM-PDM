
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <string>
#include <cmath>

#include "computeStatPDM.h"

void computeGroupVecMeanDiffStat(int numSubjects, int numFeatures, int tupelSize,
                                 int * groupLabel, double * featureValue,
                                 double * statistic, double & avgStat, double & quant95Stat)
// computes the group difference statistic for the given group Labels
//
// the group difference statistic is computed as the simple difference magnitude of the means
// numSubjects = Number of Subjects Group A + Group B
// numFeatures = Number of scalar features per Subject
// groupLabel = int array [numSubj2ects] with group Labels, Group A = 0, Group B = 1
// featureValue = double array [numSubjects*numFeatures] with scalar feature values for
//   all Subjects e.g. n-th feature of m-th subjects is at featureValue[m * numFeatures + n]
// statistic = double array [numFeatures/tupelSize] output with group difference statistics
// features are always gathered in n-tupels, here representing the x, y, z coordinate
//   (but can be used for any other n tupel)
// easy to extend to n-tupel
{

  static int curNumFeatures = numFeatures;

  static double * meanA = new double[numFeatures];
  static double * meanB = new double[numFeatures];

  if( ( curNumFeatures != numFeatures) )
    {
    delete meanA;
    delete meanB;
    meanA = new double[numFeatures];
    meanB = new double[numFeatures];
    }

  int numSubjA = 0;
  int numSubjB = 0;

  // compute mean for group A and group B for each feature
  // independent of tupel-size
  for( int feat = 0; feat < numFeatures; feat++ )
    {
    meanA[feat] = 0;
    meanB[feat] = 0;
    }
  for( int subj = 0; subj < numSubjects; subj++ )
    {
    int subjIndex = subj * numFeatures;
    if( groupLabel[subj] == GROUP_A_LABEL )
      {
      numSubjA++;
      for( int feat = 0; feat < numFeatures; feat++ )
        {
        meanA[feat] = meanA[feat] + featureValue[subjIndex + feat];
        }
      }
    else if( groupLabel[subj] == GROUP_B_LABEL )
      {
      numSubjB++;
      for( int feat = 0; feat < numFeatures; feat++ )
        {
        meanB[feat] = meanB[feat] + featureValue[subjIndex + feat];
        }
      }
    else
      {
      std::cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
      }
    }
  for( int feat = 0; feat < numFeatures; feat++ )
    {
    meanA[feat] = meanA[feat] / numSubjA;
    meanB[feat] = meanB[feat] / numSubjB;
    }

  double* diffmean = new double[tupelSize];
  // Mean distance  metric computation
  for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
    {
    for( int dim = 0; dim < tupelSize; dim++ )
      {
      diffmean[dim] = meanA[feat * tupelSize + dim] - meanB[feat * tupelSize + dim];
      }
    // (meanA - meanB) . (meanA - meanB)
    double dotVal = 0;
    for( int dim = 0; dim < tupelSize; dim++ )
      {
      dotVal = dotVal + diffmean[dim] * diffmean[dim];
      }
    statistic[feat] = sqrt(dotVal);
    }

  // compute 95th quantile
  quant95Stat = computeQuantile(numFeatures, statistic, 0.95);
  // compute mean
  avgStat = 0;
  for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
    {
    avgStat += statistic[feat];
    }
  avgStat /= numFeatures / tupelSize;

}

void computeGroupVecHotelDiffStat(int numSubjects, int numFeatures, int tupel_size,
                                  int * groupLabel, double * featureValue,
                                  double * statistic, double & avgStat, double & quant95Stat)
// computes the group difference statistic for the given group Labels
//
// the group difference statistic is computed as the Hotelling T^2 statistic
//     T^2 = (mean_A - mean_B)' Cov(A,B) (mean_A - mean_B)
//     Cov(A,B) = ( (numA - 1) Cov(A) + (numB - 1) Cov(B) ) / (numA + numB - 2)
// numSubjects = Number of Subjects Group A + Group B
// numFeatures = Number of scalar features per Subject
// groupLabel = int array [numSubj2ects] with group Labels, Group A = 0, Group B = 1
// featureValue = double array [numSubjects*numFeatures] with scalar feature values for
//   all Subjects e.g. n-th feature of m-th subjects is at featureValue[m * numFeatures + n]
// statistic = double array [numFeatures/tupelSize] output with group difference statistics
// features are always gathered in n-tupels, here representing the x, y, z coordinate
//   (but can be used for any other n tupel)
// easy to extend to n-tupel
{

  static int      curNumFeatures = numFeatures;
  static double * meanA = new double[numFeatures];
  static double * meanB = new double[numFeatures];

  if( ( curNumFeatures != numFeatures) )
    {
    delete meanA;
    delete meanB;
    meanA = new double[numFeatures];
    meanB = new double[numFeatures];
    }

  int numSubjA = 0;
  int numSubjB = 0;
  // compute mean for group A and group B for each feature
  // independent of tupel-size
  for( int feat = 0; feat < numFeatures; feat++ )
    {
    meanA[feat] = 0;
    meanB[feat] = 0;
    }
  for( int subj = 0; subj < numSubjects; subj++ )
    {
    int subjIndex = subj * numFeatures;
    if( groupLabel[subj] == GROUP_A_LABEL )
      {
      numSubjA++;
      for( int feat = 0; feat < numFeatures; feat++ )
        {
        meanA[feat] = meanA[feat] + featureValue[subjIndex + feat];
        }
      }
    else if( groupLabel[subj] == GROUP_B_LABEL )
      {
      numSubjB++;
      for( int feat = 0; feat < numFeatures; feat++ )
        {
        meanB[feat] = meanB[feat] + featureValue[subjIndex + feat];
        }
      }
    else
      {
      std::cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
      }
    }
  for( int feat = 0; feat < numFeatures; feat++ )
    {
    meanA[feat] = meanA[feat] / numSubjA;
    meanB[feat] = meanB[feat] / numSubjB;
    }

  if( tupel_size == 3 )
    {
    const int tupelSize = 3; // needs this due to hardcoded size of covariance matrix
    typedef itk::Matrix<double, tupelSize, tupelSize> CovMatType;

    static CovMatType * covarA   = new CovMatType[numFeatures / tupelSize];
    static CovMatType * covarB   = new CovMatType[numFeatures / tupelSize];
    static CovMatType * covarAB  = new CovMatType[numFeatures / tupelSize];
    static CovMatType * InvcovarAB  = new CovMatType[numFeatures / tupelSize];
    if( ( curNumFeatures != numFeatures) )
      {
      delete covarA;
      delete covarB;
      delete covarAB;
      delete InvcovarAB;
      covarA   = new CovMatType[numFeatures / tupelSize];
      covarB   = new CovMatType[numFeatures / tupelSize];
      covarAB  = new CovMatType[numFeatures / tupelSize];
      InvcovarAB  = new CovMatType[numFeatures / tupelSize];
      }
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      covarA[feat].Fill(0);
      covarB[feat].Fill(0);
      }
    for( int subj = 0; subj < numSubjects; subj++ )
      {
      int subjIndex = subj * numFeatures;
      if( groupLabel[subj] == GROUP_A_LABEL )
        {
        for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
          {
          for( int row = 0; row < tupelSize; row++ )
            {
            for( int column = 0; column < tupelSize; column++ )
              {
              int indexDiff1 =  feat * tupelSize + row;
              int indexDiff2 =  feat * tupelSize + column;
              (covarA[feat])[row][column] =  (covarA[feat])[row][column]
                + (featureValue[subjIndex + indexDiff1] - meanA[indexDiff1])
                * (featureValue[subjIndex + indexDiff2] - meanA[indexDiff2]);

              }
            }
          }
        }
      else if( groupLabel[subj] == GROUP_B_LABEL )
        {
        for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
          {
          for( int row = 0; row < tupelSize; row++ )
            {
            for( int column = 0; column < tupelSize; column++ )
              {
              int indexDiff1 =  feat * tupelSize + row;
              int indexDiff2 =  feat * tupelSize + column;
              (covarB[feat])[row][column] =  (covarB[feat])[row][column]
                + (featureValue[subjIndex + indexDiff1] - meanB[indexDiff1])
                * (featureValue[subjIndex + indexDiff2] - meanB[indexDiff2]);
              }
            }
          }
        }
      else
        {
        std::cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
        }
      }

    // compute the pooled covariance matrices [nxn] covarAB,
    // matrix is fully symmetric -> does not matter whether column or row first
    // double factorAB = 1.0 / ( numSubjA - 1 + numSubjB - 1) * ( 1.0 / numSubjA + 1.0 / numSubjB);
    double factorA = 1.0 / ( numSubjA - 1 ) / numSubjA;
    double factorB = 1.0 / ( numSubjB - 1 ) / numSubjB;
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      for( int row = 0; row < tupelSize; row++ )
        {
        for( int column = 0; column < tupelSize; column++ )
          {
          (covarAB[feat])[row][column] =
            // ((covarA[feat])[row][column] + (covarB[feat])[row][column] ) * factorAB;
            // original T^2
            (covarA[feat])[row][column] * factorA + (covarB[feat])[row][column] * factorB;
          // modified T^2 statistics that is more robust if n1!=n2 and unequal covariance
          }
        }
      }
    // invert covarAB for each covariance matrix
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      InvcovarAB[feat] = covarAB[feat].GetInverse();
      }

    double diffmean[tupelSize];
    double firstTerm[tupelSize];
    // Hotelling metric computation
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      for( int dim = 0; dim < tupelSize; dim++ )
        {
        diffmean[dim] = meanA[feat * tupelSize + dim] - meanB[feat * tupelSize + dim];
        }
      // (meanA - meanB)' * CovarAB^(-1) * (meanA - meanB)
      // first Matrix multiplication
      for( int column = 0; column < tupelSize; column++ )
        {
        firstTerm[column] = 0;
        for( int row = 0; row < tupelSize; row++ )
          {
          firstTerm[column] =
            firstTerm[column] + diffmean[row] * (InvcovarAB[feat])[row][column];
          }
        }
      // second Dot-product
      double dotVal = 0;
      for( int dim = 0; dim < tupelSize; dim++ )
        {
        dotVal = dotVal + firstTerm[dim] * diffmean[dim];
        }

      statistic[feat] = dotVal;
      }
    // compute 95th quantile
    quant95Stat = computeQuantile(numFeatures / tupelSize, statistic, 0.95);
    // compute mean
    avgStat = 0;
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      avgStat += statistic[feat];
      }
    avgStat /= numFeatures / tupelSize;
    }
  else if( tupel_size == 2 )
    {
    const int tupelSize = 2; // needs this due to hardcoded size of covariance matrix
    typedef itk::Matrix<double, tupelSize, tupelSize> CovMatType;

    static CovMatType * covarA   = new CovMatType[numFeatures / tupelSize];
    static CovMatType * covarB   = new CovMatType[numFeatures / tupelSize];
    static CovMatType * covarAB  = new CovMatType[numFeatures / tupelSize];
    static CovMatType * InvcovarAB  = new CovMatType[numFeatures / tupelSize];
    if( ( curNumFeatures != numFeatures) )
      {
      delete covarA;
      delete covarB;
      delete covarAB;
      delete InvcovarAB;
      covarA   = new CovMatType[numFeatures / tupelSize];
      covarB   = new CovMatType[numFeatures / tupelSize];
      covarAB  = new CovMatType[numFeatures / tupelSize];
      InvcovarAB  = new CovMatType[numFeatures / tupelSize];
      }
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      covarA[feat].Fill(0);
      covarB[feat].Fill(0);
      }
    for( int subj = 0; subj < numSubjects; subj++ )
      {
      int subjIndex = subj * numFeatures;
      if( groupLabel[subj] == GROUP_A_LABEL )
        {
        for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
          {
          for( int row = 0; row < tupelSize; row++ )
            {
            for( int column = 0; column < tupelSize; column++ )
              {
              int indexDiff1 =  feat * tupelSize + row;
              int indexDiff2 =  feat * tupelSize + column;
              (covarA[feat])[row][column] =  (covarA[feat])[row][column]
                + (featureValue[subjIndex + indexDiff1] - meanA[indexDiff1])
                * (featureValue[subjIndex + indexDiff2] - meanA[indexDiff2]);

              }
            }
          }
        }
      else if( groupLabel[subj] == GROUP_B_LABEL )
        {
        for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
          {
          for( int row = 0; row < tupelSize; row++ )
            {
            for( int column = 0; column < tupelSize; column++ )
              {
              int indexDiff1 =  feat * tupelSize + row;
              int indexDiff2 =  feat * tupelSize + column;
              (covarB[feat])[row][column] =  (covarB[feat])[row][column]
                + (featureValue[subjIndex + indexDiff1] - meanB[indexDiff1])
                * (featureValue[subjIndex + indexDiff2] - meanB[indexDiff2]);
              }
            }
          }
        }
      else
        {
        std::cerr << " group label " << groupLabel[subj] << " does not exist" << std::endl;
        }
      }
    // covarA = Sa * (na-1)
    // covarB = Sb * (nb-1)

    // compute the pooled covariance matrices [nxn] covarAB,
    // matrix is fully symmetric -> does not matter whether column or row first
    //    double factorAB = 1.0 / ( numSubjA - 1 + numSubjB - 1) * ( 1.0 / numSubjA + 1.0 / numSubjB);
    double factorA = 1.0 / ( numSubjA - 1 ) / numSubjA;
    double factorB = 1.0 / ( numSubjB - 1 ) / numSubjB;
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      for( int row = 0; row < tupelSize; row++ )
        {
        for( int column = 0; column < tupelSize; column++ )
          {
          (covarAB[feat])[row][column] =
            // ((covarA[feat])[row][column] + (covarB[feat])[row][column] ) * factorAB;
            // original T^2
            (covarA[feat])[row][column] * factorA + (covarB[feat])[row][column] * factorB;
          // modified T^2 statistics that is more robust if n1!=n2 and unequal covariance
          // => Seber, G.A.F., 1984, Multivariate Observations, Wiley
          // under assumption of homogenous covariance matrices identical with original T^2
          }
        }
      }
    // invert covarAB for each covariance matrix
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      InvcovarAB[feat] = covarAB[feat].GetInverse();
      }

    double diffmean[tupelSize];
    double firstTerm[tupelSize];
    // Hotelling metric computation
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      for( int dim = 0; dim < tupelSize; dim++ )
        {
        diffmean[dim] = meanA[feat * tupelSize + dim] - meanB[feat * tupelSize + dim];
        }
      // (meanA - meanB)' * CovarAB^(-1) * (meanA - meanB)
      // first Matrix multiplication
      for( int column = 0; column < tupelSize; column++ )
        {
        firstTerm[column] = 0;
        for( int row = 0; row < tupelSize; row++ )
          {
          firstTerm[column] =
            firstTerm[column] + diffmean[row] * (InvcovarAB[feat])[row][column];
          }
        }
      // second Dot-product
      double dotVal = 0;
      for( int dim = 0; dim < tupelSize; dim++ )
        {
        dotVal = dotVal + firstTerm[dim] * diffmean[dim];
        }

      statistic[feat] = dotVal;
      }
    // compute 95th quantile
    quant95Stat = computeQuantile(numFeatures / tupelSize, statistic, 0.95);
    // compute mean
    avgStat = 0;
    for( int feat = 0; feat < numFeatures / tupelSize; feat++ )
      {
      avgStat += statistic[feat];
      }
    avgStat /= numFeatures / tupelSize;
    }

}
