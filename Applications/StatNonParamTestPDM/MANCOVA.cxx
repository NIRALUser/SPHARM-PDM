#include "MANCOVA.h"

void precomputeStatisticQuantities( vnl_matrix<double>& A, vnl_matrix<double>& X, vnl_matrix<double>& precompAA,
                                    vnl_matrix<double>& precompXPseudoInverse )
{
  double zeroOutTolerance = 1e-8;

  precompAA = A * vnl_svd<double>(X.transpose() * X, zeroOutTolerance).inverse() * A.transpose();
  precompAA = vnl_svd<double>(precompAA).inverse();

  precompXPseudoInverse = vnl_svd<double>(X.transpose() * X, zeroOutTolerance).inverse() * X.transpose();

}

double computeStatisticValue( vnl_matrix<double>& A,
                              vnl_matrix<double>& B,
                              vnl_matrix<double>& C,
                              vnl_matrix<double>& X,
                              vnl_matrix<double>& Y,
                              vnl_matrix<double>& E,
                              vnl_matrix<double>& H,
                              vnl_matrix<double>& precompAA,
                              vnl_matrix<double>& precompXPseudoInverse,
                              vnl_matrix<double>& NDtmp,
                              vnl_matrix<double>& PPtmp,
                              vnl_matrix<double>& DDtmp,
                              vnl_matrix<double>& NNtmp,
                              vnl_matrix<double>& AAtmp,
                              vnl_matrix<double>& Ctmp, int computeStatisticsType )
{

  // computing the residual error

  B = precompXPseudoInverse * Y;

  NDtmp = Y - X * B;
  E = NDtmp.transpose() * NDtmp;

  // computing the hypothesis matrix

  Ctmp = A * B - C;

  H = Ctmp.transpose() * precompAA * Ctmp;

  DDtmp = vnl_svd<double>(E).inverse() * H;

  vnl_diag_matrix<std::complex<double> >::iterator it;          // Iterator for stepping through eigenvalues
  vnl_real_eigensystem                             eigs(DDtmp); // Compute the eigenvalues

  double test_statistic; // statistic to calculate

  switch( computeStatisticsType )
    {
    case WILKS: // Wilks's likelihood ratio
      test_statistic = 1;
      for( it = eigs.D.begin(); it != eigs.D.end(); ++it )
        {
        test_statistic /= (1 + it->real() );
        }
      break;
    case HOTELLING: // Hotelling
      test_statistic = 0;
      for( it = eigs.D.begin(); it != eigs.D.end(); ++it )
        {
        test_statistic += it->real();
        }
      break;
    case PILLAI: // Pillai
      test_statistic = 0;
      for( it = eigs.D.begin(); it != eigs.D.end(); ++it )
        {
        test_statistic += (it->real() ) / (1 + it->real() );
        }
      break;
    case ROY: // Roy
      test_statistic = (eigs.D.begin() )->real();
      for( it = eigs.D.begin(); it != eigs.D.end(); ++it )
        {
        if( it->real() > test_statistic )
          {
          test_statistic = it->real();
          }
        }
      break;
    default:
      std::cerr << "Error:Unknown test statistic: " << computeStatisticsType << std::endl;
      exit( -1 );
    }

  return test_statistic;
}

void doMANCOVATesting( unsigned int numSubjects, unsigned int numFeatures, unsigned int numGroupTypes,
                       unsigned int numIndependent, unsigned int testColumn, unsigned int numPerms,
                       vnl_matrix<int> * groupLabel,
                       vnl_matrix<double>  * & featureValue, bool interactionTest, unsigned int numA, unsigned int numB,
                       double significanceLevel, int computeStatisticsType, vnl_vector<double>& rawP,
                       PointsContainerPointer & meanPoints )
{
  // this is a complete re-write of the MANCOVA code (MN)

  // first create all the matrices we need (for the correct setup to
  // compute the desired test statistic)
  // afterwards we do a permutation test to see how well our real
  // values compare to the permuted ones

  // first determine all the sizes we will deal with

  std::cout << "Performing MANCOVA permutation testing." << std::endl;
  std::cout << "Using " << numPerms << " permutations." << std::endl;

  unsigned int nrGroupVariables = 0;
  if( numGroupTypes == 1 )
    {
    if( numA == 0 || numB == 0 )
      {
      nrGroupVariables = 0;  // there is only one group, so no
      // specific group variables
      }
    else
      {
      nrGroupVariables = 2;  // there is one group
      }
    }
  else
    {
    nrGroupVariables = (1 << numGroupTypes);
    }

  unsigned int p = nrGroupVariables + numIndependent;
  unsigned int n = numSubjects;
  unsigned int d = dimension;

  // set up the matrices

  vnl_matrix<double> A;
  vnl_matrix<double> B;
  vnl_matrix<double> C;
  vnl_matrix<double> X;
  vnl_matrix<double> XP;
  vnl_matrix<double> Y;
  vnl_matrix<double> H;
  vnl_matrix<double> E;

  // to compute the mean shape

  vnl_matrix<double> YAverage;

  // the A matrix

  if( interactionTest )
    {
    // matrix will be 1xp, with a one at the position of the variable
    // to be tested

    A.set_size(1, p);
    A.fill( 0 );
    A(0, nrGroupVariables + testColumn ) = 1; // test column is 0 indexed
    }
  else
    {
    // group test
    // matrix will contrast the corresponding groups
    A.set_size(numGroupTypes, p);
    A.fill( 0 );
    if( nrGroupVariables > 0 )
      {
      if( numGroupTypes > 2 )
        {
        std::cerr << "ERROR: Not yet implemented for more than two group types. ABORT." << std::endl;
        exit( -1 );
        }

      switch( numGroupTypes )
        {
        case 1:
          A(0, 0) = 1;
          A(0, 1) = -1;
          break;
        case 2:
          A(0, 0) = 1;
          A(0, 2) = -1;
          A(1, 1) = 1;
          A(1, 3) = -1;
          break;
        default: std::cerr << "ERROR: unsupported number of groups = " << numGroupTypes << ". ABORT." << std::endl;
          exit( -1 );
        }
      }
    else
      {
      std::cerr << "ERROR: Trying to run a group test, but no groups are defined. ABORT." << std::endl;
      exit( -1 );
      }
    }

  // the B matrix
  // ... will be pxd

  B.set_size( p, d ); // will contain the fitting coefficients
  B.fill( 0 );

  // the C matrix
  // ... simply a zero matrix, but with one row for interaction test
  // and two rows for group test, as many columns as spatial
  // dimensions

  if( interactionTest )
    {
    C.set_size(1, dimension);
    }
  else
    {
    if( numGroupTypes == 1 )
      {
      C.set_size(1, dimension);
      }
    else
      { // must be two now
      C.set_size(2, dimension);
      }
    }
  C.fill( 0 );

  // the X matrix
  // ... the design matrix
  // is constant in principle, but will be permuted in the permutation test

  X.set_size( n, p );
  X.fill( 0 );

  // fill up the non-permuted part

  if( nrGroupVariables > 0 ) // if we really have multiple groups
    {
    for( unsigned int iI = 0; iI < numSubjects; iI++ )
      {
      unsigned int which_group = 0;
      for( unsigned int gr = 0; gr < numGroupTypes; ++gr )
        {
        if( (*groupLabel)[iI][gr] == GROUP_A_LABEL )
          {
          which_group += (1 << gr); // this is a way to count in binary for all of the group type combinations
          }
        }
      X( iI, which_group ) = 1;
      }
    }
  // fill up the interaction info
  for( unsigned int iI = 0; iI < numSubjects; iI++ )
    {
    for( unsigned int iJ = 0; iJ < numIndependent; iJ++ )
      {
      X( iI, iJ + nrGroupVariables ) = (*featureValue)[iI][numFeatures * dimension + iJ];
      }
    }

  // XP will be used to compute the permutations

  XP.set_size( n, p );
  XP = X;

  // the Y matrix
  // ... the observation matrix
  // we do this point, by point
  // so it will change depending on the feature we look at

  Y.set_size( n, d );
  Y.fill( 0 );

  // the E matrix
  // ... the residual matrix

  E.set_size( p, p );
  E.fill( 0 );

  H.set_size( dimension, dimension );
  H.fill( 0 );

  // temporary matrices (so we only need to allocate once and can then
  // reuse)

  vnl_matrix<double> NDtmp;
  vnl_matrix<double> PPtmp;
  vnl_matrix<double> AAtmp;
  vnl_matrix<double> Ctmp;
  vnl_matrix<double> NNtmp;
  vnl_matrix<double> DDtmp;

  NDtmp.set_size( n, d );
  PPtmp.set_size( p, p );
  DDtmp.set_size( d, d );
  NNtmp.set_size( n, n );
  AAtmp.set_size( A.rows(), A.rows() );
  Ctmp.set_size( C.rows(), C.columns() );

  // compute the baseline statistic for all the surface points

  vnl_vector<double> origStatistics(numFeatures);

  vnl_matrix<double> precompAA;
  precompAA.set_size( A.rows(), A.rows() );

  vnl_matrix<double> precompXPseudoInverse;
  precompXPseudoInverse.set_size( p, n );

  precomputeStatisticQuantities( A, X, precompAA, precompXPseudoInverse );
  for( unsigned int feat = 0; feat < numFeatures; feat++ )
    {
    for( unsigned int tup = 0; tup < tupelSize; tup++ )
      {
      Y.set_column(tup, featureValue->get_column(feat * tupelSize + tup) );
      }

    origStatistics(feat) =
      computeStatisticValue( A, B, C, X, Y, E, H, precompAA, precompXPseudoInverse, NDtmp, PPtmp, DDtmp, NNtmp, AAtmp,
                             Ctmp,
                             computeStatisticsType );

    // let's also compute the mean shape while we are at it

    YAverage = X * B; // the B-matrix was computed
    // the average over all the subjects is the desired point
    // coordinate

    PointType curPoint;
    for( unsigned int dim = 0; dim < dimension; dim++ )
      {
      curPoint[dim] = 0;
      for( unsigned int sub = 0; sub < numSubjects; sub++ )
        {
        curPoint[dim] += YAverage(sub, dim);
        }
      curPoint[dim] /= numSubjects;

      }

    // this defines the current mean point
    meanPoints->InsertElement( feat, curPoint ); // write it out

    } // end feature loop

  // now let's do the permutation testing
  // we will do FDR correction later

  rawP.set_size( numFeatures );
  rawP.fill( 0 );

  // so far only implement it for the interaction test

  itk::RandomPermutation permutation(numSubjects);
  itk::RandomPermutation permutation2(numSubjects);
  double                 currentStatistic;
  for( unsigned int nP = 0; nP < numPerms; nP++ )
    {
    // get current permutation

    if( numPerms > 100 )
      {
      if( nP % (numPerms / 100) == 0 )
        {
        std::cout << 100 * nP / numPerms << "% done\r";
        std::cout.flush();
        }
      }

    permutation.Shuffle();

    // compute the current X

    if( interactionTest )
      { // we need different permutations, depending on if we do an
        // interaction test or if a group test is performed
      for( unsigned int sub = 0; sub < numSubjects; sub++ )
        {
        XP( sub, nrGroupVariables + testColumn ) = X( permutation[sub], nrGroupVariables + testColumn );
        }

      }
    else
      {
      // this is the group test
      // of importance here is that we do not permute across gender
      // currently a maximum of two group labels are supported
      //
      // if there is only one group column we simply permute labels
      //
      // if there are two group columns we only permute among
      // identical labels in column one (e.g., column 1 represents
      // gender)

      switch( numGroupTypes )
        {
        case 1: // simply permute the two first columns of X
          for( unsigned int sub = 0; sub < numSubjects; sub++ )
            {
            XP( sub, 0 ) = X( permutation[sub], 0 );
            XP( sub, 1 ) = X( permutation[sub], 1 );
            }
          break;
        case 2: // permute columns (0,2) and (1,3) respectively
          // need an extra permutation
          permutation2.Shuffle();
          for( unsigned int sub = 0; sub < numSubjects; sub++ )
            {
            XP( sub, 0 ) = X( permutation[sub], 0 );
            XP( sub, 2 ) = X( permutation[sub], 2 );

            XP( sub, 1 ) = X( permutation2[sub], 1 );
            XP( sub, 3 ) = X( permutation2[sub], 3 );

            }
          break;
        default: std::cerr << "ERROR: unsupported number of groups = " << numGroupTypes << ". ABORT." << std::endl;
          exit( -1 );
        }
      }

    precomputeStatisticQuantities( A, XP, precompAA, precompXPseudoInverse );
    // now go through all of the points with this permutation and
    // compute the statistic
    for( unsigned int feat = 0; feat < numFeatures; feat++ )
      {
      for( unsigned int tup = 0; tup < tupelSize; tup++ )
        {
        Y.set_column(tup, featureValue->get_column(feat * tupelSize + tup) );
        }

      currentStatistic =
        computeStatisticValue( A, B, C, XP, Y, E, H, precompAA, precompXPseudoInverse, NDtmp, PPtmp, DDtmp, NNtmp,
                               AAtmp,
                               Ctmp,
                               computeStatisticsType );

      // depending on the statistic type we need to look for a maximum
      // or for a minimum

      if( computeStatisticsType == ROY || computeStatisticsType == HOTELLING || computeStatisticsType == PILLAI )
        {

        if( currentStatistic > origStatistics(feat) )
          {
          rawP(feat) += 1; // we found one that was larger (for right sided)
          }
        }

      if( computeStatisticsType == WILKS )
        {
        if( currentStatistic < origStatistics(feat) )
          {
          rawP(feat) += 1; // we found one that was smaller (for left sided)
          }
        }

      }

    }
  for( unsigned int feat = 0; feat < numFeatures; feat++ )
    {
    rawP(feat) /= numPerms;
    }

  std::cout << "Completed MANCOVA permutation testing." << std::endl;

}
