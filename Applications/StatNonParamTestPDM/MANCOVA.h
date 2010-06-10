#ifndef MANCOVA_H
#define MANCOVA_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>
#include <vnl/algo/vnl_real_eigensystem.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include <itkImageRandomNonRepeatingConstIteratorWithIndex.h>

#include "miscMANCOVA.h"

#include "constants.h"
#include "newTypes.h"


void precomputeStatisticQuantities( vnl_matrix<double>& A, vnl_matrix<double>& X, vnl_matrix<double>& precompAA, vnl_matrix<double>& precompXPseudoInverse );

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
			      vnl_matrix<double>& Ctmp, int computeStatisticsType );


void doMANCOVATesting( unsigned int numSubjects, unsigned int numFeatures, unsigned int numGroupTypes, unsigned int numIndependent, unsigned int testColumn, unsigned int numPerms, vnl_matrix<int> * groupLabel, vnl_matrix<double>  * &featureValue, bool interactionTest, unsigned int numA, unsigned int numB, double significanceLevel, int computeStatisticsType, vnl_vector<double>& rawP, PointsContainerPointer &meanPoints  );


#endif
