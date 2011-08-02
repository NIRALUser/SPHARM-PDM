/***********************************
*** Written by Marc Niethammer  ***
*** 05/2009                      ***
*** Modified by Beatriz Paniagua ***
*** 09/2009          ***
************************************/

#ifndef NONMANCOVA_H
#define NONMANCOVA_H

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_math.h>
#include <vnl/vnl_vector.h>

#include "miscMANCOVA.h"

#include <itkImageRandomNonRepeatingConstIteratorWithIndex.h>

#include "constants.h"
#include "newTypes.h"

#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>

void do_ScalarInteractionTest( unsigned int numSubjects, unsigned int numFeatures, unsigned int testColumn,
                               vnl_matrix<double>  * & featureValue,  PointsContainerPointer & meanPoints,
                               vtkPolyDataNormals *& MeshNormals,
                               vnl_vector<double>& spearmanRhoDist, vnl_vector<double>& spearmanRhoPro,
                               vnl_vector<double>& spearmanRhoProPval,
                               vnl_vector<double>& spearmanRhoDistPval, vnl_vector<double>& pearsonRhoDist,
                               vnl_vector<double>& pearsonRhoPro,
                               vnl_vector<double>& pearsonRhoProPval, vnl_vector<double>& pearsonRhoDistPval,
                               int correlationType, bool useParametricP,
                               unsigned int numPerms,
                               std::string outbase );

#endif
