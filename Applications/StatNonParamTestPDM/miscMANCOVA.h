#ifndef MISCMANCOVA_H
#define MISCMANCOVA_H

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>
#include <vnl/vnl_vector_fixed.h>
#include <vnl/vnl_vector.h>

#include "constants.h"
#include "newTypes.h"

#include <string>
#include <vector.h>
#include <fstream>
#include <iostream>

#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataToitkMesh.h"

#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>

#include <itkImageRandomNonRepeatingConstIteratorWithIndex.h>

#include "boost/filesystem.hpp"

// contains a number of helper-functions for the MANCOVA shape
// analysis code

vnl_vector<double> TiedRankCalc(vnl_vector<double>& data);

vnl_vector<double> fdrCorrection( vnl_vector<double>& rawP, double fdrLevel , double & fdrThresh );

vnl_vector<double> bonferroniCorrection( vnl_vector<double>& rawP );

void output_vector(vnl_vector<double> data, std::string outbase, std::string toAppend);
void output_matrix(vnl_matrix<double> data, std::string outbase, std::string toAppend);

void load_MeshList_file( std::string filename,int surfaceColumn, unsigned int numIndependent, 
			 unsigned int numGroupTypes, vector<int> columnGroupTypes,vector<int> columnIndependent,bool scaleOn, int scaleColumn,
			 unsigned int &numSubjects, unsigned int &numA, unsigned int &numB,
		    vnl_matrix<int>* &groupLabel, vnl_matrix<double>* &scaleFactor, 
			 std::string* &meshFileNames,  vnl_matrix<double> * &indValue, bool computeScaleFactorFromVolumes );


void load_Meshes( bool scaleOn, unsigned int numSubjects, 
		  unsigned int numIndependent, 
		  vnl_matrix<double>* &scaleFactor, 
		  vnl_matrix<double> * &indValue,
		  std::string* &meshFileNames, 
		  unsigned int & numFeatures, vnl_matrix<double>  * &featureValue,
		  MeshType::Pointer & surfaceMesh, MeshSpatialObjectType::Pointer & SOMesh);

void load_KWMreadableInputFile( bool scaleOn, unsigned int numSubjects, 
		  unsigned int numIndependent, 
		  vnl_matrix<double>* &scaleFactor, 
		  vnl_matrix<double> * &indValue,
		  std::string* &meshFileNames, 
		  unsigned int & numFeatures, vnl_matrix<double>  * &featureValue,
		  MeshType::Pointer & surfaceMesh, MeshSpatialObjectType::Pointer & SOMesh);

void write_SubjectPoints( std::string fileName, unsigned int sub, 
			  vnl_matrix<double>  * &featureValue, 
			  unsigned int numFeatures,
			  MeshType::Pointer & surfaceMesh, 
			  MeshSpatialObjectType::Pointer & SOMesh );

void write_SurfaceProperties(  std::string outbase,  PointsContainerPointer & meanPoints,  
			       PointsContainerPointer & meanPointsA,  
			       PointsContainerPointer & meanPointsB,   
			       vnl_matrix<double>& diffVectors,  
			       vnl_vector<double>& normProjections, 
			       vnl_vector<double>& normDistProjections, 
			       vnl_matrix<double>& zScoresProjected,
			       vnl_matrix<double>& zScores,
			       bool writeOutZScores,
			       MeshType::Pointer & surfaceMesh, 
			       MeshSpatialObjectType::Pointer & SOMesh, 
			       std::string* &meshFileNames , int KWMreadableInputFile);

void compute_SurfaceProperties( PointsContainerPointer & meanPoints, PointsContainerPointer & meanPointsA, PointsContainerPointer & meanPointsB, vnl_matrix<double> & diffVectors, vnl_vector<double>& normProjections, vnl_vector<double>& normDistProjections, vnl_matrix<double>& zScores, vnl_matrix<double>& zScoresProjected, unsigned int numSubjects, unsigned int numA, unsigned int numB, unsigned int numFeatures, vnl_matrix<int> * &groupLabel, vnl_matrix<double>* &featureValue, vtkPolyDataNormals *& MeshNormals, MeshType::Pointer & surfaceMesh, std::string outbase, std::string* &meshFileNames, int KWMreadableInputFile);

double computePearsonCorrelation( vnl_vector<double>& x, vnl_vector<double>& y );
double computePearsonCorrelationWithP( vnl_vector<double>& x, vnl_vector<double>& y, double &dP, double &dPP, double &dPN );

#endif
