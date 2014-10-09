#include <iostream>
#include <fstream>
#include <vector>
#include <itkPoint.h>
#include "argio.hh"
#include "SphericalHarmonicPolynomial.h"

#include "itkMeshSpatialObject.h"
#include "itkSpatialObjectWriter.h"
#include "itkSpatialObjectReader.h"
#include "itkMesh.h"

#include "SphericalHarmonicSpatialObject.h"
#include "SphericalHarmonicCoefficientFileReader.h"
#include "SphericalHarmonicMeshSource.h"

using namespace std ;
using namespace neurolib ;

void openFiles(ofstream *output, int nPts) ;
void Compute (SphericalHarmonicSpatialObject::CoefListType *coefs, double *icos, int nPts) ;
void ComputeRidge (SphericalHarmonicSpatialObject::CoefListType *coefs, double *icos, int nPts) ;
void NormalizeVector (double *x) ;

static int degree, level ;
static bool firstDer, secondDer, unit, normal, first, second, kappa, c, s, princDirs, princCurve, h, k, gradKappa, caustic ;
static bool parFile ;
static bool phiThetaFile ;
static char* baseName ;

// takes in a spharm coefficient file, and generates the evaluations of the first and second derivatives
int main(int argc, const char **argv)
{
  // get the arguments
  char *coefFileName;
  
  coefFileName = ipGetStringArgument(argv, "-coef", NULL);  
  baseName = ipGetStringArgument(argv, "-output", NULL);  
  firstDer = ipExistsArgument (argv, "-d");
  secondDer = ipExistsArgument (argv, "-dd");
  degree = ipGetIntArgument (argv, "-deg", 12) ;
  level = ipGetIntArgument (argv, "-level", -1) ;
  unit = ipExistsArgument ( argv, "-unit") ;
  normal = ipExistsArgument (argv, "-norm") ;
  first = ipExistsArgument (argv, "-I") ;
  second = ipExistsArgument (argv, "-II") ;
  kappa = ipExistsArgument (argv, "-kappa") ;
  gradKappa = ipExistsArgument (argv, "-gradKappa") ;
  caustic = ipExistsArgument (argv, "-caustic") ;
  c = ipExistsArgument (argv, "-c") ;
  s = ipExistsArgument (argv, "-s") ;
  h = ipExistsArgument (argv, "-h") ;
  k = ipExistsArgument (argv, "-k") ;

  princDirs = ipExistsArgument (argv, "-princ") ;
  princCurve = ipExistsArgument (argv, "-prinCurve") ;  
  parFile = ipExistsArgument (argv, "-parFile") ;
  phiThetaFile = ipExistsArgument (argv, "-phiTheta") ;
  // make sure the arguments are valid
  if ( !coefFileName || !baseName || ipExistsArgument(argv, "-help") )
  {
    cout << "Usage:" << endl ;
    cout << "SPHARMTool -coef CoefficientFileName -output OutputFileBaseName [-level subdivisionLevel] [-deg spharmDegree] [-d] [-dd] [-normal] [-unit] [-I] [-II] [-kappa] [-c] [-s] [-h] [-k] [-princ] [-prinCurve] [-parFile] [-phiThetaFile]" << endl ;
    return 0 ;
  }
  
  // read the coefficients in
  SphericalHarmonicSpatialObject::CoefListType coeflist;
  SphericalHarmonicCoefficientFileReader::Pointer reader = SphericalHarmonicCoefficientFileReader::New();
  reader->SetFileName(coefFileName);
  reader->Update();
  reader->GetOutput(coeflist);
  
  // Create an itkMesh
  SphericalHarmonicMeshSource::Pointer meshsrc = SphericalHarmonicMeshSource::New();
  meshsrc->SetCoefs(coeflist);
  meshsrc->SetDegree(degree);
  if ( level > 0 )
  {
    meshsrc->SetLevel (level) ;
  }
  meshsrc->Update();
  typedef SphericalHarmonicMeshSource::OutputMeshType MeshType;
  MeshType* mesh;
  mesh = meshsrc->GetOutput();
  double *icos = meshsrc->GetPrecomputedIcosList();

  // Create the mesh Spatial Object
  typedef itk::MeshSpatialObject<MeshType>        MeshSOType;
  typedef MeshType::MeshTraits                    MeshTrait;
  MeshSOType::Pointer meshSO = MeshSOType::New();
  meshSO->SetMesh(mesh);
  
  // Writing the .meta file
  typedef itk::SpatialObjectWriter<3,float,MeshTrait> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(meshSO);
  char meshName[350] ;
  sprintf (meshName, "%s.meta", baseName) ;
  writer->SetFileName(meshName);
  writer->Update();  
  
  ComputeRidge (&coeflist, icos, mesh->GetNumberOfPoints()) ;
  Compute (&coeflist, icos, mesh->GetNumberOfPoints()) ;
  return 0 ;
}

void ComputeRidge (SphericalHarmonicSpatialObject::CoefListType *coefs, double *icos, int nPts) 
{
  neurolib::SphericalHarmonicPolynomial<3> spharm ;
  spharm.SetDegree (degree) ;
  spharm.SetCoefs (*coefs) ;

  std::ofstream out, causticFile ;
  out.open ( "ridge.txt" ) ;
  out << "NUMBER_OF_POINTS=" << 5000 << endl ;
  out << "DIMENSION=3" << endl << "TYPE=Curve" << endl ;

  causticFile.open ("caustic.txt") ;
  causticFile << "NUMBER_OF_POINTS=" << 5000 << endl ;
  causticFile << "DIMENSION=3" << endl << "TYPE=Curve" << endl ;

  double point[3], uv[2] ;
  int family = 1 ;
  spharm.GetOneRidgePoint (0, degree, family, uv, point) ;
  out << point[0] << " " << point[1] << " " << point[2] << endl ;
  std::cout << point[0] << " " << point[1] << " " << point[2] << endl ;

  for ( int i = 0 ; i < 5000 ; i++ )
  {
    spharm.FollowRidge (0, degree, family, uv, uv, point, 0.005) ;
    out << point[0] << " " << point[1] << " " << point[2] << endl ;
    spharm.GetCaustic (0, degree, family, uv[0], uv[1], point) ;
    causticFile << point[0] << " " << point[1] << " " << point[2] << endl ;
    
  }
  out.close () ;
  causticFile.close () ;
}

void Compute (SphericalHarmonicSpatialObject::CoefListType *coefs, double *icos, int nPts) 
{
  const int nFiles = 23 ;
  ofstream output[nFiles] ;

  openFiles (output, nPts) ;

  neurolib::SphericalHarmonicPolynomial<3> spharm ;
  spharm.SetDegree (degree) ;
  spharm.SetCoefs (*coefs) ;

  double x[6] ;
  int i;

  try
  {
    double theta, phi ;
    for ( i=0; i < nPts ; i++ )
    {
      
      phi = icos[2*i] ;
      theta = icos[2*i+1] ;
    
      if ( phiThetaFile )
      {
        output[16] << theta << endl ;
      }
      if ( firstDer ) 
      {
        spharm.EvaluateDTheta (0, degree, phi, theta, x) ;
        if ( unit ) 
        {
          NormalizeVector (x) ;
        }
        output[0] << x[0] << " " << x[1] << " " << x[2] << endl ;
        
        spharm.EvaluateDPhi (0, degree, phi, theta, x) ;
        if ( unit ) 
        {
          NormalizeVector (x) ;
        }
        output[1] << x[0] << " " << x[1] << " " << x[2] << endl ;
      }
      if ( secondDer ) 
      {
        spharm.EvaluateDDThetaTheta (0, degree, phi, theta, x) ;
        if ( unit )
        {
          NormalizeVector (x) ;
        }
        output[2] << x[0] << " " << x[1] << " " << x[2] << endl ;
        
        spharm.EvaluateDDThetaPhi (0, degree, phi, theta, x) ;
        if ( unit )
        {
          NormalizeVector (x) ;
        }
        output[3] << x[0] << " " << x[1] << " " << x[2] << endl ;
          
        spharm.EvaluateDDPhiPhi (0, degree, phi, theta, x) ;
        if ( unit )
        {
          NormalizeVector (x) ;
        }
        output[4] << x[0] << " " << x[1] << " " << x[2] << endl ;
        
      }
      if ( normal )
      {
        spharm.GetUnitNormal (0, degree, phi, theta, x) ;
        output[5] << x[0] << " " << x[1] << " " << x[2] << endl ;
      }
      if ( first )
      {
        spharm.EvaluateFirstFundamentalForm (0, degree, phi, theta, x) ;
        output[6] << x[0] << " " << x[1] << " " << x[2] << endl ;
      }
      if ( second ) 
      {
        spharm.EvaluateSecondFundamentalForm (0, degree, phi, theta, x) ;
        output[7] << x[0] << " " << x[1] << " " << x[2] << endl ;
      }
      if ( kappa ) 
      {
        spharm.GetPrincipalCurvatures (0, degree, phi, theta, x) ;
        output[8] << x[0] << endl ;
        output[9] << x[1] << endl ;
      }
      if ( c ) 
      {
        x[0] = spharm.ComputeC (0, degree, phi, theta) ;
        output[10] << x[0] << endl ;
      }
      if ( s ) 
      {
        x[0] = spharm.ComputeS (0, degree, phi, theta) ;
        output[11] << x[0] << endl ;
      }
      if ( h ) 
      {
        x[0] = spharm.ComputeMeanCurvature (0, degree, phi, theta) ;
        output[17] << x[0] << endl ;
      }
      if ( k ) 
      {
        x[0] = spharm.ComputeGaussianCurvature (0, degree, phi, theta) ;
        output[18] << x[0] << endl ;
      }
      if ( princDirs )
      {
        spharm.GetPrincipalDirections (0, degree, phi, theta, x) ;
        output[12] << x[0] << " " << x[1] << " " << x[2] << endl ;
        output[13] << x[3] << " " << x[4] << " " << x[5] << endl ;
      }
      if ( gradKappa )
      {
        spharm.GetGradKappa (0, degree, 0, phi, theta, x) ;
        output[19] << x[0] << " " << x[1] << " " << x[2] << endl ;
        spharm.GetGradKappa (0, degree, 1, phi, theta, x) ;
        output[20] << x[0] << " " << x[1] << " " << x[2] << endl ;
      }
      
      if ( parFile )
      {
        x[0] = cos ( phi ) * sin ( theta ) ;
        x[1] = sin ( phi ) * sin ( theta ) ;
        x[2] = cos ( theta ) ;
        output[15] << x[0] << " " << x[1] << " " << x[2] << endl ;
        
      }

      if ( caustic )
      {
        spharm.GetCaustic (0, degree, 0, phi, theta, x) ;
        output[21] << x[0] << " " << x[1] << " " << x[2] << endl ;
        spharm.GetCaustic (0, degree, 1, phi, theta, x) ;
        output[22] << x[0] << " " << x[1] << " " << x[2] << endl ;
      }

    }
  }
  catch(SphericalHarmonicPolynomialException ex)
  {
    throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, ex.GetDescription());

  }

  if ( princCurve )
  {
    double curve[30000] ;
    spharm.ComputePrincipalCurve (0, degree, 0.01, 0.1, curve, 0.005, 10000) ;
    for ( i= 0 ; i < 10000 ; i++ )
    {
      output[14] << curve[3*i] << " " << curve[3*i+1] << " " << curve[3*i+2] << endl ;
    }
  }
  

  for ( i = 0 ; i < nFiles ; i++ )
  {
    output[i].close () ;
  }
}

void NormalizeVector (double *x)
{
  double length ;
  length = sqrt ( x[0] * x[0] + x[1] * x[1] + x[2] * x[2] ) ;
  if ( length != 0 )
  {
    x[0] /= length ;
    x[1] /= length ;
    x[2] /= length ;
  }
  
}

void openFiles(ofstream *output, int nPts) 
{
  char derName11[350], derName12[350], derName21[350], derName22[350], derName23[350], normName[350];
  char firstName[350], secondName[350], kappa1Name[350], kappa2Name[350], cName[350], sName[350], hName[350], kName[350] ;
  char princ1Name[350], princ2Name[350], curveName[350] ;
  char parFileName[350] ;
  char phiThetaFileName[350] ;
  char gradKappaFileName1[350], gradKappaFileName2[350] ;
  char causticFileName1[350], causticFileName2[350] ;

  sprintf (derName11, "%sDTheta.txt", baseName) ;
  sprintf (derName12, "%sDPhi.txt", baseName) ;
  sprintf (derName21, "%sDDThetaTheta.txt", baseName) ;
  sprintf (derName22, "%sDDThetaPhi.txt", baseName) ;
  sprintf (derName23, "%sDDPhiPhi.txt", baseName) ;
  sprintf (normName,  "%sNormal.txt", baseName) ;
  sprintf (firstName, "%sFirst.txt", baseName) ;
  sprintf (secondName, "%sSecond.txt", baseName) ;
  sprintf (kappa1Name, "%sKappa1.txt", baseName) ;
  sprintf (kappa2Name, "%sKappa2.txt", baseName) ;
  sprintf (cName, "%sC.txt", baseName) ;
  sprintf (sName, "%sS.txt", baseName) ;
  sprintf (hName, "%sH.txt", baseName) ;
  sprintf (kName, "%sK.txt", baseName) ;
  sprintf (princ1Name, "%sPrinc1.txt", baseName) ;
  sprintf (princ2Name, "%sPrinc2.txt", baseName) ;
  sprintf (curveName, "%sCurve.txt", baseName) ;
  sprintf (parFileName, "%sParametrization.par", baseName) ;
  sprintf (phiThetaFileName, "%sPhiTheta.txt", baseName) ;
  sprintf (gradKappaFileName1 ,"%sGradKappa1.txt", baseName ) ;
  sprintf (gradKappaFileName2 ,"%sGradKappa2.txt", baseName ) ;
  sprintf (causticFileName1 ,"%sCaustic1.txt", baseName ) ;
  sprintf (causticFileName2 ,"%sCaustic2.txt", baseName ) ;


  if ( firstDer ) 
  {
    output[0].open ( derName11 ) ;
    output[1].open ( derName12 ) ;

    output[0] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[0] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
    
    output[1] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[1] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
  }
  if ( secondDer ) 
  {
    output[2].open ( derName21 ) ;
    output[3].open ( derName22 ) ;
    output[4].open ( derName23 ) ;
  
    output[2] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[2] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
    
    output[3] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[3] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
    
    output[4] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[4] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
  }
  if ( normal )
  {
    output[5].open ( normName ) ;
  
    output[5] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[5] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
    
  }
  if ( first )
  {
    output[6].open ( firstName ) ;
  
    output[6] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[6] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
  }
  if ( second )
  {
    output[7].open ( secondName ) ;
  
    output[7] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[7] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
  }
  if ( kappa )
  {
    output[8].open ( kappa1Name ) ;
    output[9].open ( kappa2Name ) ;
    
    output[8] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[8] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
    
    output[9] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[9] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
  }
  if ( c )
  {
    output[10].open ( cName ) ;
  
    output[10] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[10] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
  }
  if ( s )
  {
    output[11].open ( sName ) ;
  
    output[11] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[11] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
  }
  if ( h )
  {
    output[17].open ( hName ) ;
  
    output[17] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[17] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
  }
  if ( k )
  {
    output[18].open ( kName ) ;
  
    output[18] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[18] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
  }
  if ( princDirs )
  {
    output[12].open ( princ1Name ) ;
    output[13].open ( princ2Name ) ;
    
    output[12] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[12] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
    
    output[13] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[13] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
  }
  if ( princCurve )
  {
    output[14].open ( curveName ) ;
  
    output[14] << "NUMBER_OF_POINTS=" << 10000 << endl ;
    output[14] << "DIMENSION=3" << endl << "TYPE=Curve" << endl ;
  }
  if ( parFile )
  {
    output[15].open ( parFileName ) ;
  }
  if ( phiThetaFile )
  {
    output[16].open ( phiThetaFileName ) ;
    output[16] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[16] << "DIMENSION=1" << endl << "TYPE=Scalar" << endl ;
  }
  if ( gradKappa )
  {
    output[19].open ( gradKappaFileName1 ) ;
    output[20].open ( gradKappaFileName2 ) ;
    
    output[19] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[19] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
    
    output[20] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[20] << "DIMENSION=3" << endl << "TYPE=Vector" << endl ;
  }
  if ( caustic )
  {
    output[21].open ( causticFileName1 ) ;
    output[22].open ( causticFileName2 ) ;
    
    output[21] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[21] << "DIMENSION=3" << endl << "TYPE=Curve" << endl ;
    
    output[22] << "NUMBER_OF_POINTS=" << nPts << endl ;
    output[22] << "DIMENSION=3" << endl << "TYPE=Curve" << endl ;
  }
  
}
