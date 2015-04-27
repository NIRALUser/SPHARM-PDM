/*
 * compute the spherical harmonics coefficient spatial object from a parametrized surface mesh
 *
 * author:  Martin Styner
 *
 */

#ifdef _MSC_VER
#pragma warning ( disable: 4786 )
#pragma warning ( disable: 4284 )
#endif

#include <iostream>
#include <vector>
#include <string>
#include <string.h>
#include <itkMatrix.h>
#include <itkVector.h>
#include <vnl/algo/vnl_symmetric_eigensystem.h>
#include "ParametricMeshToSPHARMSpatialObjectFilter.h"

#include "SphericalHarmonicSpatialObject.h"
#include "SphericalHarmonicMeshSource.h"

extern "C" {
#include "f2c.h"
#include "clapack.h"
int sgels_(char *trans, integer *m, integer *n, integer *
           nrhs, real *a, integer *lda, real *b, integer *ldb, real *work, integer *lwork, integer *info);

}
// SOMEHOW not all lapack routine have made it into vxl's netlib directory (WHY?)
// and the ones for solving linear systems seem to be missing

namespace neurolib
{

ParametricMeshToSPHARMSpatialObjectFilter::ParametricMeshToSPHARMSpatialObjectFilter()
{
  this->SetNumberOfRequiredInputs(2);
  this->SetNumberOfRequiredInputs(1);

  SpatialObjectType::Pointer output
    = dynamic_cast<SpatialObjectType *>(this->MakeOutput(0).GetPointer() );
  this->SetNumberOfIndexedOutputs(1);
#ifdef _MSC_VER
  this->SetNthOutput(0, output.GetPointer() );
#else
  this->itk::ProcessObject::SetNthOutput(0, output.GetPointer() );
#endif

  m_FlipParametrizationIndex = 0;
  m_Degree = 12;
  m_leg = NULL;
  m_flatCoeffs = NULL;
  m_GenerateDataRegion = 0;
  m_GenerateDataNumberOfRegions = 0;
  m_FlipTemplate = NULL;
  m_ParaEllipseAlignment = true;
}

ParametricMeshToSPHARMSpatialObjectFilter::~ParametricMeshToSPHARMSpatialObjectFilter()
{
  if( m_leg )
    {
    delete m_leg;
    }
}

ParametricMeshToSPHARMSpatialObjectFilter::DataObjectPointer
ParametricMeshToSPHARMSpatialObjectFilter::MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType /* idx */)
{
  SpatialObjectType::Pointer SO = SpatialObjectType::New();

  return static_cast<itk::DataObject *>(SO.GetPointer() );
}

void
ParametricMeshToSPHARMSpatialObjectFilter::SetOutput(SpatialObjectType *output)
{
  itkWarningMacro(
    <<
    "SetOutput(): This method is slated to be removed from ITK.  Please use GraftOutput() in possible combination with DisconnectPipeline() instead." );
  this->ProcessObject::SetNthOutput(0, output);
}

void
ParametricMeshToSPHARMSpatialObjectFilter::GraftOutput(SpatialObjectType *graft)
{
  SpatialObjectType * output = this->GetOutput();

  if( output && graft )
    {
    // copy the meta-information
    output->CopyInformation( graft );
    }
}

void
ParametricMeshToSPHARMSpatialObjectFilter::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();
}

ParametricMeshToSPHARMSpatialObjectFilter::InputMeshType *
ParametricMeshToSPHARMSpatialObjectFilter::GetInputSurfaceMesh()
{
#ifdef _MSC_VER
  return dynamic_cast<InputMeshType *>(this->GetInput(0) );
#else
  return dynamic_cast<InputMeshType *>(this->itk::ProcessObject::GetInput(0) );
#endif
};

ParametricMeshToSPHARMSpatialObjectFilter::InputMeshType *
ParametricMeshToSPHARMSpatialObjectFilter::GetInputParametrizationMesh()
{
#ifdef _MSC_VER
  return dynamic_cast<InputMeshType *>(this->GetInput(1) );
#else
  return dynamic_cast<InputMeshType *>(this->itk::ProcessObject::GetInput(1) );
#endif
};

void
ParametricMeshToSPHARMSpatialObjectFilter::SetInputParametrizationMesh(InputMeshType * mesh)
{
#ifdef _MSC_VER
  this->SetNthInput(1, const_cast<InputMeshType *>( mesh ) );
#else
  this->itk::ProcessObject::SetNthInput(1, const_cast<InputMeshType *>( mesh ) );
#endif
};

void
ParametricMeshToSPHARMSpatialObjectFilter::SetInputSurfaceMesh(InputMeshType * mesh)
{
#ifdef _MSC_VER
  this->SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
#else
  this->itk::ProcessObject::SetNthInput(0, const_cast<InputMeshType *>( mesh ) );
#endif
};

ParametricMeshToSPHARMSpatialObjectFilter::SpatialObjectType *
ParametricMeshToSPHARMSpatialObjectFilter::GetOutput()
{
  if( GetNumberOfOutputs() < 1 )
    {
    return 0;
    }
#ifdef _MSC_VER
  return dynamic_cast<SpatialObjectType *>(this->GetOutput(0) );
#else
  return dynamic_cast<SpatialObjectType *>(this->itk::ProcessObject::GetOutput(0) );
#endif
};

ParametricMeshToSPHARMSpatialObjectFilter::SpatialObjectType *
ParametricMeshToSPHARMSpatialObjectFilter::GetOutput(unsigned int idx)
{
  if( GetNumberOfOutputs() < idx )
    {
    return 0;
    }

#ifdef _MSC_VER
  if ( _MSC_VER < 1500 )
  {
    return dynamic_cast<SpatialObjectType *>(this->GetOutput(idx) ) ;
  }
  else
  {
    return dynamic_cast<SpatialObjectType *>(this->itk::ProcessObject::GetOutput(idx) ) ;
  }
#else
  return dynamic_cast<SpatialObjectType *>(this->itk::ProcessObject::GetOutput(idx) );
#endif


};

// Subroutines to flip parameterization

int ParametricMeshToSPHARMSpatialObjectFilter::flip_ll(double n)
{
  return (int ) sqrt( (double) (n - 1) );
}

int ParametricMeshToSPHARMSpatialObjectFilter::flip_mm(double n)
{
  double ll = flip_ll(n);

  return (int) (n - ll * ll - 1);
}

int ParametricMeshToSPHARMSpatialObjectFilter::flip_m2(double n)
{
  double ll = flip_ll(n);

  return (int) ( (n - ll * ll) / 2);
}

double  ParametricMeshToSPHARMSpatialObjectFilter::flip_fu0(double n)
{
  double retval;

  if( flip_mm(n) == 0 )
    {
    retval = pow(-1.0, flip_ll(n) );
    }
  else if( flip_mm(n) % 2 == 1 )
    {
    retval = pow(-1.0, flip_ll(n) + flip_m2(n) );
    }
  else
    {
    retval = pow(-1.0, flip_ll(n) + flip_m2(n) + 1);
    }
  return retval;
}

double  ParametricMeshToSPHARMSpatialObjectFilter::flip_fu1(double n)
{
  double retval;

  if( flip_mm(n) == 0 )
    {
    retval = pow(-1.0, flip_ll(n) );
    }
  else if( flip_mm(n) % 2 == 1 )
    {
    retval = pow(-1.0, flip_ll(n) );
    }
  else
    {
    retval = pow(-1.0, flip_ll(n) + 1);
    }
  return retval;
}

double  ParametricMeshToSPHARMSpatialObjectFilter::flip_fu2(double n)
{
  return pow(-1.0, flip_m2(n) );
}

double  ParametricMeshToSPHARMSpatialObjectFilter::flip_reflect( double n)
{
  double retval;

  if( flip_mm(n) == 0 )
    {
    retval = 1.0;
    }
  else if( flip_mm(n) % 2 == 1 )
    {
    retval =  1.0;
    }
  else
    {
    retval = -1.0;
    }
  return retval;
}

void ParametricMeshToSPHARMSpatialObjectFilter::flipu0(const CoefListType * incoef, CoefListType * outcoef)
{
  outcoef->clear();
  for( unsigned int i = 0; i < incoef->size(); i++ )
    {
    CoefType coef = (*incoef)[i];
    double   fac = flip_fu0( i + 1);
    for( int dim = 0; dim < 3; dim++ )
      {
      coef[dim] = coef[dim] * fac;
      }
    outcoef->push_back(coef);
    }
}

void ParametricMeshToSPHARMSpatialObjectFilter::flipu1(const CoefListType * incoef, CoefListType * outcoef)
{
  outcoef->clear();
  for( unsigned int i = 0; i < incoef->size(); i++ )
    {
    CoefType coef = (*incoef)[i];
    double   fac = flip_fu1( i + 1);
    for( int dim = 0; dim < 3; dim++ )
      {
      coef[dim] = coef[dim] * fac;
      }
    outcoef->push_back(coef);
    }
}

void ParametricMeshToSPHARMSpatialObjectFilter::flipu2(const CoefListType * incoef, CoefListType * outcoef)
{
  outcoef->clear();
  for( unsigned int i = 0; i < incoef->size(); i++ )
    {
    CoefType coef = (*incoef)[i];
    double   fac = flip_fu2( i + 1);
    for( int dim = 0; dim < 3; dim++ )
      {
      coef[dim] = coef[dim] * fac;
      }
    outcoef->push_back(coef);
    }
}

void ParametricMeshToSPHARMSpatialObjectFilter::reflect(const CoefListType * incoef, CoefListType * outcoef)
{
  outcoef->clear();
  for( unsigned int i = 0; i < incoef->size(); i++ )
    {
    CoefType coef = (*incoef)[i];
    double   fac = flip_reflect( i + 1);
    for( int dim = 0; dim < 3; dim++ )
      {
      coef[dim] = coef[dim] * fac;
      }
    outcoef->push_back(coef);
    }
}

double ParametricMeshToSPHARMSpatialObjectFilter::FlipAxesDistance(const CoefListType *coef1,
                                                                   const CoefListType * coef2)
{
  double diff[3];
  double sum = 0;
  for( int coefIndex = 1; coefIndex < 4; coefIndex++ )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      diff[dim] = ( (*coef1)[coefIndex])[dim] - ( (*coef2)[coefIndex])[dim];
      sum += diff[dim] * diff[dim];
      }
    }
  return sum;
}

double ParametricMeshToSPHARMSpatialObjectFilter::CoefDistance(const CoefListType *coef1, const CoefListType * coef2)
{
  double diff[3];
  double sum = 0;

  unsigned int size = coef1->size();

  if( size > coef2->size() )
    {
    size = coef2->size();
    }
  for( unsigned int coefIndex = 1; coefIndex < size; coefIndex++ )
    {
    for( int dim = 0; dim < 3; dim++ )
      {
      diff[dim] = ( (*coef1)[coefIndex])[dim] - ( (*coef2)[coefIndex])[dim];
      sum += diff[dim] * diff[dim];
      }
    }
  return sum;
}

void ParametricMeshToSPHARMSpatialObjectFilter::EllipseScaleAlign(const CoefListType * incoef, CoefListType * outcoef)
{
  outcoef->clear();
  // normalize size relative to volume = 4/3 Pi a b c
  double a =
    sqrt( ( (*incoef)[1])[0] * ( (*incoef)[1])[0] + ( (*incoef)[1])[1] * ( (*incoef)[1])[1] + ( (*incoef)[1])[2]
          * ( (*incoef)[1])[2]);
  double b =
    sqrt( ( (*incoef)[2])[0] * ( (*incoef)[2])[0] + ( (*incoef)[2])[1] * ( (*incoef)[2])[1] + ( (*incoef)[2])[2]
          * ( (*incoef)[2])[2]);
  double c =
    sqrt( ( (*incoef)[3])[0] * ( (*incoef)[3])[0] + ( (*incoef)[3])[1] * ( (*incoef)[3])[1] + ( (*incoef)[3])[2]
          * ( (*incoef)[3])[2]);
  double vol_factor = pow(a * b * c * 4 / 3 * M_PI, 1.0 / 3.0);

  EllipseAlign(incoef, outcoef, vol_factor);
}

void ParametricMeshToSPHARMSpatialObjectFilter::EllipseAlign(const CoefListType * incoef, CoefListType * outcoef)
{
  EllipseAlign(incoef, outcoef, 1.0);
}

void ParametricMeshToSPHARMSpatialObjectFilter::EllipseAlign(const CoefListType * incoef, CoefListType * outcoef,
                                                             double scale)
{
  outcoef->clear();
  // scale
  for( unsigned int i = 0; i < incoef->size(); i++ )
    {
    CoefType coef = (*incoef)[i];
    for( int dim = 0; dim < 3; dim++ )
      {
      coef[dim] = coef[dim] / scale;
      }
    outcoef->push_back(coef);
    }
  // translation
  for( int dim = 0; dim < 3; dim++ )
    {
    ( (*outcoef)[0])[dim] = 0;
    }

  double a =
    sqrt( ( (*outcoef)[1])[0] * ( (*outcoef)[1])[0] + ( (*outcoef)[1])[1] * ( (*outcoef)[1])[1] + ( (*outcoef)[1])[2]
          * ( (*outcoef)[1])[2]);
  double b =
    sqrt( ( (*outcoef)[2])[0] * ( (*outcoef)[2])[0] + ( (*outcoef)[2])[1] * ( (*outcoef)[2])[1] + ( (*outcoef)[2])[2]
          * ( (*outcoef)[2])[2]);
  double c =
    sqrt( ( (*outcoef)[3])[0] * ( (*outcoef)[3])[0] + ( (*outcoef)[3])[1] * ( (*outcoef)[3])[1] + ( (*outcoef)[3])[2]
          * ( (*outcoef)[3])[2]);
  // rotation axes
  CoefType e1 =  (*outcoef)[2];
  for( int dim = 0; dim < 3; dim++ )
    {
    e1[dim] = e1[dim] / b;
    }
  CoefType e2 =  (*outcoef)[3];
  for( int dim = 0; dim < 3; dim++ )
    {
    e2[dim] = e2[dim] / c;
    }
  CoefType e3 =  (*outcoef)[1];
  for( int dim = 0; dim < 3; dim++ )
    {
    e3[dim] = e3[dim] / a;
    }
  // std::cout << "c" << (*incoef)[1][0] << "," <<  (*incoef)[1][1] << "," << (*incoef)[1][2] << std::endl;
  // std::cout << "c" << (*incoef)[2][0] << "," <<  (*incoef)[2][1] << "," << (*incoef)[2][2] << std::endl;
  // std::cout << "c" << (*incoef)[3][0] << "," <<  (*incoef)[3][1] << "," << (*incoef)[3][2] << std::endl;
  // std::cout <<  std::endl;
  // rotNormTRS[coef_] := coef(x,y,z).{{e1,e2,e3}; (1x3) (3x3) ->(1x3)
  // x * e11 + y * e12 + z * e13 ...
  for( unsigned int i = 0; i < outcoef->size(); i++ )
    {
    CoefType coef = (*outcoef)[i];
    CoefType Tcoef = (*outcoef)[i];
    Tcoef[0] = coef[0] * e1[0] + coef[1] * e1[1] + coef[2] * e1[2];
    Tcoef[1] = coef[0] * e2[0] + coef[1] * e2[1] + coef[2] * e2[2];
    Tcoef[2] = coef[0] * e3[0] + coef[1] * e3[1] + coef[2] * e3[2];
    (*outcoef)[i] = Tcoef;
    }

  // std::cout << "c" << (*outcoef)[1][0] << "," <<  (*outcoef)[1][1] << "," << (*outcoef)[1][2] << std::endl;
  // std::cout << "c" << (*outcoef)[2][0] << "," <<  (*outcoef)[2][1] << "," << (*outcoef)[2][2] << std::endl;
  // std::cout << "c" << (*outcoef)[3][0] << "," <<  (*outcoef)[3][1] << "," << (*outcoef)[3][2] << std::endl;
  // std::cout <<  std::endl;
}

/** adjust coefficients such that first order ellipsoid axis build a right hand coordinate system*/
void ParametricMeshToSPHARMSpatialObjectFilter::MakeEllipsoidAxisRightHand(CoefListType * coef)
{
  double a =
    sqrt(
      ( (*coef)[1])[0] * ( (*coef)[1])[0] + ( (*coef)[1])[1] * ( (*coef)[1])[1] + ( (*coef)[1])[2] * ( (*coef)[1])[2]);
  double b =
    sqrt(
      ( (*coef)[2])[0] * ( (*coef)[2])[0] + ( (*coef)[2])[1] * ( (*coef)[2])[1] + ( (*coef)[2])[2] * ( (*coef)[2])[2]);
  double c =
    sqrt(
      ( (*coef)[3])[0] * ( (*coef)[3])[0] + ( (*coef)[3])[1] * ( (*coef)[3])[1] + ( (*coef)[3])[2] * ( (*coef)[3])[2]);

  itk::Vector<double, 3> e1, e2, e3;
  for( int dim = 0; dim < 3; dim++ )
    {
    e1[dim] = ( (*coef)[2])[dim] / b;
    }
  for( int dim = 0; dim < 3; dim++ )
    {
    e2[dim] = ( (*coef)[3])[dim] / c;
    }
  for( int dim = 0; dim < 3; dim++ )
    {
    e3[dim] = ( (*coef)[1])[dim] / a;
    }

  itk::Vector<double, 3> crossVec = itk::CrossProduct(e1, e2);
  double                 dotProduct = e3 * crossVec;

  //   std::cout <<" ellipsoid axes: (" << e1[0] << "," << e1[1] << ","<< e1[2] << "), ";
  //   std::cout << "(" << e2[0] << "," << e2[1] << ","<< e2[2] << "), ";
  //   std::cout << "(" << e3[0] << "," << e3[1] << ","<< e3[2] << ")" << std:: endl;

  if( dotProduct < 0 )   // ellipsoid axes span left handed coordinate system
    {
    FlipAlign(coef, 7); // flip x - axis
    std::cout << "right handed first order ellipsoid axes enforced" << std::endl;
    }

}

/** Flip the parametrization for optimal match with FlipTemplate */
void ParametricMeshToSPHARMSpatialObjectFilter::BestFlip(CoefListType * coef, const CoefListType * flipTemplate)
{
  CoefListType alignFlipTemplate, alignFitFlipTemplate;
  CoefListType alignCoef, alignCoef1, alignCoef2, coefTemplate;

  EllipseScaleAlign(flipTemplate, &alignFlipTemplate);
  EllipseScaleAlign(coef, &alignCoef);
  alignFitFlipTemplate = alignFlipTemplate;

  const int numPoss = 8;
  double    possParaMat[3 * numPoss] = {1, 1, 1,     -1, 1, 1,     1, -1, 1,     1, 1, -1,
                                        -1, -1, 1,    -1, 1, -1,    1, -1, -1,   -1, -1, -1};
  // possParaMat test all mirrorings of the coordinates rather than the parameter space
  // parameter space mirrorings are testing in BestFlipAlign

  double minDist = 100000;
  int    minIndex = -1;
  int    minDistFlip = -1;
  for( int index = 0; index < numPoss; index++ )
    {
    // flip only the 1-4 coefficient -> the ones of the first order ellipsoid
    for( int matDim = 0; matDim < 3; matDim++ )
      {
      CoefType coefEntry = alignFlipTemplate[matDim + 1];
      for( int dim = 0; dim < 3; dim++ )
        {
        coefEntry[dim] = coefEntry[dim] * possParaMat[index * 3 + matDim];
        }
      alignFitFlipTemplate[matDim + 1] = coefEntry;
      }
    alignCoef1 = alignCoef;
    int flipID = BestFlipAlign(&alignCoef1, &alignFitFlipTemplate);
    // BestFlipAlign will only test the first order ellipsoid axis for best alignment

    EllipseScaleAlign(&alignCoef1, &alignCoef2);

    double dist = CoefDistance(&alignCoef2, &alignFlipTemplate); // distance to non-mirrored template
    // std::cout << "Dist " << flipID << "," << index <<":" << dist << std::endl;
    if( dist < minDist )
      {
      minDist = dist;
      minIndex = index;
      minDistFlip = flipID;
      }
    }

  std::cout << "BestFlip: " << minDistFlip << "," << minIndex << std::endl;

  FlipAlign(coef, minDistFlip);
}

void ParametricMeshToSPHARMSpatialObjectFilter::FlipAlign(CoefListType * coef, int  flipIndex)
{
  CoefListType curflipCoef;
  CoefListType curflipCoef2;

  curflipCoef = *coef;

  if( flipIndex == 1 )
    {
    flipu0(coef, &curflipCoef);
    }
  else if( flipIndex == 2 )
    {
    flipu1(coef, &curflipCoef);
    }
  else if( flipIndex == 3 )
    {
    flipu2(coef, &curflipCoef);
    }
  else if( flipIndex == 4 )
    {
    reflect(coef, &curflipCoef);
    }
  else if( flipIndex == 5 )
    {
    reflect(coef, &curflipCoef2);
    flipu0(&curflipCoef2, &curflipCoef);
    }
  else if( flipIndex == 6 )
    {
    reflect(coef, &curflipCoef2);
    flipu1(&curflipCoef2, &curflipCoef);
    }
  else if( flipIndex == 7 )
    {
    reflect(coef, &curflipCoef2);
    flipu2(&curflipCoef2, &curflipCoef);
    }
  (*coef) = curflipCoef;
}

/** Flip the parametrization for optimal match with FlipTemplate,
    returns the index of the best flip function */
int ParametricMeshToSPHARMSpatialObjectFilter::BestFlipAlign(CoefListType * coef, const CoefListType * flipTemplate)
{
  // check all functions for improved match: Identity, flipu0, flipu1, flipu2, reflect,
  //  flipu0[reflect[]], flipu1[reflect[]], flipu2[reflect[]]

  // Identity
  CoefListType curflipCoef;
  CoefListType curflipCoef2;
  CoefListType bestFlipCoef = *coef;
  double       bestFlipMatch = FlipAxesDistance(coef, flipTemplate);
  double       curMatch = bestFlipMatch;
  int          bestFlipIndex = 0; // identity

  // std::cout  << "BestFlip Vals: " << curMatch << ",";  std::cout.flush();

  flipu0(coef, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 1;
    }
  // std::cout << curMatch << ",";   std::cout.flush();

  flipu1(coef, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 2;
    }
  // std::cout << curMatch << ",";   std::cout.flush();

  flipu2(coef, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 3;
    }
  // std::cout << curMatch << ",";   std::cout.flush();

  reflect(coef, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 4;
    }
  // std::cout << curMatch << ",";   std::cout.flush();

  reflect(coef, &curflipCoef2);
  flipu0(&curflipCoef2, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 5;
    }
  // std::cout << curMatch << ",";   std::cout.flush();

  reflect(coef, &curflipCoef2);
  flipu1(&curflipCoef2, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 6;
    }
  // std::cout << curMatch << ",";   std::cout.flush();

  reflect(coef, &curflipCoef2);
  flipu2(&curflipCoef2, &curflipCoef);
  curMatch = FlipAxesDistance(&curflipCoef, flipTemplate);
  if( curMatch < bestFlipMatch )
    {
    bestFlipMatch = curMatch; bestFlipCoef = curflipCoef;
    bestFlipIndex = 7;
    }
  // std::cout << curMatch << " --> " << bestFlipMatch << std::endl;

  (*coef) = bestFlipCoef;

  return bestFlipIndex;
}

ParametricMeshToSPHARMSpatialObjectFilterLegendre::ParametricMeshToSPHARMSpatialObjectFilterLegendre(const int l_in) :
  SQRT2(sqrt( (double) 2) ),  L(l_in)
{
  int    i, j, l, m, sgn;
  double n1, n2;

  // compute normalizing factors
  a = new double[(L + 1) + L * (L + 1) / 2];
  for( i = 0, l = 0; l <= L; l++ )
    {
    for( m = 0; m <= l; m++ )
      {
      a[i++] = A(l, m);
      }
    }
  // create associated Legendre polynomials
  // allocate array to array of polynomials
  plm = (double * * *) malloc( (L + 1) * sizeof(void *) ); // new double**[L+1];
  assert(plm != 0);
  for( l = 0; l <= L; l++ )
    {
    // allocate array of polynomials
    plm[l] = (double * *) malloc( (l + 1) * sizeof(void *) ); // new double*[l+1];
    assert(plm[l] != 0);
    for( m = 0; m <= l; m++ )
      {
      // allocate array for coeff's of polynomial
      (plm[l])[m] = (double *) malloc( (l - m + 2 + 1) * sizeof(double) );
      assert( (plm[l])[m] != 0);
      ( (plm[l])[m])[l - m + 1] = 0.0;
      ( (plm[l])[m])[l - m + 2] = 0.0;
      }
    }
  // P(0,0)
  ( (plm[0])[0])[0] = 1.0; ( (plm[0])[0])[1] = 0.0;  ( (plm[0])[0])[2] = 0.0;
  // P(1,0)
  ( (plm[1])[0])[0] = 0.0; ( (plm[1])[0])[1] = -1.0; // *(-1)^m !!!!!!!
  ( (plm[1])[0])[2] = 0.0; ( (plm[1])[0])[3] = 0.0;
  // P(1,1)
  ( (plm[1])[1])[0] = 1.0; ( (plm[1])[1])[1] = 0.0; ( (plm[1])[1])[2] = 0.0;
  for( l = 2; l <= L; l++ )
    {
    // compute coefficients of Legendre polynomials (recursive)
    n1 = (double)(l + l - 1); n2 = (double)(l - 1);
    // cout << "l" << l << ":";
    for( j = l; j > 0; j-- )
      {
      ( (plm[l])[0])[j] = (n1 * ( (plm[l - 1])[0])[j - 1] - n2 * ( (plm[l - 2])[0])[j])
        / (double)l;
      // cout << " " << ((plm[l])[0])[j];
      }
    ( (plm[l])[0])[0] = -n2 * ( (plm[l - 2])[0])[0] / (double)l;
    // cout << " " << ((plm[l])[0])[0];
    sgn = 1;
    for( m = 1; m <= l; m++ )
      {
      // differentiate Legendre polynomials m times
      // cout << " m" << m << ":";
      sgn *= -1;
      for( j = l - m; j >= 0; j-- )
        {
        ( (plm[l])[m])[j] = sgn * (j + 1) * ( (plm[l])[m - 1])[j + 1];
        // cout << " " << ((plm[l])[m])[j];
        }
      }
    // cout << endl;
    }
}

ParametricMeshToSPHARMSpatialObjectFilterLegendre::~ParametricMeshToSPHARMSpatialObjectFilterLegendre()
{
  int l, m;
  for( l = 0; l <= L; l++ )
    {
    for( m = 0; m <= l; m++ )
      {
      free( (plm[l])[m]);
      }
    free(plm[l]);
    }
  free(plm);
  delete [] a;
  // cout << "destructor called!"<<endl;
}

double ParametricMeshToSPHARMSpatialObjectFilterLegendre::A(int l, int m)
{
  int    i;
  double j;

  // (m >=0) && (l >= m)
  j = 1.0;
  for( i = l + m; i > l - m; i-- )
    {
    j *= (double)i;
    }
  return sqrt( (double) ( (l + l + 1) / j) / (4 * M_PI) );
}

double ParametricMeshToSPHARMSpatialObjectFilterLegendre::P(int l, int m, float z)
{
  int    i;
  double f, p1o, p_1, p_2;

  if( m == 0 )
    {
    if( l == 0 )
      {
      return 1;
      }
    else
      {
      p_1 = z; p_2 = 1;
      for( i = 2; i <= l; i++ )
        {
        p1o = p_1;
        p_1 = ( (i + i - 1) * z * p_1 - (i - 1) * p_2) / double(i);
        p_2 = p1o;
        }
      return p_1;
      }
    }
  else
    {
    f = ( (plm[l])[m])[l - m];
    for( i = l - m - 1; i >= 0; i-- )
      {
      f = z * f + ( (plm[l])[m])[i];
      }
    return f;
    }
}

// Tabelliert P_m^m, P_{m+1}^m, ... P_l^m  in  p_ptr[0]... p_ptr[l-m].
// Das heisst p_ptr[ll-m] = P_ll^m,   fuer m <= ll <= l.
// Obiges bezieht sich auf den urspruenglich uebergebenen p_ptr.
void plgndr_row(const int l, const int m,
                const double x,
                const double somx2,             // somx2=sqrt((1.0-x)*(1.0+x));
                double *p_ptr)
{
  double pmm = 1.0;

  if( m > 0 )
    {
    double fact = 1.0;
    for( int i = 1; i <= m; i++ )
      {
      pmm *= -fact * somx2;
      fact += 2.0;
      }
    }
  *p_ptr++ = pmm;
  if( l > m )
    {
    double pmmp1 = x * (2 * m + 1) * pmm;
    *p_ptr++ = pmmp1;
    for( int ll = (m + 2); ll <= l; ll++ )
      {
      double pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
      *p_ptr++ = pll;
      pmm = pmmp1;
      pmmp1 = pll;
      }
    }
}

/*
inline void plgndr_row(const int l, const int m,
                       const double x,
                       double *p_ptr)
{
  plgndr_row(l,m,x sqrt((1.0-x)*(1.0+x)), p_ptr);
}
*/
/*********************************************************/

void ParametricMeshToSPHARMSpatialObjectFilterLegendre::ylm(int l, int m, float x, float y, float z,
                                                            double & re, double & im)
{
  int    i;
  double p,  p_ptr[256];

  switch( m )
    {
    case 0: re = 1.0; im = 0.0;
      break;
    case 1: re = x; im = y; // sincos(atan2(y, x), &im, &re);
      break;
    default:
      //
      re = (1.0 - z) * (1.0 + z);
      for( p = re, i = 1; i < m; i++ )
        {
        p *= re;
        }
      p = sqrt( (double) p);
      if( p < 1e-20 )
        {
        p = 0;
        }
      else
        {
        p = 1 / p;
        }
      im = sin(m * atan2(y, x) );
      re = cos(m * atan2(y, x) );
      re *= p; im *= p;
      /*
      for (re= x, im= y, i= 1; i<m; i++) {
        re1= re*x - im*y;
        im= re*y + im*x;
        re= re1;
      }
      */
    }
  plgndr_row(l, m, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), p_ptr);
  p = p_ptr[l - m];
  i = l * (l + 1) / 2 + m;
  re *= p * a[i] * SQRT2; im *= p * a[i] * SQRT2;
}

void ParametricMeshToSPHARMSpatialObjectFilterLegendre::Ylm(int l, int m, float x, float y, float z,
                                                            double & re, double & im)
{

  if( m < 0 )
    {
    ylm(l, -m, x, y, z, re, im);
    if( m % 2 == 0 )
      {
      im = -im;
      }
    else
      {
      re = -re;
      }
    }
  else
    {
    ylm(l, m, x, y, z, re, im);
    }
  // cout << x << "," << y << "," << z << ": ";
  // cout << re << " " << im << ": ";
}

void ParametricMeshToSPHARMSpatialObjectFilterLegendre::Ylm(int l, int m, float theta, float phi,
                                                            double & re, double & im)
{
  const int     MAX_L = 10000;
  static double p[MAX_L + 1];
  double        z;

  z = cos(theta);
  // sincos(m*phi, &im, &re);
  im = sin(m * phi);
  re = cos(m * phi);
  plgndr_row(l, m, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), p);
  re *= a[l * (l + 1) / 2 + m] * p[l - m] * SQRT2;
  im *= a[l * (l + 1) / 2 + m] * p[l - m] * SQRT2;
}

int
ParametricMeshToSPHARMSpatialObjectFilter::Get_BaseVal( PointsContainerPointer paraPoints, float * & A, int & nvert,
                                                        int & degree)
{
  nvert = paraPoints->Size();

  double * vert = new double[nvert * 3];
  for( int i = 0; i < nvert; i++ )
    {
    PointType curPoint =  paraPoints->GetElement(i);
    vert[3 * i + 0] = curPoint[0]; vert[3 * i + 1] = curPoint[1]; vert[3 * i + 2] = curPoint[2];
    }

  double re, im, z;

  if( m_leg )
    {
    delete m_leg;
    }
  m_leg = new ParametricMeshToSPHARMSpatialObjectFilterLegendre(m_Degree);

  double * plm = new double[(m_Degree + 1) * (m_Degree + 1) * 3];

  degree = m_leg->L + 1;
  degree *= degree; // fuer re und im
  A = new float[degree * nvert];
  for( int ind = 0, i = 0; i < nvert; i++, ind += 3 )
    {
    z = vert[ind + 2];

    plgndr_row(m_leg->L, 0, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), plm);
    for( int l = 0; l <= m_leg->L; l++ )
      {
      A[i + l * l * nvert] = m_leg->a[l * (l + 1) / 2] * plm[l];      // assign Y_l^0
      }
    plgndr_row(m_leg->L, 1, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), plm);

    im = sin(atan2(vert[ind + 1], vert[ind]) );
    re = cos(atan2(vert[ind + 1], vert[ind]) );
    for( int l = 1; l <= m_leg->L; l++ )
      {
      int j = l * l + 1;
      A[i + j * nvert] = m_leg->a[l * (l + 1) / 2 + 1] * plm[l - 1] * re;
      A[i + j * nvert + nvert] = m_leg->a[l * (l + 1) / 2 + 1] * plm[l - 1] * im;
      }
    for( int m = 2; m <= m_leg->L; m++ )
      {
      im = sin(m * atan2(vert[ind + 1], vert[ind]) );
      re = cos(m * atan2(vert[ind + 1], vert[ind]) );
      // double p= (1.0-z)*(1.0+z);
      // p= sqrt(pow(p, m));
      // if (p<1e-20) p= 0; else p= 1/p;
      // re*= p; im*= p;
      plgndr_row(m_leg->L, m, z, sqrt( (double) (1.0 - z) * (1.0 + z) ), plm);
      for( int l = m; l <= m_leg->L; l++ )
        {
        int j = l * l + m + m - 1;
        A[i + j * nvert] = m_leg->a[l * (l + 1) / 2 + m] * plm[l - m] * re;
        A[i + j * nvert + nvert] = m_leg->a[l * (l + 1) / 2 + m] * plm[l - m] * im;
        }
      }
    }

  delete plm;
  delete vert;

  return 0;
}

void
ParametricMeshToSPHARMSpatialObjectFilter::ComputeCoeffs()
{
  float * A;
  integer numRH, m, n;

  InputMeshType::Pointer surfMesh = GetInputSurfaceMesh();
  InputMeshType::Pointer paraMesh = GetInputParametrizationMesh();

  PointsContainerPointer surfPoints = surfMesh->GetPoints();
  PointsContainerPointer paraPoints = paraMesh->GetPoints();

  if( surfPoints->Size() != paraPoints->Size() )
    {
    throw ParametricMeshToSPHARMSpatialObjectFilterException(
            __FILE__, __LINE__,
            "ParametricMeshToSPHARMSpatialObjectFilter: surface and para do not have same number of points");
    }

  int m_int;
  int n_int;
  Get_BaseVal(paraPoints, A, m_int, n_int);
  // std::cout << "m " << m_int << ", n " << n_int << std::endl;
  delete A;
  A = NULL;
  if( m_flatCoeffs )
    {
    delete m_flatCoeffs;
    }
  m_flatCoeffs = new double[n_int * 3];
  for( int curDim = 0; curDim < 3; curDim++ )
    {
    Get_BaseVal(paraPoints, A, m_int, n_int);
    m = m_int;
    n = n_int;

    numRH = 1;
    float * obj = new float[m * numRH];
    for( int i = 0; i < m; i++ )
      {
      PointType curPoint = surfPoints->GetElement(i);
      obj[i] = curPoint[curDim];
      }

    char    trans[20] = "N";
    int     fac = n * m;
    integer workSize = fac * 2;
    integer info;
    float * work = new float[workSize];

    integer lda = m;
    integer ldb = m;

    sgels_(trans, &m, &n, &numRH, A, &lda, obj, &ldb, work, &workSize, &info);

    delete A;
    A = NULL;
    delete work;
    work = NULL;

    int curElem = 0;
    for( int j = 0, l = 0; l <= (int) m_Degree; l++ )
      {
      // x Re Yl0
      m_flatCoeffs[curDim + curElem * 3] = obj[j];
      curElem++;
      // x Im Yl0 = 0
      for( j++, m = 1; m <= l; m++, j++ )
        {
        // x Re Ylm
        m_flatCoeffs[curDim + curElem * 3] = obj[j];
        curElem++;
        // x Im Ylm
        j++;
        m_flatCoeffs[curDim + curElem * 3] = obj[j];
        curElem++;
        }
      }
    delete obj;
    }
  m_coeffs.clear();
  SpatialObjectType::ScalarType elem[3];
  for( int i = 0; i < n; i++ )
    {
    // xyz Re Yl0
    elem[0] = m_flatCoeffs[i * 3 + 0];
    elem[1] = m_flatCoeffs[i * 3 + 1];
    elem[2] = m_flatCoeffs[i * 3 + 2];

    CoefType coef = elem;
    m_coeffs.push_back(coef);
    }

}

void
ParametricMeshToSPHARMSpatialObjectFilter::RotateParametricMesh()
{

  typedef vnl_matrix<double>        InputMatrixType;
  typedef itk::Matrix<double, 3, 3> EigenVectorMatrixType;
  typedef itk::Vector<double, 3>    vector3DType;

  InputMatrixType MuvwMatrix(3, 3);
  MuvwMatrix.fill(0);
  MuvwMatrix(0, 1) = -0.345494;
  MuvwMatrix(1, 2) = -0.345494;
  MuvwMatrix(2, 0) = 0.488603;

  InputMatrixType coefMatrix(3, 3);
  for( int row = 0; row < 3; row++ )
    {
    for( int col = 0; col < 3; col++ )
      {
      coefMatrix(row, col) = m_flatCoeffs[(row + 1) * 3 + col];
      }
    }
  InputMatrixType origAxes(3, 3);
  origAxes = MuvwMatrix * coefMatrix;
  InputMatrixType origAxesT(3, 3);
  origAxesT = origAxes.transpose();

  InputMatrixType resMatrix(3, 3);
  resMatrix = origAxes * origAxesT;

  vnl_symmetric_eigensystem<double> EigenSystem(resMatrix);

  double ev[3];
  for( int dim = 0; dim < 3; dim++ )
    {
    ev[dim] = EigenSystem.get_eigenvalue(dim);
    }

  // std::cout << "values " << ev[0] <<"," << ev[1] <<"," << ev[2] <<"," << std::endl;

  int sortIndex[3] = {0, 1, 2};
  if( ev[1] > ev[0] && ev[1] > ev[2] )
    {
    sortIndex[0] = 1;
    if( ev[2] > ev[0] )
      {
      sortIndex[1] = 2;
      sortIndex[2] = 0;
      }
    else
      {
      sortIndex[1] = 0;
      sortIndex[2] = 2;
      }
    }
  else if( ev[2] > ev[0] && ev[2] > ev[1] )
    {
    sortIndex[0] = 2;
    if( ev[0] > ev[1] )
      {
      sortIndex[1] = 0;
      sortIndex[2] = 1;
      }
    else
      {
      sortIndex[1] = 1;
      sortIndex[2] = 0;
      }
    }
  else
    {
    sortIndex[0] = 0;
    if( ev[2] > ev[1] )
      {
      sortIndex[1] = 2;
      sortIndex[2] = 1;
      }
    else
      {
      sortIndex[1] = 1;
      sortIndex[2] = 2;
      }
    }

  //  std::cout << "values " << ev[sortIndex[0]] <<"," << ev[sortIndex[1]] <<"," << ev[sortIndex[2]] <<"," << std::endl;

  EigenVectorMatrixType RotT;
  for( int row = 0; row < 3; row++ )
    {
    vector3DType vec;
    vec.SetVnlVector(EigenSystem.get_eigenvector(sortIndex[row]) );
    for( int col = 0; col < 3; col++ )
      {
      RotT(row, col) = vec[col];
      }
    }

  InputMeshType::Pointer paraMesh = GetInputParametrizationMesh();
  PointsContainerPointer paraPoints = paraMesh->GetPoints();

  InputMeshType::Pointer RotParametricMesh = InputMeshType::New();
  PointsContainerPointer points = PointsContainer::New();

  PointType curPoint;
  for( unsigned int i = 0; i < paraPoints->Size(); i++ )
    {
    curPoint =  paraPoints->GetElement(i);
    vector3DType vert;
    vert[0] = curPoint[0]; vert[1] = curPoint[1]; vert[2] = curPoint[2];
    vert = RotT * vert;
    double vertex[3];
    vertex[0] = vert[0]; vertex[1] = vert[1]; vertex[2] = vert[2];
    points->InsertElement(i, PointType(vertex) );
    }
  RotParametricMesh->SetPoints(points); // connectivity is lost, but that is ok
  SetInputParametrizationMesh(RotParametricMesh);
}

ParametricMeshToSPHARMSpatialObjectFilter::CoefListType *
ParametricMeshToSPHARMSpatialObjectFilter::GetEllipseAlignCoef()
{
  CoefListType * coef = new CoefListType;

  EllipseAlign(&m_coeffs, coef);

  return coef;
}

ParametricMeshToSPHARMSpatialObjectFilter::CoefListType *
ParametricMeshToSPHARMSpatialObjectFilter::GetEllipseScaleAlignCoef()
{
  CoefListType * coef = new CoefListType;

  EllipseScaleAlign(&m_coeffs, coef);

  return coef;
}

/** Generate the data */
void
ParametricMeshToSPHARMSpatialObjectFilter::GenerateData()
{

  if( this->GetNumberOfInputs() < 2 )
    {
    throw ParametricMeshToSPHARMSpatialObjectFilterException(
            __FILE__, __LINE__,
            "ParametricMeshToSPHARMSpatialObjectFilter: 2 surfaces are needed, a surface and a parameterization mesh");
    }

  SpatialObjectType::Pointer SPHARMSO = SpatialObjectType::New();

  // std::cout << "raw coefs " << std::endl;
  ComputeCoeffs();
  SPHARMSO->SetCoefs(m_coeffs);

  if( m_ParaEllipseAlignment )
    {
    // Determine rotation matrix & rotate ParametricMesh
    // std::cout << "rotate para " << std::endl;
    RotateParametricMesh();
    m_coeffs.clear();

    // Recompute the coefficients
    // std::cout << "right hand" << std::endl;
    ComputeCoeffs();
    MakeEllipsoidAxisRightHand(&m_coeffs);
    SPHARMSO->SetCoefs(m_coeffs);
    }

  if( m_FlipTemplate.IsNotNull() )
    {
    CoefListType flipTemplateCoeffs;
    m_FlipTemplate->GetCoefs(flipTemplateCoeffs);
    BestFlip(&m_coeffs, &flipTemplateCoeffs);
    }

  // User-specified final flip
  if( (m_FlipParametrizationIndex > 0 && m_FlipParametrizationIndex <= 7) )
    {
    std::cout << "final flip " << m_FlipParametrizationIndex << std::endl;
    FlipAlign(&m_coeffs, m_FlipParametrizationIndex);
    }

  // std::cout << "set coeffs" << std::endl;
  SPHARMSO->SetCoefs(m_coeffs);

#ifdef _MSC_VER
  this->SetNthOutput(0, SPHARMSO.GetPointer() );
#else
  this->itk::ProcessObject::SetNthOutput(0, SPHARMSO.GetPointer() );
#endif

}

/** PrintSelf */
void
ParametricMeshToSPHARMSpatialObjectFilter::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  //   os << indent
  //      << "ObjectValue: "
  //      << static_cast<itk::NumericTraits<unsigned char>::PrintType>(m_ObjectValue)
  //      << std::endl;

}

}
