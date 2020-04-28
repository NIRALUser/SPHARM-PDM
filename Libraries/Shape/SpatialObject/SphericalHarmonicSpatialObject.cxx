/*=========================================================================

  Author: Christine Xu

=========================================================================*/

#include "SphericalHarmonicSpatialObject.h"
#include "SphericalHarmonicMeshSource.h"

#include <math.h>
#include <iostream>

namespace neurolib
{

/** Constructor */
SphericalHarmonicSpatialObject::SphericalHarmonicSpatialObject()
{
  this->SetTypeName("SphericalHarmonicSpatialObject");

  m_Count = 0;
  m_Harmonic = 0;

  m_Subdiv = 20;
  // m_Theta = 0;
  // m_Phi = 0;
  // m_Vertex = 0;
  // m_Triangs = 0;

  m_CoefsMeshSpatialObject = MeshSpatialObjectType::New();

}

/** Destructor */
SphericalHarmonicSpatialObject::~SphericalHarmonicSpatialObject()
{
  if( !m_Coefs.empty() )
    {
    m_Coefs.clear();
    }

  m_Count = 0;
  m_Harmonic = 0;

  // Run time error will be generated if the following is included!!
  // if(m_CoefsMeshSpatialObject.IsNotNull())
  //  m_CoefsMeshSpatialObject->Delete();

}

void SphericalHarmonicSpatialObject::PrintSelf(std::ostream& os, itk::Indent indent) const
{

  Superclass::PrintSelf( os, indent );
  os << "Count: " << m_Count << std::endl;
  os << "Harmonic: " << m_Harmonic << std::endl;

  os << "Coefficient List: " << std::endl;
  CoefListType::const_iterator iter = m_Coefs.begin();
  while( iter != m_Coefs.end() )
    {
    CoefType elem = *iter;
    os << elem[0] << "  " << elem[1] << "  " << elem[2] << std::endl;
    iter++;
    }

}

void SphericalHarmonicSpatialObject::SetCoefs(SphericalHarmonicSpatialObject::CoefListType& coeflist)
{
  m_Coefs = coeflist;
  m_Count = m_Coefs.size(); // Each item in the vector is a 3D point.
  m_Harmonic = (int) floor(sqrt(double(m_Count)) ) - 1;
  // std::cout<<"m_Count = " <<m_Count <<std::endl;
  ComputeHiddenMeshSpatialObject();
}

void SphericalHarmonicSpatialObject::GetCoefs(SphericalHarmonicSpatialObject::CoefListType& coeflist) const
{
  coeflist = m_Coefs;
}

/** p_coefsMeshSpatialObject is the hidden mesh to provide results
for the three geometric functions (IsInside,  ComputeBoundingBox).
Pay attention to the variable 'depth'?? */
bool SphericalHarmonicSpatialObject::ValueAtInObjectSpace( const PointType & /* point */,
                                              double & /* value */,
                                              unsigned int /* depth */,
                                              const std::string & /* name */) const
{
  return false;
}

bool SphericalHarmonicSpatialObject::IsEvaluableAtInObjectSpace(const PointType & /* point */, unsigned int /* depth */,
                             const std::string & /* name */) const
{
  return false;
}

bool SphericalHarmonicSpatialObject::IsInsideInObjectSpace(const PointType & point) const
{
  if( this->GetMyBoundingBoxInObjectSpace()->IsInside(point) )
    {
    return m_CoefsMeshSpatialObject->IsInsideInObjectSpace(point);
    }
  return false;
}

void SphericalHarmonicSpatialObject::ComputeMyBoundingBox() 
{
  itkDebugMacro( "Computing SphericalHarmonicSpatialObject object space bounding box" );
  
  
  m_CoefsMeshSpatialObject->Update();

  PointType ptMin, ptMax;
  ptMin = m_CoefsMeshSpatialObject->GetMyBoundingBoxInObjectSpace()->GetMinimum();
  ptMax = m_CoefsMeshSpatialObject->GetMyBoundingBoxInObjectSpace()->GetMaximum();


  this->GetModifiableMyBoundingBoxInObjectSpace()->SetMinimum(ptMin);
  this->GetModifiableMyBoundingBoxInObjectSpace()->SetMaximum(ptMax);
  
}

void SphericalHarmonicSpatialObject::ComputeHiddenMeshSpatialObject()
{
  SphericalHarmonicMeshSource::Pointer m_Meshsrc = SphericalHarmonicMeshSource::New();

  m_Meshsrc->SetCoefs(m_Coefs);
  m_Meshsrc->SetLevel(m_Subdiv);
  m_Meshsrc->Update();
  m_CoefsMesh = m_Meshsrc->GetOutput();
  m_CoefsMeshSpatialObject->SetMesh(m_CoefsMesh);
  m_CoefsMeshSpatialObject->ComputeBoundingBox();
  this->Update();
}

} // end namespace neurolib
