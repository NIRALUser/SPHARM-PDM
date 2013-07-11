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
  this->SetDimension(SphericalHarmonicSpatialObjectDimension); // Dimension 3
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
bool SphericalHarmonicSpatialObject::ValueAt( const PointType & /* point */,
                                              double & /* value */,
                                              unsigned int /* depth */,
                                              char * /* name */) const
{
  return false;
}

bool SphericalHarmonicSpatialObject::IsEvaluableAt(const PointType & /* point */,
                                                   unsigned int /* depth */,
                                                   char * /* name */) const
{
  return false;
}

bool SphericalHarmonicSpatialObject::IsInside(const PointType & point, unsigned int depth, char *name) const
{
  itkDebugMacro( "Checking the point [" << point << "] is inside the SphericalHarmonicSpatialObject" );
  if( name == NULL )
    {
    if( IsInside(point) )
      {
      return true;
      }
    }
  else if( strstr(typeid(Self).name(), name) )
    {
    if( IsInside(point) )
      {
      return true;
      }
    }
  return Superclass::IsInside(point, depth, name);
}

bool SphericalHarmonicSpatialObject::IsInside(const PointType & point) const
{
  if( !this->GetIndexToWorldTransform()->GetInverse(const_cast<TransformType *>(this->GetInternalInverseTransform() ) ) )
    {
    return false;
    }

  PointType transformedPoint = this->GetInternalInverseTransform()->TransformPoint(point);
  // std::cout<<" transformedPoint = [ " << transformedPoint[0] <<" , " <<transformedPoint[1] <<" , "
  // <<transformedPoint[2] << " ] "<<std::endl;

  this->ComputeLocalBoundingBox();

  if( this->GetBounds()->IsInside(point) )
    {
    return m_CoefsMeshSpatialObject->IsInside(transformedPoint);
    }
  return false;
}

bool SphericalHarmonicSpatialObject::ComputeLocalBoundingBox() const
{
  itkDebugMacro( "Computing SphericalHarmonicSpatialObject local bounding box" );

  if( this->GetBoundingBoxChildrenName().empty()
      || strstr(typeid(Self).name(), this->GetBoundingBoxChildrenName().c_str() ) )
    {
    m_CoefsMeshSpatialObject->ComputeLocalBoundingBox();

    PointType ptMin, ptMax;
    ptMin = m_CoefsMeshSpatialObject->GetBoundingBox()->GetMinimum();
    ptMin = this->GetIndexToWorldTransform()->TransformPoint(ptMin);
    ptMax = m_CoefsMeshSpatialObject->GetBoundingBox()->GetMaximum();
    ptMax = this->GetIndexToWorldTransform()->TransformPoint(ptMax);

    const_cast<BoundingBoxType *>(this->GetBounds() )->SetMinimum(ptMin);
    const_cast<BoundingBoxType *>(this->GetBounds() )->SetMaximum(ptMax);
    }
  return true;
}

void SphericalHarmonicSpatialObject::ComputeHiddenMeshSpatialObject()
{
  SphericalHarmonicMeshSource::Pointer m_Meshsrc = SphericalHarmonicMeshSource::New();

  m_Meshsrc->SetCoefs(m_Coefs);
  m_Meshsrc->SetLevel(m_Subdiv);
  m_Meshsrc->Update();
  m_CoefsMesh = m_Meshsrc->GetOutput();
  m_CoefsMeshSpatialObject->SetMesh(m_CoefsMesh);
  m_CoefsMeshSpatialObject->ComputeLocalBoundingBox();

}

} // end namespace neurolib
