#ifndef __namicSphericalHarmonicSpatialObject_h

#define __namicSphericalHarmonicSpatialObject_h

/** \class SphericalHarmonicSpatialObject
 *
 *  \brief This class provides 3D shape representation using spherical harmonic coefficients.
 *
 *  \author Christine Xu
 */

#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkSpatialObject.h"
#include "itkMeshSpatialObject.h"

#include <vector>

#define SphericalHarmonicSpatialObjectDimension 3

namespace neurolib
{

class SphericalHarmonicSpatialObject : public itk::SpatialObject<SphericalHarmonicSpatialObjectDimension>
{
public:

  typedef SphericalHarmonicSpatialObject                                  Self;
  typedef itk::SpatialObject<SphericalHarmonicSpatialObjectDimension>     Superclass;
  typedef itk::SmartPointer<Self>                                         Pointer;
  typedef itk::SmartPointer<const Self>                                   ConstPointer;
  typedef itk::SmartPointer<Superclass>                                   SuperclassPointer;
  typedef double                                                          ScalarType;
  typedef itk::Point<ScalarType, SphericalHarmonicSpatialObjectDimension> CoefType;
  typedef std::vector<CoefType>                                           CoefListType;
  typedef Superclass::PointType                                           PointType;
  typedef Superclass::TransformType                                       TransformType;
  typedef Superclass::BoundingBoxType                                     BoundingBoxType;

  /** Traits for Mesh, DefaultDynamicMeshTraits<TPixelType, VPointDimension, VMaxTopologicalDimension, TCoordRep> */
  typedef itk::DefaultDynamicMeshTraits<float, SphericalHarmonicSpatialObjectDimension,
                                        SphericalHarmonicSpatialObjectDimension, double> MeshTrait;
  typedef itk::Mesh<float, SphericalHarmonicSpatialObjectDimension,
                    MeshTrait>                                                           MeshType;
  typedef itk::MeshSpatialObject<MeshType>
  MeshSpatialObjectType;
  typedef MeshSpatialObjectType::Pointer
  MeshSpatialObjectPointer;

  // typedef typename SphericalHarmonicMeshSource MeshSourceType;

  itkStaticConstMacro(NumberOfDimension, unsigned int, SphericalHarmonicSpatialObjectDimension);

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Method for creation through the object factory. */
  itkTypeMacro( SphericalHarmonicSpatialObject, SpatialObject );

  /** Get the count */
  itkGetConstReferenceMacro(Count, int);

  /** Get the harmonic */
  itkGetConstReferenceMacro(Harmonic, int);

  /** Get the subdiv */
  itkGetConstReferenceMacro(Subdiv, int);

  /** Get the theta */
  // itkGetConstReferenceMacro(Theta, int);

  /** Get the phi */
  // itkGetConstReferenceMacro(Phi, int);

  /** Get the vertex */
  // itkGetConstReferenceMacro(Vertex, int);

  /** Get the triangs */
  // itkGetConstReferenceMacro(Triangs, int);

  /** Get the coefficients */
  void GetCoefs(CoefListType& coeflist) const; // Deep copy of the m_Coefs

  void SetCoefs(CoefListType& coeflist);

  /** Derived spatial geometric functions from the superclass SpatialObject,
      resulting directly from calling corresponding functions of p_coefsMeshSpatialObject. */
  virtual bool  ValueAt(const PointType & point, double & value, unsigned int depth = 0, char *name = NULL) const;

  virtual bool  IsEvaluableAt(const PointType & point, unsigned int depth = 0, char *name = NULL) const;

  virtual bool  IsInside(const PointType & point, unsigned int depth, char *name = NULL) const;

  virtual bool  IsInside(const PointType & point) const;

  // virtual bool  ComputeBoundingBox () const;
  virtual bool ComputeLocalBoundingBox() const;

  virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

/** The four vitual geometric functions calls this function whenever the b_Dirty is true. */
  void ComputeHiddenMeshSpatialObject();

protected:

  SphericalHarmonicSpatialObject();
  virtual ~SphericalHarmonicSpatialObject();

  /** Set functions */
  itkSetMacro(Count, int);
  itkSetMacro(Harmonic, int);
  // itkSetMacro(Theta, int);
  // itkSetMacro(Phi, int);
  // itkSetMacro(Vertex, int);
  // itkSetMacro(Triangs, int);

  /** List of coefficients. */
  CoefListType m_Coefs;
  /** Number of the x,y,z coefficients in total. */

  int m_Count;
  /** Maximal degree of the spherical harmonics. */

  int m_Harmonic;

  /** Level of icosahedral subdivision of the background MeshSpatialObject(p_coefsMeshSpatialObject). */
  int m_Subdiv;
  /** Number of theta */
  // int m_Theta;
  /** Number of phi */
  // int m_Phi;
  /** Number of vertices */
  // int m_Vertex;
  /** Number of triangles */
  // int m_Triangs;
  MeshType::Pointer m_CoefsMesh;
  /** This MeshSpatialObject is hidden in the background
      and is used for providing results for the above geometric functions. */
  MeshSpatialObjectPointer m_CoefsMeshSpatialObject;

};

} // end namespace neurolib

#endif // __namicSphericalHarmonicSpatialObject_h
