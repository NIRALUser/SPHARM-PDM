#ifndef _namicSphericalHarmonicMeshSource_h
#define _namicSphericalHarmonicMeshSource_h

/** \class SphericalHarmonicMeshSource
 *
 *  \brief This class provides convenient class for generating 3D meshes represented by Spherical Harmonics descriptor.
 *
 *  \author Christine Xu
 */

#include "SphericalHarmonicPolynomial.h"
#include "SphericalHarmonicSpatialObject.h"

#include "itkMeshSource.h"
#include "itkMesh.h"

namespace neurolib
{

class SphericalHarmonicMeshSourceException : public itk::ExceptionObject 
{
public:
  /** Run-time information. */
  itkTypeMacro( SphericalHarmonicMeshSourceException, ExceptionObject );
  
  /** Constructor. */
  SphericalHarmonicMeshSourceException(const char *file, unsigned int line, 
                           const char* message = "Error in generating meshes from Spherical Harmonics") : 
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  SphericalHarmonicMeshSourceException(const std::string &file, unsigned int line, 
                           const char* message = "Error in generating meshes from Spherical Harmonics") : 
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }
};

typedef SphericalHarmonicSpatialObject::MeshType  SphericalHarmonicMeshSourceMeshType;
//default dimension for mesh is 3
class SphericalHarmonicMeshSource : public itk::MeshSource< SphericalHarmonicMeshSourceMeshType>
{
public:
  typedef SphericalHarmonicMeshSource Self; 
  typedef SphericalHarmonicMeshSourceMeshType MeshType; 
  typedef itk::MeshSource<MeshType>  Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;
  
  /** Method for creation through the object factory. */
  itkNewMacro(Self);  

  /** Run-time type information (and related methods). */
  itkTypeMacro(SphericalHarmonicMeshSource, MeshSource);
  
  typedef SphericalHarmonicSpatialObject::ScalarType ScalarType;
  typedef SphericalHarmonicSpatialObject::CoefType CoefType;
  typedef SphericalHarmonicSpatialObject::CoefListType CoefListType;
  
  typedef double Point3[3];
  
  void SetCoefs(CoefListType& coeflist);
  typedef SphericalHarmonicSpatialObject::MeshTrait MeshTrait;
  typedef MeshType::PointsContainerPointer PointsContainerPointer;
  typedef MeshType::PointsContainer   PointsContainer;
  typedef MeshType::PointType PointType;
  typedef MeshType::CellsContainerPointer CellsContainerPointer;
  typedef MeshType::CellsContainer CellsContainer;
  typedef MeshType::CellType CellType;  
  
  itkGetConstReferenceMacro(Level, unsigned int);  
  itkSetMacro(Level, unsigned int);
  
  itkGetConstReferenceMacro(Degree, unsigned int);  
  void SetDegree(unsigned int d);
  
  itkGetConstReferenceMacro(FromL, unsigned int);  
  itkSetMacro(FromL, unsigned int);
  
  itkGetConstReferenceMacro(ToL, unsigned int);  
  itkSetMacro(ToL, unsigned int);
  
  itkGetConstReferenceMacro(Dimension, unsigned int);

  MeshType * GetOutputParaMesh()
    { return outputParaMesh; }

  double * GetPrecomputedIcosList ()
    { return m_icos ; }
  
protected:
  SphericalHarmonicMeshSource();
  ~SphericalHarmonicMeshSource();
  
  void GenerateData();
  
  void set_up_icosahedron_triangs(Point3* all_vert,
        Point3* all_triangs,
        int subdiv,
        int n_vert,
        int n_phi,
        int n_theta,
        double *icos,
        int n_triangs,
        int *triangs);        
  void interpol_vert(int n_phi,
        int n_theta,
        Point3 *mesh,
        int n_vertex,
        double *icos,
        Point3 *vertex);
  int mod(int x, int y);
  void interpol_2d(int n_phi,
      int n_theta,
      int xd,
      int xu,
      int yd,
      int yu,
      double ksi,
      double eta,
      Point3 *mesh,
      double *xi,
      double *yi,
      double *zi);

private:
  MeshType::Pointer outputParaMesh;
  
  unsigned int m_Dimension;//dimension of the output mesh

  unsigned int m_Level; //subdivision level
    
  unsigned int m_Degree; //degree of the coefficients
  CoefListType m_Coefs;
  
  unsigned int m_FromL;
  unsigned int m_ToL;

  double *m_icos ;  // icosahedron
};

}

#endif
