/** \class ParametricMeshToSPHARMSpatialObjectFilter.h
 *
 *  \brief This class provides a filter for the conversion of a parametrized Surface Mesh
 *         into a set of SPHARM coefficients.
 *
 *  \author Martin Styner
 */
#ifndef __namicParametricMeshToSPHARMSpatialObjectFilter_h
#define __namicParametricMeshToSPHARMSpatialObjectFilter_h

#include <itkMesh.h>
#include <itkTriangleCell.h>
#include <itkDefaultDynamicMeshTraits.h>
#include <itkProcessObject.h>
#include <itkVector.h>

#include "SphericalHarmonicSpatialObject.h"

/** \class ParametricMeshToSPHARMSpatialObjectFilter
 *
 *
 * \par
 * This class construct a 3D SPHARM (spherical harmonics descriptor) coefficients from
 * a parametrized mesh surface (2 meshes, a surface mesh and its corresponding
 * spherical parametrization mesh). The input of this filter is commonly computed using
 * BinaryMask3DEqualAreaParametricMeshSource.
 *
 * \par PARAMETERS
 *   FlipTemplate : Needed for the parameter alignment based on the first order ellipsoid (1 order SPHARM
 *   coefficient). When establishing correspondence with SPHARM, then it is essential to provide the same
 *   FlipTemplate to all surfaces that are compared.
 *
 * \par REFERENCE
 *  1. C. Brechbuhler, G. Gerig, and O. Kubler: Parametrization of closed surfaces for 3-D shape description,
 *     CVGIP: Image Under., vol. 61, pp.154 170, 1995.
 *
 * \par INPUT
 * The input should be 2 surface meshes: 1) a surface mesh representing the object
 * 2) a surface mesh representing the spherical parametrization. The parametrization is
 * optimal for encoding the SPHARM coefficients when it is preserving the area ratio of surface elements
 * from the object surface to the parametrization surface.
 *
 *  */

namespace neurolib
{

class ParametricMeshToSPHARMSpatialObjectFilterException : public itk::ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro( ParametricMeshToSPHARMSpatialObjectFilterException, ExceptionObject );

  /** Constructor. */
  ParametricMeshToSPHARMSpatialObjectFilterException(const char *file, unsigned int line,
                                                     const char* message = "Error in generating SPHARM coeffs") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  ParametricMeshToSPHARMSpatialObjectFilterException(const std::string & file, unsigned int line,
                                                     const char* message = "Error in generating SPHARM coeffs") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

};

class ParametricMeshToSPHARMSpatialObjectFilterLegendre
{
  const double       SQRT2;
  double * * *       plm;
  double       A(int l, int m);

  double       P(int l, int m, float z);

  void              ylm(int l, int m, float x, float y, float z, double & re, double & im);

public:
  const int L;
  double *  a;
  void              Ylm(int l, int m, float x, float y, float z, double & re, double & im);

  void              Ylm(int l, int m, float theta, float phi, double & re, double & im);

  ParametricMeshToSPHARMSpatialObjectFilterLegendre(const int L);
  ~ParametricMeshToSPHARMSpatialObjectFilterLegendre();
};

class ITK_EXPORT  ParametricMeshToSPHARMSpatialObjectFilter :
  public itk::ProcessObject
{
public:

  /** Standard "Self" typedef. */
  typedef ParametricMeshToSPHARMSpatialObjectFilter Self;
  typedef itk::ProcessObject                        Superclass;
  typedef itk::SmartPointer<Self>                   Pointer;
  typedef itk::SmartPointer<const Self>             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ParametricMeshToSPHARMSpatialObjectFilter, ProcessObject);

  /** Hold on to the type information specified by the template parameters. */
  typedef  itk::DataObject::Pointer                 DataObjectPointer;
  typedef  neurolib::SphericalHarmonicSpatialObject SpatialObjectType;
  typedef  SpatialObjectType::MeshType              InputMeshType;
  typedef  InputMeshType::Pointer                   InputMeshPointer;
  typedef  InputMeshType::MeshTraits                InputMeshTrait;
  typedef  InputMeshTrait::PointType                MeshPointType;
  typedef  InputMeshTrait::PixelType                MeshPixelType;

  /** Some convenient typedefs. */
  typedef  InputMeshType::CellTraits             CellTraits;
  typedef  InputMeshType::PointsContainerPointer PointsContainerPointer;
  typedef  InputMeshType::PointsContainer        PointsContainer;
  typedef  InputMeshType::CellsContainerPointer  CellsContainerPointer;
  typedef  InputMeshType::CellsContainer         CellsContainer;
  typedef  InputMeshType::PointType              PointType;
  typedef  InputMeshType::CellType               CellType;
  typedef  itk::TriangleCell<CellType>           TriangleType;

  typedef SpatialObjectType::CoefListType CoefListType;
  typedef SpatialObjectType::CoefType     CoefType;

  /** returns the Surface Mesh */
  virtual InputMeshType *  GetInputSurfaceMesh();

  /** sets the Surface Mesh */
  virtual void SetInputSurfaceMesh(InputMeshType * mesh);

  /** returns the Parametrization Mesh. Spherical vertices correspond to Surface vertices of same index */
  virtual InputMeshType * GetInputParametrizationMesh();

  /** sets the Parametrization Mesh. Spherical vertices correspond to Surface vertices of same index */
  virtual void SetInputParametrizationMesh(InputMeshType * mesh);

  SpatialObjectType * GetOutput();

  SpatialObjectType * GetOutput(unsigned int idx);

  CoefListType * GetEllipseAlignCoef();

  CoefListType * GetEllipseScaleAlignCoef();

  /** Degree of spherical harmonics */
  itkGetConstReferenceMacro(Degree, unsigned int);
  itkSetMacro(Degree, unsigned int);

  /** Whether to do the ellipse alignment of the parametrization */
  itkGetMacro(ParaEllipseAlignment, bool);
  itkSetMacro(ParaEllipseAlignment, bool);

  /** Fliptemplate in order to align the parametrization, without a flipTemplate it is difficult to
      compare different SPHARM objects, as the correspondence can be flipped at any of the axis (since
      the first order ellipsoid is ambiguosly defined */
  SpatialObjectType * GetFlipTemplate()
  {
    return m_FlipTemplate;
  }
  void SetFlipTemplate( SpatialObjectType * flipTemplate)
  {
    m_FlipTemplate = flipTemplate;
  }

  /** SetParametrizationFlipIndex allows an optional additional flipping of the parametrization.
      This should rarely be necessary, but is available for completeness.
      There are 7 possibilties and thus flipIndex has to be a value between 1 and 7
      0 = no additional flip, 1 = flip along x & y axis, 2 = flip along y & z, 3 = flip along x & z,
      4 = flip along x, 5 = flip along y,  6 = flip along x & y & z, 7 = flip at z
      where y is the smallest, x is the second smallest and z is the long axis of the ellipsoid
   */
  void SetParametrizationFlipIndex(int flipIndex)
  {
    m_FlipParametrizationIndex  = flipIndex;
  }

  using itk::ProcessObject::SetOutput;

  /** Set the mesh output of this process object. This call is slated
   * to be removed from ITK. You should GraftOutput() and possible
   * DataObject::DisconnectPipeline() to properly change the output. */
  void SetOutput(SpatialObjectType *output);

  /** Graft the specified DataObject onto this ProcessObject's output. */
  virtual void GraftOutput(SpatialObjectType *output);

  /** Make a DataObject of the correct type to used as the specified
   * output.  */
  using itk::ProcessObject::MakeOutput;

  virtual itk::ProcessObject::DataObjectPointer MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx) ITK_OVERRIDE;

  void GenerateOutputInformation() ITK_OVERRIDE
  {
  }
  void GenerateData() ITK_OVERRIDE;

protected:
  ParametricMeshToSPHARMSpatialObjectFilter();
  ~ParametricMeshToSPHARMSpatialObjectFilter() ITK_OVERRIDE;
  void PrintSelf(std::ostream& os, itk::Indent indent) const ITK_OVERRIDE;

  /**  Since all DataObjects should be able to set the requested region in
   * unstructured form, just copy output->RequestedRegion all inputs. */
  void GenerateInputRequestedRegion() ITK_OVERRIDE;

  /** compute best possible flip of coordinates and parameterization */
  int  BestFlipAlign(CoefListType * coef, const CoefListType * flipTemplate);

  void BestFlip(CoefListType * coef, const CoefListType * flipTemplate);

  /** adjust coefficients such that first order ellipsoid axis build a right hand coordinate system*/
  void MakeEllipsoidAxisRightHand(CoefListType * coef);

  double FlipAxesDistance(const CoefListType *coef1, const CoefListType * coef2);

  double CoefDistance(const CoefListType *coef1, const CoefListType * coef2);

  void FlipAlign(CoefListType * coef,  int  flipIndex);

  void EllipseAlign(const CoefListType * incoef, CoefListType * outcoef);

  void EllipseScaleAlign(const CoefListType * incoef, CoefListType * outcoef);

  void EllipseAlign(const CoefListType * incoef, CoefListType * outcoef, double scale);

  void ComputeCoeffs();

  /** Determine rotation matrix & rotate ParametricMesh */
  void RotateParametricMesh();

private:

  // flip normalization
  int flip_ll(double n);

  int flip_mm(double n);

  int flip_m2(double n);

  double flip_fu0(double n);

  double flip_fu1(double n);

  double flip_fu2(double n);

  double flip_reflect(double n);

  void flipu0(const CoefListType * incoef, CoefListType * outcoef);

  void flipu1(const CoefListType * incoef, CoefListType * outcoef);

  void flipu2(const CoefListType * incoef, CoefListType * outcoef);

  void reflect(const CoefListType * incoef, CoefListType * outcoef);

  int Get_BaseVal(PointsContainerPointer paraPoints, float * & A, int & m, int & n);

  unsigned int m_Degree;
  CoefListType m_coeffs;
  double *     m_flatCoeffs;

  bool m_ParaEllipseAlignment;
  int  m_FlipParametrizationIndex;

  SpatialObjectType::Pointer m_FlipTemplate;

  ParametricMeshToSPHARMSpatialObjectFilterLegendre * m_leg;

  /** Used by streaming: The requested region of the output being processed
   * by the execute method. Set in the GenerateInputRequestedRegion method. */
  int m_GenerateDataRegion;
  int m_GenerateDataNumberOfRegions;
};

}

#endif
