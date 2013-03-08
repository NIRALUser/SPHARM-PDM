/** \class BinaryMask3DEqualAreaParametricMeshSource
 *
 *  \brief This class provides a mesh source for both a surface mesh and the corresponding Spherical Parametrization Mesh
 *
 *  \author Martin Styner
 */
#ifndef __namicBinaryMask3DEqualAreaParametricMeshSource_h
#define __namicBinaryMask3DEqualAreaParametricMeshSource_h

#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageToMeshFilter.h>
#include <itkTriangleCell.h>
#include <itkDefaultDynamicMeshTraits.h>
#include <itkProcessObject.h>
#include <math.h>

#include <itkIndex.h>

/** \class BinaryMask3DEqualAreaParametricMeshSource
 *
 *
 * \par
 * This class construct a 3D mesh surface and its corresponding Spherical parametrization
 * based on a binary mask. The binary mask has thus to be of spherical topology.
 * The parametrization is an optimization process that leads to area preservation
 * and minimization of angular distortion.
 *
 *
 * \par PARAMETERS
 * ObjectValue: The ObjectValue parameter is used to identify the object. In most applications,
 *  pixels in the object region are assigned to "1", so the default value of ObjectValue is
 *  set to "1"
 * NumberOfIterations: The number of iterations for the computation of the parametrization
 *  for small object (<1000 vertices) this can be in the order of 100-200 (have a cup of tea)
 *  for larger object (>1000 vertices) the default (500) or up to 1000 iterations is better
 *    (visit the closest coffe place)
 *  for hug objects (>100'000 vertices) this should be at least 5000 (go home and let it run overnight...)
 *
 * InitParametricMesh: Optional initialization, can be used to continue on an earlier run of this filter
 *
 * \par REFERENCE
 *  1. C. Brechbuhler, G. Gerig, and O. Kubler: Parametrization of closed surfaces for 3-D shape description,
 *     CVGIP: Image Under., vol. 61, pp.154 170, 1995.
 *  2. C. Brechbuhler: Description and Analysis of 3-D Shapes by Parametrization of Closed Surfaces, 1995,
 *     Diss., IKT/BIWI, ETH Z¨urich, ISBN 3-89649-007-9.
 *
 * \par INPUT
 * The input should be a 3D binary image of spherical topology with isotropic spacing
 * Non-isotropic spacing will lead to parametrization that does NOT preserve the aral
 *
 *  */

namespace itk
{

class BinaryMask3DEqualAreaParametricMeshSourceException : public ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro(BinaryMask3DEqualAreaParametricMeshSourceException, ExceptionObject );

  /** Constructor. */
  BinaryMask3DEqualAreaParametricMeshSourceException(const char *file, unsigned int line, const char* message =
                                                       "Error in generating parametrized mesh from Binary Mask") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  BinaryMask3DEqualAreaParametricMeshSourceException(const std::string & file, unsigned int line, const char* message =
                                                       "Error in generating parametrized mesh from Binary Mask") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

};

#define EqualAreaParametricMeshSourceMeshTypeDimension 3

typedef DefaultDynamicMeshTraits<float, EqualAreaParametricMeshSourceMeshTypeDimension,
                                 EqualAreaParametricMeshSourceMeshTypeDimension,
                                 double> EqualAreaParametricMeshSourceMeshTrait;
typedef Mesh<float, EqualAreaParametricMeshSourceMeshTypeDimension,
             EqualAreaParametricMeshSourceMeshTrait>
BinaryMask3DEqualAreaParametricMeshSourceMeshType;

template <class TInputImage>
class ITK_EXPORT  BinaryMask3DEqualAreaParametricMeshSource :
  public ImageToMeshFilter<TInputImage, BinaryMask3DEqualAreaParametricMeshSourceMeshType>
{
public:

  /** Standard "Self" typedef. */
  typedef BinaryMask3DEqualAreaParametricMeshSource                                         Self;
  typedef ImageToMeshFilter<TInputImage, BinaryMask3DEqualAreaParametricMeshSourceMeshType> Superclass;
  typedef SmartPointer<Self>                                                                Pointer;
  typedef SmartPointer<const Self>                                                          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BinaryMask3DEqualAreaParametricMeshSource, ImageToMeshFilter);

  /** Hold on to the type information specified by the template parameters. */
  typedef BinaryMask3DEqualAreaParametricMeshSourceMeshType OutputMeshType;
  typedef typename OutputMeshType::Pointer                  OutputMeshPointer;
  typedef typename OutputMeshType::MeshTraits               OutputMeshTrait;
  typedef typename OutputMeshTrait::PointType               MeshPointType;
  typedef typename OutputMeshTrait::PixelType               MeshPixelType;

  /** Some convenient typedefs. */
  typedef typename OutputMeshType::Pointer                MeshPointer;
  typedef typename OutputMeshType::CellTraits             CellTraits;
  typedef typename OutputMeshType::PointsContainerPointer PointsContainerPointer;
  typedef typename OutputMeshType::PointsContainer        PointsContainer;
  typedef typename OutputMeshType::CellsContainerPointer  CellsContainerPointer;
  typedef typename OutputMeshType::CellsContainer         CellsContainer;
  typedef typename OutputMeshType::PointType              PointType;
  typedef typename OutputMeshType::CellType               CellType;
  typedef TriangleCell<CellType>                          TriangleType;

  typedef CovariantVector<double, 2> doubleVector;
  typedef CovariantVector<int, 2>    intVector;

  /** Input Image Type Definition. */
  typedef TInputImage                           InputImageType;
  typedef typename InputImageType::Pointer      InputImagePointer;
  typedef typename InputImageType::ConstPointer InputImageConstPointer;
  typedef typename InputImageType::PixelType    InputPixelType;
  typedef typename InputImageType::SpacingType  InputSpacingType;
  typedef typename InputImageType::PointType    InputOriginType;
  typedef typename InputImageType::SizeType     InputSizeType;

  /** Type definition for the classified image index type. */
  typedef typename InputImageType::IndexType InputImageIndexType;

  itkGetConstReferenceMacro(ObjectValue, InputPixelType);
  itkSetMacro(ObjectValue, InputPixelType);

  /** Euler Number of mesh  */
  itkGetConstReferenceMacro(EulerNum,  int);
  itkSetMacro(EulerNum, int);

  /** Number of Iterations of Parametrization optimization */
  itkGetConstReferenceMacro(NumberOfIterations, unsigned int);
  itkSetMacro(NumberOfIterations, unsigned int);

  using itk::ProcessObject::SetInput;
  /** accept the input image */
  virtual void SetInput( const InputImageType * inputImage );

  /** returns the Surface Mesh: same as GetOutput(0) */
  virtual  OutputMeshType *  GetSurfaceMesh();

  /** returns the Parametrization Mesh: same as GetOutput(1). Forms a sphere,
      Spherical vertices correspond to Surface vertices of same index */
  virtual OutputMeshType * GetParametrizationMesh();

  /** Parametrization Mesh as starting value for optimization */
  virtual OutputMeshType * GetInitParametricMesh()
  {
    return m_InitParametricMesh.GetPointer();
  }

  virtual void SetInitParametricMesh( OutputMeshType * in)
  {
    m_InitParametricMesh = in; m_InitParametricMeshSet = true;
  };
  virtual void ClearInitParametricMesh()
  {
    m_InitParametricMesh = NULL; m_InitParametricMeshSet = false;
  };
protected:
  BinaryMask3DEqualAreaParametricMeshSource();
  ~BinaryMask3DEqualAreaParametricMeshSource();
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

  virtual void GenerateOutputInformation()
  {
  };                                          // do nothing
private:

  typedef Image<unsigned char, 3>            CharImageType;
  typedef typename CharImageType::RegionType CharImageRegion;
  typedef typename CharImageType::IndexType  CharImageIndex;
  typedef ImageRegionIterator<CharImageType> CharImageIterator;

  BinaryMask3DEqualAreaParametricMeshSource(const Self &); // purposely not implemented
  void operator=(const Self &);                            // purposely not implemented

  unsigned int      m_NumberOfIterations;
  int               m_EulerNum;
  InputPixelType    m_ObjectValue;
  OutputMeshPointer m_InitParametricMesh;
  bool              m_InitParametricMeshSet;
};
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "BinaryMask3DEqualAreaParametricMeshSource.txx"
#endif

#endif
