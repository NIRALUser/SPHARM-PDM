#ifndef _itkMesh3DProcrustesAlignFilter_h
#define _itkMesh3DProcrustesAlignFilter_h

#include "itkProcessObject.h"
#include "itkAffineTransform.h"
#include "itkTransformMeshFilter.h"

#include <fstream>

#if !defined(M_PI)
#define M_PI 3.14159265358979323846264338327950288   /* pi */
#endif

namespace itk
{
/** \class Mesh3DProcrustesAlignFilter
*   \brief Class for groupwise aligning a set of 3D meshes.
*
* All input meshes must have the same number of points (with corresponding
* indices).
* Default mode: Iterative Generalized Procrustes alignment with initialization
* from first object
* Option: UseInitialAverageOn/Off() Use  average structure as initialization (off by default)
*   Only appropriate if objects are already pre-aligned
* Option: UseSingleIterationOn/Off() Only run one iteration (off by default)
*
* All output meshes will be centered at the origin and scaled to
* match a mean shape with norm 1.
* Option: UseNormalizationOn/Off() disables normalization of scaling and centering (on by default)
* Option: UseScalingOn/Off() disables scaling matching with Procrustes (does not disable scaling normalization of the mean shape) (on by default)
*
* GetTransform() can be used to query the employed transformation for each
* input mesh.
*
* \author Tobias Heimann. Division Medical and Biological Informatics,
*         German Cancer Research Center, Heidelberg, Germany.
* change  Martin Styner, UNC support for single template and initialization with average
* TODO: Enable/Disable Normalization of centering and scaling to origin
*
*/
template <class TInputMesh, class TOutputMesh>
class Mesh3DProcrustesAlignFilter : public ProcessObject
{

public:

  /** Standard typedefs. */
  typedef Mesh3DProcrustesAlignFilter Self;
  typedef ProcessObject               Superclass;
  typedef SmartPointer<Self>          Pointer;
  typedef SmartPointer<const Self>    ConstPointer;

  /** Convenient typedefs. */
  typedef TInputMesh                                                         InputMeshType;
  typedef TOutputMesh                                                        OutputMeshType;
  typedef typename InputMeshType::Pointer                                    InputMeshPointer;
  typedef typename OutputMeshType::Pointer                                   OutputMeshPointer;
  typedef typename InputMeshType::PointType                                  InputPointType;
  typedef typename OutputMeshType::PointType                                 OutputPointType;
  typedef typename OutputMeshType::PointsContainer                           OutputPointsContainer;
  typedef typename OutputMeshType::PointsContainerPointer                    OutputPointsContainerPointer;
  typedef DataObject::Pointer                                                DataObjectPointer;
  typedef typename OutputMeshType::CoordRepType                              CoordRepType;
  typedef vnl_matrix<CoordRepType>                                           MatrixType;
  typedef AffineTransform<CoordRepType, 3>                                   TransformType;
  typedef typename TransformType::Pointer                                    TransformPointer;
  typedef TransformMeshFilter<OutputMeshType, OutputMeshType, TransformType> TransformMeshType;
  typedef typename TransformMeshType::Pointer                                TransformMeshPointer;
  typedef std::vector<TransformMeshPointer>                                  TransformMeshArray;
  typedef typename TransformType::OffsetType                                 TranslationType;
  typedef typename TransformType::MatrixType                                 RotationType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(Mesh3DProcrustesAlignFilter, ProcessObject);

  /** Sets the number of input meshes that have to be aligned.
  * Call this before any other methods.
  */
  void SetNumberOfInputs( unsigned int num );

  using itk::ProcessObject::SetInput;

  /** Sets one input mesh (starting with idx=0). */
  void SetInput( unsigned int idx, InputMeshPointer mesh );

  /** Returns one of the input meshes (starting with idx=0). */
  InputMeshType * GetInput( unsigned int idx )
  {
    return static_cast<InputMeshType *>(this->ProcessObject::GetInput( idx ) );
  }

  /** Gets one transformed output mesh (starting with idx=0). */
  OutputMeshType * GetOutput(unsigned int idx);

  /** Gets the transformed mean of all aligned meshes. */
  itkGetConstObjectMacro( Mean, OutputMeshType );

  /** Returns the transform used to align a mesh (starting with idx=0). */
  TransformType * GetTransform( unsigned int idx )
  {
    return m_MeshTransform[idx]->GetTransform();
  }

  // bp2009
  TranslationType GetRotationDegrees( TransformType* transform )
  {
    TranslationType rotationXYZ = transform->GetTranslation();
    RotationType    rotation = transform->GetMatrix();
    double          rotX, cosY, rotY, rotZ;

    rotY = (-1) * (asin(rotation(0, 2) ) ); // C++ functions need radians
    rotY = rotY * ( 180.0 / M_PI );         // I should have degrees here
    cosY = cos(rotY * M_PI / 180.0);        // I have radians here
    rotX = acos( (rotation(2, 2) ) / cosY);
    rotX = rotX * ( 180.0 / M_PI ); // I should have degrees here
    rotZ = acos( (rotation(0, 0) ) / cosY);
    rotZ = rotZ * ( 180.0 / M_PI ); // I should have degrees here
    // std::cout << rotX << std::endl;
    // std::cout << rotY << std::endl;
    // std::cout << rotZ << std::endl;

    rotationXYZ.SetElement(0, rotX);
    rotationXYZ.SetElement(1, rotY);
    rotationXYZ.SetElement(2, rotZ);

    return rotationXYZ;
  }

  /** - IO Function to get a file with the Rotation Matrix + Translational Vector. */
  void PrintTransform( TransformType* transform )
  {
    TranslationType trans = transform->GetTranslation();
    TranslationType rots = GetRotationDegrees(transform);
    RotationType    rotation = transform->GetMatrix();
    // std::cout << "Shape " << i << std::endl;
    std::ofstream output;

    output.open("Transform.info", std::ios::out);
    output << "ROTATION MATRIX" << std::endl;
    for( unsigned dim1 = 0; dim1 < 3; dim1++ )
      {
      for( unsigned dim2 = 0; dim2 < 3; dim2++ )
        {
        output << "Position " << dim1 << " " << dim2 << " : " << rotation(dim1, dim2) << std::endl;
        }
      }

    output << "TRANSLATION VECTOR" << std::endl;
    for( unsigned dim1 = 0; dim1 < 3; dim1++ )
      {
      output << "Position " << dim1 << " : " << trans.GetElement(dim1) << std::endl;
      }

    output << "ROTATION VECTOR" << std::endl;
    for( unsigned dim1 = 0; dim1 < 3; dim1++ )
      {
      output << "Position " << dim1 << " : " << rots.GetElement(dim1) << std::endl;
      }

    output.close();
  }

  // bp2009

  /** Get/Set the convergence value that determines when the iterative
  * calculation of alignment is stopped. The smaller the value, the higher
  * the accuracy (the default should be sufficient for all applications).
  */
  itkGetConstMacro( Convergence, double );
  itkSetMacro(Convergence, double );

  /** Creates an ouput object. */
  using itk::ProcessObject::MakeOutput;

  virtual itk::ProcessObject::DataObjectPointer
  MakeOutput(itk::ProcessObject::DataObjectPointerArraySizeType idx);

  /** Normalization with Scaling on (default)*/
  void SetUseScalingOn()
  {
    this->SetUseScaling(true);
  }

  /** Normalization with Scaling off */
  void SetUseScalingOff()
  {
    this->SetUseScaling(false);
  }

  /** Set/Get whether or not the filter will use Scaling normalization or not */
  itkSetMacro(UseScaling, bool);
  itkGetMacro(UseScaling, bool);

  /** use average as reference for initialization */
  void SetUseInitialAverageOn()
  {
    this->SetUseInitialAverage(true);
  }

  /** do not use average as reference for initialization. The first surface will be used instead.(default) */
  void SetUseInitialAverageOff()
  {
    this->SetUseInitialAverage(false);
  }

  /** Set/Get whether or not the average is used as initialization */
  itkSetMacro(UseInitialAverage, bool);
  itkGetMacro(UseInitialAverage, bool);

  /** Set/Get whether or not the centering and Scaling to a shape with norm 1 is performed */
  void SetUseNormalizationOn()
  {
    this->SetUseNormalization(true);
  }

  void SetUseNormalizationOff()
  {
    this->SetUseNormalization(false);
  }

  itkSetMacro(UseNormalization, bool);
  itkGetMacro(UseNormalization, bool);

  /** Only do one iteration */
  void SetUseSingleIterationOn()
  {
    this->SetUseSingleIteration(true);
  }

  /** Iterate until convergence (default)*/
  void SetUseSingleIterationOff()
  {
    this->SetUseSingleIteration(false);
  }

  /** Set/Get whether or not only one iteration should be run */
  itkSetMacro(UseSingleIteration, bool);
  itkGetMacro(UseSingleIteration, bool);

  /** Set/Get wheter rotation should be performed or not */
  void SetAlignRotationOn()
  {
    this->SetAlignRotation(true);
  }

  void SetAlignRotationOff()
  {
    this->SetAlignRotation(false);
  }

  itkSetMacro(AlignRotation, bool );
  itkGetMacro(AlignRotation, bool );
protected:

  /** Standard constructor. */
  Mesh3DProcrustesAlignFilter();

  /** Standard destructor. */
  ~Mesh3DProcrustesAlignFilter();

  /** Performs the alignment. */
  virtual void GenerateData();

  /** Calculates the center coordinates of the specified input mesh. */
  TranslationType GetMeshCenter( unsigned int idx );

  /** Uses the current transformations to calculate a new mean.*/
  void CalculateMean();

  /** Returns the best transform to fit input mesh idx, translated to
   * zero origin by using the values in m_Center, to the given targetMesh.
   * targetMesh has to be centered at zero origin as well!
   */
  TransformPointer GetProcrustesMatch( unsigned int idx, OutputMeshPointer targetMesh, TranslationType targetCenter );

private:

  /** When the difference between two consecutive mean calculations gets
  * samller than m_Convergence, the alignment is terminated. */
  double m_Convergence;
  /** Holds the transforms for all input meshes.*/
  TransformMeshArray m_MeshTransform;
  /** Holds the center coordinates for all input meshes.*/
  std::vector<TranslationType> m_Center;
  /** The mean of all transformed meshes. */
  OutputMeshPointer m_Mean;
  TranslationType   m_MeanCenter;
  /** The mean of all transformed meshes from the preceding iteration. */
  OutputMeshPointer m_OldMean;

  /** Scaling on/off */
  bool m_UseScaling;

  /** Normalization on/off: Centering and Scaling to Unit Length */
  bool m_UseNormalization;

  /** Use Average for initialization on/off */
  bool m_UseInitialAverage;

  /** Only run one single iteration on/off */
  bool m_UseSingleIteration;

  /** Perform rotation on/off */
  bool m_AlignRotation;

};

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMesh3DProcrustesAlignFilter.txx"
#endif

#endif
