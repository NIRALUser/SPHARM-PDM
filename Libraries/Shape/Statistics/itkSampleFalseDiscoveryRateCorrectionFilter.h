/*=========================================================================

  SampleFalseDiscoveryRateCorrectionFilter

=========================================================================*/
#ifndef __SampleFalseDiscoveryRateCorrectionFilter_h
#define __SampleFalseDiscoveryRateCorrectionFilter_h

#include "itkListSample.h"
#include "itkSampleAlgorithmBase.h"

namespace neurolib{ 
namespace Statistics{
  
/** \class SampleFalseDiscoveryRateCorrectionFilter
 * \brief This filter corrects a sample vector of p-values using false discovery rate (FDR)
 *
 * This filter corrects a sample vector of p-values using false discovery rate.
 * The procedure has been taken from the 
 * Thresholding of Statistical Maps in Functional Neuroimaging Using the False Discovery Rate
 *    CR Genovese, NA Lazar, T Nichols, NeuroImage 15 2002
 * Usage: Set the acceptable false discovery rate Q (default 5%, SetFalseDicoveryRate) and the
 * number of steps N (default 1000, SetNumberOfSteps). The corrected p-values are determined at
 * the lowest level q_i where they pass the FDR threshold rate q_i = Q * (N - i )/ N 
 *
 * \sa 
 */

class SampleMeasurementVectorSizeException : public itk::ExceptionObject 
{
public:

  /** Run-time information. */
  itkTypeMacro( SampleMeasurementVectorSizeException , ExceptionObject );
  
  /** Constructor. */
  SampleMeasurementVectorSizeException(const char *file, unsigned int line, 
                           const char* message = "Bad Measurement Vector Size") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  SampleMeasurementVectorSizeException(const std::string &file, unsigned int line, 
                           const char* message = "Bad Measurement Vector Size") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }
};
  
template< class TSample >
class SampleFalseDiscoveryRateCorrectionFilter :
  public itk::Statistics::SampleAlgorithmBase< TSample >
{
public:
  /** Standard class typedefs. */
  typedef SampleFalseDiscoveryRateCorrectionFilter Self ;
  typedef itk::Statistics::SampleAlgorithmBase< TSample > Superclass ;
  typedef itk::SmartPointer<Self> Pointer ;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** Standard Macros */
  itkTypeMacro(SampleFalseDiscoveryRateCorrectionFilter, itk::Statistics::SampleAlgorithmBase);
  itkNewMacro(Self) ;
  

  typedef typename TSample::MeasurementVectorType MeasurementVectorType ;
  typedef typename MeasurementVectorType::ValueType 	MeasurementType;
  typedef typename Superclass::InputSampleType InputSampleType ;
  typedef typename itk::Statistics::ListSample< MeasurementVectorType > OutputSampleType ;

  /** Returns the correct p-value data in a ListSample object */
  typename OutputSampleType::Pointer GetOutput() ;

  /** Get/Set the highest acceptable false discovery rate */
  itkGetConstReferenceMacro(MaximumFalseDiscoveryRate, MeasurementType);  
  itkSetMacro(MaximumFalseDiscoveryRate, MeasurementType);
  
  /** Get/Set number of steps for computing the correction, represents the numerical resolution
      for the correction */
  itkGetConstReferenceMacro(NumberOfSteps, unsigned int);  
  itkSetMacro(NumberOfSteps, unsigned int);
  
protected:
  SampleFalseDiscoveryRateCorrectionFilter() ;
  virtual ~SampleFalseDiscoveryRateCorrectionFilter() ;
  void PrintSelf(std::ostream& os, itk::Indent indent) const;

  /** Do the FDR correction */
  void GenerateData() ;

private:
  typename OutputSampleType::Pointer m_Output ;

  /** The highest acceptable false discovery rate */
  MeasurementType m_MaximumFalseDiscoveryRate;
  
  /** Number of steps for computing the correction, represents the numerical resolution
      for the correction */
  unsigned int m_NumberOfSteps;

  
  
} ; // end of class
    
} // end of namespace Statistics 
} // end of namespace neurolib

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSampleFalseDiscoveryRateCorrectionFilter.txx"
#endif

#endif

