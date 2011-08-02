/*=========================================================================

  SampleFalseDiscoveryRateCorrectionFilter.txx

  Martin Styner
 This filter corrects a sample vector of p-values using false discovery rate.
 The procedure has been taken from the
 Thresholding of Statistical Maps in Functional Neuroimaging Using the False Discovery Rate
     CR Genovese, NA Lazar, T Nichols, NeuroImage 15 2002
  Usage: Set the acceptable false discovery rate Q (default 5%, SetFalseDicoveryRate) and the
  number of steps N (default 1000, SetNumberOfSteps). The corrected p-values are determined at
  the lowest level q_i where they pass the FDR threshold rate q_i = Q * (N - i )/ N


=========================================================================*/
#ifndef __itkSampleFalseDiscoveryRateCorrectionFilter_txx
#define __itkSampleFalseDiscoveryRateCorrectionFilter_txx

#include <algorithm>

namespace neurolib
{
namespace Statistics
{

template <class TSample>
SampleFalseDiscoveryRateCorrectionFilter<TSample>
::SampleFalseDiscoveryRateCorrectionFilter()
{
  m_Output = OutputSampleType::New();
  m_MaximumFalseDiscoveryRate = 0.05;
  m_NumberOfSteps = 1000;
}

template <class TSample>
SampleFalseDiscoveryRateCorrectionFilter<TSample>
::~SampleFalseDiscoveryRateCorrectionFilter()
{
}

template <class TSample>
void
SampleFalseDiscoveryRateCorrectionFilter<TSample>
::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "FDR:MaximumFalseDiscoveryRate: " << m_MaximumFalseDiscoveryRate << std::endl;
  os << indent << "FDR:NumberOfSteps: " << m_NumberOfSteps << std::endl;
  os << indent << "FDR:Output: " << m_Output << std::endl;
}

template <class TSample>
typename SampleFalseDiscoveryRateCorrectionFilter<TSample>::OutputSampleType::Pointer
SampleFalseDiscoveryRateCorrectionFilter<TSample>
::GetOutput()
{
  return m_Output;
}

// VS6 cannot compile the following method without the following pragma
#if defined(_MSC_VER)
#pragma inline_depth(0)
#endif

template <class TSample>
void
SampleFalseDiscoveryRateCorrectionFilter<TSample>
::GenerateData()
{
  // typename InputSampleType::ConstIterator iterInput = this->GetInputSample()->Begin() ;
  // typename InputSampleType::ConstIterator endInput = this->GetInputSample()->End() ;
  typename Superclass::ConstIterator iterInput = this->GetInputSample()->Begin();
  typename Superclass::ConstIterator endInput = this->GetInputSample()->End();

  if( this->GetMeasurementVectorSize() > 1 )
    {
    // pvalues are always single dimensional, thus raise error is this is not the case
    throw SampleMeasurementVectorSizeException(__FILE__, __LINE__, "Measurement vector size != 1!");
    }

  m_Output->SetMeasurementVectorSize( this->GetMeasurementVectorSize() );
  m_Output->Clear();

  // Copy the Samples and initialize output to all 1
  typedef typename std::vector<MeasurementType> sortSampleType;
  sortSampleType sortSample;
  sortSample.clear();

  // std::cout << "number of samples " << this->GetInputSample()->Size()  << std::endl;

  MeasurementVectorType mv;
  while( iterInput != endInput )
    {
    sortSample.push_back( (iterInput.GetMeasurementVector() )[0] );  // take only the first element of the vector
    mv[0] = 1;
    m_Output->PushBack(mv);   // initialize output to 1
    ++iterInput;
    }

  // Sort the Copied Samples
  std::sort(sortSample.begin(), sortSample.end() );

  // compute Threshold
  // Threshold = i / N * q

  double       incrementRate = m_MaximumFalseDiscoveryRate / m_NumberOfSteps;
  unsigned int sampleSize = sortSample.size();
  for( double currentRate =  m_MaximumFalseDiscoveryRate; currentRate > 0; currentRate -= incrementRate )
    {
    MeasurementType thresholdValue = 0;
    for( unsigned int index = 0; index < sampleSize; index++ )
      {
      double currentValue = sortSample[index];
      // as the FDR correction start at 'index' 1, we have to add one to the C++ index
      if( currentValue <= (double) (index + 1) / sampleSize * currentRate )
        {
        thresholdValue = currentValue;
        }
      }
    // std::cout << currentRate << "/" << thresholdValue << ",";

    // Threshold the input data
    typename OutputSampleType::ConstIterator iterOutput = m_Output->Begin();

    iterInput = this->GetInputSample()->Begin();
    while( iterInput != endInput )
      {
      if( thresholdValue >= (iterInput.GetMeasurementVector() )[0] )
        {
        m_Output->SetMeasurement(iterOutput.GetInstanceIdentifier(), 0, currentRate);
        }
      ++iterOutput;
      ++iterInput;
      }
    }

  // std::cout << std::endl;

  // debug output
  /*
  {
    for (unsigned int index = 0 ; index < sampleSize ; index++)
    {
      std::cout << sortSample[index] << ", ";
    }
    std::cout << std::endl;

    typename OutputSampleType::ConstIterator iterOutput = m_Output->Begin();
    for (unsigned int index = 0 ; index < sampleSize ; index++)
    {
      std::cout << (iterOutput.GetMeasurementVector())[0] << ", ";
      ++iterOutput;
    }
    std::cout << std::endl;
  }
  */
}

} // end of namespace Statistics
} // end of namespace neurolib

#endif
