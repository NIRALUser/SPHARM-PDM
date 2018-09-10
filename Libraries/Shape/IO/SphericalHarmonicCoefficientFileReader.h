#ifndef __namicSphericalHarmonicCoefficientFileReader_h
#define __namicSphericalHarmonicCoefficientFileReader_h

/** \class SphericalHarmonicCoefficientFileReader
 *
 *  \brief This class reads a coefficient file and returns a list of spherical harmonic coefficients.
 *
 *  \author Christine Xu
 */

#include "SphericalHarmonicSpatialObject.h"

#include <string.h>

namespace neurolib
{

class SphericalHarmonicCoefficientFileReaderException : public itk::ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro( ImageFileReaderException, ExceptionObject );

  /** Constructor. */
  SphericalHarmonicCoefficientFileReaderException(const char *file, unsigned int line,
                                                  const char* message = "Error in IO") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  SphericalHarmonicCoefficientFileReaderException(const std::string & file, unsigned int line,
                                                  const char* message = "Error in IO") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

};

class SphericalHarmonicCoefficientFileReader : public itk::Object
{
public:
  /** SmartPointer typedef support */
  typedef SphericalHarmonicCoefficientFileReader Self;
  typedef itk::SmartPointer<Self>                Pointer;

  typedef SphericalHarmonicSpatialObject::ScalarType   ScalarType;
  typedef SphericalHarmonicSpatialObject::CoefType     CoefType;
  typedef SphericalHarmonicSpatialObject::CoefListType CoefListType;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  typedef itk::Object Superclass;
  itkTypeMacro(SphericalHarmonicCoefficientFileReader, Object);

  /** Load a coef file. */
  void Update(void);

  /** Set the filename  */
  itkSetStringMacro(FileName);

  /** Get the filename */
  itkGetStringMacro(FileName);

  void GetOutput(CoefListType& coeflist)
  {
    coeflist = m_Coefs;
  }

protected:
  std::string m_FileName;

  SphericalHarmonicCoefficientFileReader();
  virtual ~SphericalHarmonicCoefficientFileReader() ITK_OVERRIDE;
private:
  CoefListType m_Coefs;

};

} // end namespace neurolib

#endif
