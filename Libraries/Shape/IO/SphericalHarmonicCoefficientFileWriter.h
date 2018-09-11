#ifndef __namicSphericalHarmonicCoefficientFileWriter_h
#define __namicSphericalHarmonicCoefficientFileWriter_h

/** \class SphericalHarmonicCoefficientFileReader
 *
 *  \brief This class writes out a list of spherical harmonic coefficients.
 *
 *  \author Christine Xu
 */

#include "SphericalHarmonicSpatialObject.h"

namespace neurolib
{

class SphericalHarmonicCoefficientFileWriterException : public itk::ExceptionObject
{
public:
  /** Run-time information. */
  itkTypeMacro( SphericalHarmonicCoefficientFileWriterException, ExceptionObject );

  /** Constructor. */
  SphericalHarmonicCoefficientFileWriterException(const char *file, unsigned int line,
                                                  const char* message = "Error in IO") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

  /** Constructor. */
  SphericalHarmonicCoefficientFileWriterException(const std::string & file, unsigned int line,
                                                  const char* message = "Error in IO") :
    ExceptionObject(file, line)
  {
    SetDescription(message);
  }

};

class SphericalHarmonicCoefficientFileWriter : public itk::Object
{
public:
  /** SmartPointer typedef support */
  typedef SphericalHarmonicCoefficientFileWriter Self;
  typedef itk::SmartPointer<Self>                Pointer;

  typedef SphericalHarmonicSpatialObject::ScalarType   ScalarType;
  typedef SphericalHarmonicSpatialObject::CoefType     CoefType;
  typedef SphericalHarmonicSpatialObject::CoefListType CoefListType;

  /** Method for creation through the object factory */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  typedef itk::Object Superclass;
  itkTypeMacro(SphericalHarmonicCoefficientFileWriter, Object);

  /** Write out a coef file. */
  void Update(void);

  /** Set the filename  */
  itkSetStringMacro(FileName);

  /** Get the filename */
  itkGetStringMacro(FileName);

  void SetInput(CoefListType& coeflist)
  {
    m_Coefs = coeflist;
  }

protected:
  std::string m_FileName;

  SphericalHarmonicCoefficientFileWriter();
  virtual ~SphericalHarmonicCoefficientFileWriter() ITK_OVERRIDE;
private:
  CoefListType m_Coefs;

};

} // end namespace neurolib

#endif
