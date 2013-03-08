/*=========================================================================

  Author: Christine Xu

=========================================================================*/

#include "SphericalHarmonicCoefficientFileReader.h"

#include <stdio.h>

namespace neurolib
{

SphericalHarmonicCoefficientFileReader::SphericalHarmonicCoefficientFileReader()
{
  m_FileName = "";
}

SphericalHarmonicCoefficientFileReader::~SphericalHarmonicCoefficientFileReader()
{

}

void SphericalHarmonicCoefficientFileReader::Update()
{
  if( !m_Coefs.empty() )
    {
    m_Coefs.clear();
    }

  int count;

  if( m_FileName == "" )
    {
    throw SphericalHarmonicCoefficientFileReaderException(__FILE__, __LINE__, "FileName must be specified");
    }

  FILE* file = fopen(m_FileName.c_str(), "r");

  if( file == NULL )
    {
    throw SphericalHarmonicCoefficientFileReaderException(__FILE__, __LINE__, "Coef file couldn't be open");
    }

  const char* format_3_real = " { %lf , %lf , %lf } , ";

  /** Test for format exception. */
  if( fscanf(file, " { %d,", &count) != 1 )
    {
    throw SphericalHarmonicCoefficientFileReaderException(__FILE__, __LINE__, "Error in the header of the Coef file");
    }
  for( int i = 0; i < count; i++ )
    {

    ScalarType elem[3];

    if( fscanf(file, format_3_real, &(elem[0]), &(elem[1]), &(elem[2]) ) != 3 )
      {
      throw SphericalHarmonicCoefficientFileReaderException(__FILE__, __LINE__, "Error in the coefficient format");
      }

    CoefType coef = elem;
    m_Coefs.push_back(coef);
    }

  if( fscanf(file, "}") != 0 )
    {
    throw SphericalHarmonicCoefficientFileReaderException(__FILE__, __LINE__,
                                                          "Error in the end of the Coef file: '}' expected");
    }

  fclose(file);
}

} // end namespace neurolib
