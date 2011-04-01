/*=========================================================================

  Author: Christine Xu

=========================================================================*/

#include "SphericalHarmonicCoefficientFileWriter.h"

#include <stdio.h>

namespace neurolib

{

SphericalHarmonicCoefficientFileWriter::SphericalHarmonicCoefficientFileWriter()
{
  m_FileName = "";
}

SphericalHarmonicCoefficientFileWriter::~SphericalHarmonicCoefficientFileWriter()
{
  
}

void SphericalHarmonicCoefficientFileWriter::Update()
{
  
  if ( m_FileName == "" )
  {
    throw SphericalHarmonicCoefficientFileWriterException(__FILE__, __LINE__, "FileName must be specified");
  }
  
  FILE* file = fopen(m_FileName.c_str(), "w");
  
  if(file == NULL)
  {
    throw SphericalHarmonicCoefficientFileWriterException(__FILE__, __LINE__, "Coef file couldn't be created");
  }
  
  if(m_Coefs.empty())
  {
    fclose(file);
    return;
  }
  
  int count = m_Coefs.size();
  
  fprintf(file, "{ %d,", count);
  
  CoefListType::const_iterator iter = m_Coefs.begin();
  while( iter != m_Coefs.end())
  {
    CoefType elem = *iter;
    fprintf(file, "{%lf, %lf, %lf}", elem[0], elem[1], elem[2]);
    iter++;
    if(iter != m_Coefs.end())
      fprintf(file, ",\n");
  }
  
  fprintf(file, "}");
  
  fclose(file);
}

}//end namespace neurolib
