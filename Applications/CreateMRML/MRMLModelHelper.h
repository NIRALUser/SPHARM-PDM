#ifndef _MRMLModelHelper_h
#define _MRMLModelHelper_h

#include "MRMLColorableHelper.h"

class MRMLModelHelper : public MRMLColorableHelper
{
   public:
      std::string GetType()
      {
         return "Model" ;
      }
};

#endif
