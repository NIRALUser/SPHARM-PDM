#ifndef _MRMLTransformHelper_h
#define _MRMLTransformHelper_h

#include "MRMLNodeHelper.h"
#include <vector>

class MRMLTransformHelper : public MRMLNodeHelper
{
   public:
      std::string GetType()
      {
         return "Transform" ;
      }
      std::vector< double > GetTransform()
      {
         return m_Transform ;
      }
      int SetTransform( std::vector< double > transform )
      {
        if( transform.size() != 12 )
        {
          return 1 ;
        }
        m_Transform.clear() ;
        m_Transform = transform ;
        return 0 ;
      }
      void ClearTransform()
      {
        m_Transform.clear() ;
      }
   private:
      std::vector< double > m_Transform ;
};

#endif
