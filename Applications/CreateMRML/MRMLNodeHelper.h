#ifndef _MRMLNodeHelper_h
#define _MRMLNodeHelper_h

#include <string>

#define CORRECTVALUEMACRO( x , error )\
if( x >= 0.0 && x <= 1.0 )\
{\
  error = 0 ;\
}\
else\
{\
  error = 1 ;\
}\

class MRMLNodeHelper
{
   public:
      std::string GetFileName()
      {
         return m_FileName ;
      }
      void SetFileName( std::string fileName )
      {
         m_FileName = fileName ;
      }
      std::string GetParentName()
      {
         return m_ParentName ;
      }
      void SetParentName( std::string parentName )
      {
         m_ParentName = parentName ;
      }
      void SetNodeName( std::string nodeName )
      {
         m_NodeName = nodeName ;
      }
      std::string GetNodeName()
      {
         return m_NodeName ;
      }
      virtual std::string GetType() = 0 ;
      virtual void Print()
      {
         std::cout << "Type: " << this->GetType() << std::endl ;
         std::cout << "FileName: " << m_FileName << std::endl ;
         std::cout << "ParentName: " << m_ParentName << std::endl ;
         std::cout << "NodeName: " << m_NodeName << std::endl ;
      }
      virtual ~MRMLNodeHelper() {}
   protected:
      std::string m_FileName ;
      std::string m_ParentName ;
      std::string m_NodeName ;
};

#endif
