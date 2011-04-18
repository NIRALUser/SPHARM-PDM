#ifndef _CreateMRMLSceneHelper_h
#define _CreateMRMLSceneHelper_h

#include <vtkMRMLScene.h>
#include <vtkMRMLLinearTransformNode.h>
#include <vtkMRMLTransformStorageNode.h>
#include <vtkMRMLScalarVolumeNode.h>
#include <vtkMRMLVolumeArchetypeStorageNode.h>
#include <vtkMRMLScalarVolumeDisplayNode.h>
#include <vtkMRMLDiffusionWeightedVolumeDisplayNode.h>
#include <vtkMRMLVectorVolumeNode.h>
#include <vtkMRMLVectorVolumeDisplayNode.h>
#include <vtkMRMLDiffusionTensorVolumeNode.h>
#include <vtkMRMLDiffusionTensorVolumeDisplayNode.h>
#include <vtkMRMLDiffusionWeightedVolumeNode.h>
#include <vtkMRMLNRRDStorageNode.h>
#include <vtkMRMLModelStorageNode.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLDisplayNode.h>
#include <vtkMRMLColorTableNode.h>
#include <vtkMRMLColorTableStorageNode.h>
#include <vtkMRMLFiducialListStorageNode.h>
#include <vtkMRMLFiducialListNode.h>
#include <itksys/SystemTools.hxx>
#include <itkTransformFileWriter.h>
#include <itkAffineTransform.h>

#include "MRMLModelHelper.h"
#include "MRMLTransformHelper.h"
#include "MRMLVolumeHelper.h"
#include "MRMLFiducialHelper.h"
#include <vector>
#include <string>

class CreateMRMLSceneHelper
{
   public:
      CreateMRMLSceneHelper() ;
      ~CreateMRMLSceneHelper() ;
      void SetOutputSceneName( std::string name ) ;
      void SetInputSceneName( std::string name ) ;
      int Write() ;
      void SetInputs( std::vector< MRMLNodeHelper* > arg ) ;
//      void AddInput( InputClass* arg ) ;
   private:
      std::string m_SceneName ;
      std::string m_InputSceneName ;
      std::vector< MRMLNodeHelper* > m_Arguments ;
      vtkMRMLScene* m_Scene ;
      std::string RemoveExtension( std::string input ) ;
      int CheckDoublons( ) ;
      int CheckDoublonsWithScene( ) ;
      void PrintArguments() ;
      void PrintSceneNodes() ;
      int AddVolume( MRMLVolumeHelper *volume ) ;
      int SetParentNode( vtkMRMLTransformableNode *child , const char* parentName ) ;
      int AddTransform( MRMLTransformHelper *input ) ;
      int AddColorable( MRMLColorableHelper* colorable ,
                        vtkMRMLDisplayNode* dnode ,
                        vtkMRMLStorageNode* snode ,
                        vtkMRMLDisplayableNode* inode
                      ) ;
      int AddModel( MRMLModelHelper *model ) ;
      int AddFiducial( MRMLFiducialHelper *fiducial ) ;
};

#endif
