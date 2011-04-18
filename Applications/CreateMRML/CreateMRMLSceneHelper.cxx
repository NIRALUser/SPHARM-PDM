#include "CreateMRMLSceneHelper.h"

CreateMRMLSceneHelper::CreateMRMLSceneHelper()
{
   m_Scene = vtkMRMLScene::New() ;
}

CreateMRMLSceneHelper::~CreateMRMLSceneHelper()
{
   m_Scene->Delete() ;
}

void CreateMRMLSceneHelper::SetOutputSceneName( std::string name )
{
   m_SceneName = name ;
}

void CreateMRMLSceneHelper::SetInputSceneName( std::string name )
{
  m_InputSceneName = name ;
}

void CreateMRMLSceneHelper::SetInputs( std::vector< MRMLNodeHelper* > arg )
{
   m_Arguments = arg ;
}

int CreateMRMLSceneHelper::Write()
{
  if( CheckDoublons( ) )
  {
    std::cerr << "Error: Some node names appear multiple times. Please verify your input information" << std::endl ;
    PrintArguments() ;
    return 1 ;
  }
  if( !m_InputSceneName.empty() )
  {
    m_Scene->SetURL( m_InputSceneName.c_str() ) ;
    m_Scene->Import() ;
    if( CheckDoublonsWithScene() )
    {
      std::cerr << "Error: Some node names appear inserted are already present in the given scene. Please verify your input information" << std::endl ;
      std::cerr << "Inserted nodes:" << std::endl ;
      PrintArguments() ;
      std::cerr << "Existing nodes in scene:" << std::endl ;
      PrintSceneNodes() ;
      return 1 ;
    }
  }
   //Create scene
   for( unsigned int i = 0 ; i < m_Arguments.size() ; i++ )
{
   if( !m_Arguments[ i ]->GetType().compare( "Transform" ) )
   {
      MRMLTransformHelper *transform= dynamic_cast< MRMLTransformHelper*>( m_Arguments[ i ] ) ;
      if( !transform )
      {
         std::cerr << "Problem dynamic casting MRMLTransformHelper: This should never happen!!!!" << std::endl ;
         return 1 ;
      }
      if( AddTransform( transform ) )
      {
         return 1 ;
      }
   }
   else if( !m_Arguments[ i ]->GetType().compare( "Model" ) )
   {
      MRMLModelHelper *model= dynamic_cast< MRMLModelHelper*>( m_Arguments[ i ] ) ;
      if( !model )
      {
         std::cerr << "Problem dynamic casting MRMLModelHelper: This should never happen!!!!" << std::endl ;
         return 1 ;
      }
      if( AddModel( model ) )
      {
         return 1 ;
      }
   }
   else if( !m_Arguments[ i ]->GetType().compare( "Volume" ) )
   {
      MRMLVolumeHelper *volume = dynamic_cast< MRMLVolumeHelper*>( m_Arguments[ i ] ) ;
      if( !volume )
      {
         std::cerr << "Problem dynamic casting MRMLVolumeHelper: This should never happen!!!!" << std::endl ;
         return 1 ;
      }
      if( AddVolume( volume ) )
      {
         return 1 ;
      }
   }
   else if( !m_Arguments[ i ]->GetType().compare( "FiducialList" ) )
   {
      MRMLFiducialHelper *fiducial = dynamic_cast< MRMLFiducialHelper*>( m_Arguments[ i ] ) ;
      if( !fiducial )
      {
         std::cerr << "Problem dynamic casting MRMLFiducialHelper: This should never happen!!!!" << std::endl ;
         return 1 ;
      }
      if( AddFiducial( fiducial ) )
      {
         return 1 ;
      }
    }
    else
    {
      std::cerr << "An error occured while writing the scene. Unexpected node found. This should never happen" << std::endl ;
      return 1 ;
    }
  }
  m_Scene->SetURL( m_SceneName.c_str() ) ;
  m_Scene->Commit() ;
  return 0 ;
}


void CreateMRMLSceneHelper::PrintArguments()
{
  for( unsigned int i = 0 ; i < m_Arguments.size() ; i++ )
  {
    if( m_Arguments[ i ]->GetNodeName().empty() )
    {
      continue ;
    }
    std::cerr << m_Arguments[ i ]->GetFileName() << " : " << m_Arguments[ i ]->GetNodeName() << std::endl ;
  }
}

void CreateMRMLSceneHelper::PrintSceneNodes()
{
  for( int i = 0 ; i < m_Scene->GetNumberOfNodes() ; i++ )
  {
    std::cerr << m_Scene->GetNthNode( i )->GetName() << std::endl ;
  }
}


//Check that each node name is unique
int CreateMRMLSceneHelper::CheckDoublons()
{
   for( unsigned int i = 0 ; i < m_Arguments.size() - 1 ; i++ )
   {
      for( unsigned int j = i + 1 ; j < m_Arguments.size() ; j++ )
      {
         if( !m_Arguments[ i ]->GetNodeName().empty() )
         {
           if( !m_Arguments[ i ]->GetNodeName().compare( m_Arguments[ j ]->GetNodeName() ) )
           {
              std::cerr << "This node's name appears multiple times: "<< m_Arguments[ i ]->GetNodeName() << std::endl ;
              return 1 ;
           }
         }
      }
   }
   return 0 ;
}

int CreateMRMLSceneHelper::CheckDoublonsWithScene()
{
  for( int i = 0 ; i < m_Scene->GetNumberOfNodes() ; i++ )
  {
    for( unsigned int j = 0 ; j < m_Arguments.size() ; j++ )
    {
      if( !m_Arguments[ j ]->GetNodeName().empty() )
      {
        if( !m_Arguments[ j ]->GetNodeName().compare( m_Scene->GetNthNode( i )->GetName() ) )
        {
          std::cerr << "This node's name appears multiple times: "<< m_Arguments[ i ]->GetNodeName() << std::endl ;
          return 1 ;
        }
      }
    }
  }
  return 0 ;
}


int CreateMRMLSceneHelper::AddFiducial( MRMLFiducialHelper *fiducial )
{
   bool output = 0 ;
   if( fiducial->GetId().empty() )
   {
      return 1 ;
   }
   vtkMRMLFiducialListNode *inode = vtkMRMLFiducialListNode::New() ;
   if( !fiducial->GetNodeName().empty() )
   {      
      inode->SetName( fiducial->GetNodeName().c_str() ) ;
   }
   m_Scene->AddNode( inode ) ;
   std::vector< float > xyz = fiducial->GetPosition() ;
   inode->AddFiducialWithLabelXYZSelectedVisibility( fiducial->GetLabel().c_str() ,
         xyz[ 0 ] ,
         xyz[ 1 ] ,
         xyz[ 2 ] ,
         fiducial->IsSelected() ,
         fiducial->IsVisible()
        ) ;
   std::vector< float > wxyz = fiducial->GetOrientation() ;
   inode->SetNthFiducialOrientation( 0 , wxyz[ 0 ] , wxyz[ 1 ] , wxyz[ 2 ] , wxyz[ 3 ] ) ;
   inode->SetNthFiducialID( 0 , fiducial->GetId().c_str() ) ;
   inode->SetColor( fiducial->GetR() , fiducial->GetG() , fiducial->GetB() ) ;
   inode->SetTextScale( fiducial->GetTextScale() ) ;
   std::vector< float > sc = fiducial->GetSelectedColor() ;
   inode->SetSelectedColor( sc[ 0 ] , sc[ 1 ] , sc[ 2 ] ) ;
   if( SetParentNode( inode , fiducial->GetParentName().c_str() ) )
   {
      output = 1 ;
   }
   inode->Delete() ;
   return output ;
}

int CreateMRMLSceneHelper::AddVolume( MRMLVolumeHelper *volume )
{
   vtkMRMLStorageNode *snode ;
   vtkMRMLScalarVolumeDisplayNode* dnode ;
   vtkMRMLVolumeNode *inode ;
   if( !volume->GetVolumeType().compare( "scalar" ) )
   {
      snode = vtkMRMLVolumeArchetypeStorageNode::New() ;
      dnode = vtkMRMLScalarVolumeDisplayNode::New() ;
      vtkMRMLScalarVolumeNode *vnode = vtkMRMLScalarVolumeNode::New() ;
      if( volume->IsLabelMap() )
      {
         vnode->LabelMapOn() ;
      }
      inode = vnode ;
   }
   else if( !volume->GetVolumeType().compare( "DWI" ) )
   {
      snode = vtkMRMLNRRDStorageNode::New() ;
      dnode = vtkMRMLDiffusionWeightedVolumeDisplayNode::New() ;
      inode = vtkMRMLDiffusionWeightedVolumeNode::New() ;
   }
   else if( !volume->GetVolumeType().compare( "DTI" ) )
   {
      snode = vtkMRMLNRRDStorageNode::New() ;
      dnode = vtkMRMLDiffusionTensorVolumeDisplayNode::New() ;
      inode = vtkMRMLDiffusionTensorVolumeNode::New() ;
   }
   else if( !volume->GetVolumeType().compare( "vector" ) )
   {
      snode = vtkMRMLVolumeArchetypeStorageNode::New() ;
      dnode = vtkMRMLVectorVolumeDisplayNode::New() ;
      inode = vtkMRMLVectorVolumeNode::New() ;
   }
   else
   {
      std::cerr << "Error: Volume type unknown" << std::endl ;
      return EXIT_FAILURE ;
   }
   return AddColorable( volume , dnode , snode , inode ) ;
}

int CreateMRMLSceneHelper::AddModel( MRMLModelHelper *model )
{
   vtkMRMLModelStorageNode *snode = vtkMRMLModelStorageNode::New() ;
   vtkMRMLModelDisplayNode *dnode = vtkMRMLModelDisplayNode::New() ;
   vtkMRMLModelNode *inode = vtkMRMLModelNode::New() ;
   return AddColorable( model , dnode , snode , inode ) ;
}

std::string
CreateMRMLSceneHelper::RemoveExtension( std::string input )
{
  std::string output ;
  if( input.compare( "" ) )
  {
    output = input ;
    //remove folder
    std::size_t pos = output.find_last_of( "/\\" ) ;
    if( pos != std::string::npos )
    {
      output = output.substr( pos + 1 ) ;
    }
    //remove extension
    pos = output.find_last_of( "." ) ;
    if( pos != std::string::npos )
    {
      output = output.substr( 0 , pos ) ;
    }
  }
return output ;
}

int CreateMRMLSceneHelper::AddColorable( MRMLColorableHelper* colorable ,
                  vtkMRMLDisplayNode* dnode ,
                  vtkMRMLStorageNode* snode ,
                  vtkMRMLDisplayableNode* inode
                )
{
   int output = 0 ;
   snode->SetFileName( colorable->GetFileName().c_str() ) ;
   m_Scene->AddNode( snode ) ;
   dnode->SetOpacity( colorable->GetOpacity() ) ;
   double colorRGB[ 3 ] ;
   colorRGB[ 0 ] = colorable->GetR() ;
   colorRGB[ 1 ] = colorable->GetG() ;
   colorRGB[ 2 ] = colorable->GetB() ;
   dnode->SetColor( colorRGB ) ;
   dnode->SetScalarVisibility( 1 ) ;
   if( !strcmp( colorable->GetColorString() , "" ) )
   {
      dnode->SetAndObserveColorNodeID( colorable->GetColor() ) ;
   }
   else
   {
      vtkMRMLColorTableStorageNode *ctsnode = vtkMRMLColorTableStorageNode::New() ;
      ctsnode->SetFileName( colorable->GetColorString() ) ;
      vtkMRMLColorTableNode *ctnode = vtkMRMLColorTableNode::New() ;
      m_Scene->AddNode( ctsnode ) ;
      m_Scene->AddNode( ctnode ) ;
      ctnode->SetName( RemoveExtension( colorable->GetColorString() ).c_str() ) ;
      ctnode->SetAndObserveStorageNodeID( ctsnode->GetID() ) ;
      dnode->SetAndObserveColorNodeID( ctnode->GetID() ) ;
      ctnode->Delete() ;
      ctsnode->Delete() ;
   }
   m_Scene->AddNode( dnode ) ;
   if( !colorable->GetNodeName().empty() )
   {
      inode->SetName( colorable->GetNodeName().c_str() ) ;
   }
   if( !colorable->GetActiveScalarName().empty() )
   {
      dnode->SetActiveScalarName( colorable->GetActiveScalarName().c_str() ) ;
   }
   inode->SetAndObserveDisplayNodeID( dnode->GetID() ) ;
   inode->SetAndObserveStorageNodeID( snode->GetID() ) ;
   if( SetParentNode( inode , colorable->GetParentName().c_str() ) )
   {
      output = 1 ;
   }
   if( !output )
   {
      m_Scene->AddNode( inode ) ;
   }
   inode->Delete() ;
   dnode->Delete() ;
   snode->Delete() ;
   return output ;
}

int CreateMRMLSceneHelper::AddTransform( MRMLTransformHelper *input )
{
   int output = 0 ;
   vtkMRMLTransformStorageNode* snode = vtkMRMLTransformStorageNode::New() ;
   snode->SetFileName( input->GetFileName().c_str() ) ;
   m_Scene->AddNode( snode ) ;
   vtkMRMLLinearTransformNode* lnode = vtkMRMLLinearTransformNode::New() ;
   if( !input->GetNodeName().empty() )
   {
      lnode->SetName( input->GetNodeName().c_str() ) ;
   }
   lnode->SetAndObserveStorageNodeID( snode->GetID() ) ;
   if( !input->GetTransform().empty() )
   {
/*     vtkMatrix4x4* matrix = vtkMatrix4x4::New() ;
     std::vector< double > vec = input->GetTransform() ;
     for( int i = 0 ; i < 4 ; i++ )
     {
       for( int j = 0 ; j < 4 ; j++ )
       {
         matrix->SetElement( i , j , vec[ i*4 + j ] ) ;
       }
     }
     lnode->ApplyTransform( matrix ) ;
     matrix->Delete() ;*/
     //check if transform file exists
     std::vector< std::string > vecstring ;
     std::string path = itksys::SystemTools::GetFilenamePath( m_SceneName ) ;
     itksys::SystemTools::ConvertToUnixSlashes( path ) ;
     path += "/" + input->GetFileName() ;
     //if not, writes transform file
     if( !itksys::SystemTools::FileExists( path.c_str() , true ) )
     {
       std::vector< double > vec = input->GetTransform() ;
       typedef itk::AffineTransform< double , 3 > AffineTransformType ;
       AffineTransformType::Pointer affine = AffineTransformType::New() ;
       itk::Array< double > array( 12 ) ;
       for( unsigned int i = 0 ; i < array.size() ; i++ )
       {
         array[ i ] = vec[ i ] ;
       }
       affine->SetParameters( array ) ;
       typedef itk::TransformFileWriter WriterType ;
       WriterType::Pointer transformWriter = WriterType::New() ;
       transformWriter->SetFileName( path.c_str() ) ;
       transformWriter->AddTransform( affine ) ;
       transformWriter->Update() ;
     }
   }
   if( SetParentNode( lnode ,  input->GetParentName().c_str() ) )
   {
      output = 1 ;
   }
   if( !output )
   {
      m_Scene->AddNode( lnode ) ;
   }
   snode->Delete() ;
   lnode->Delete() ;
   return  output ;
}

int CreateMRMLSceneHelper::SetParentNode( vtkMRMLTransformableNode *child , const char* parentName )
{
   if( !strcmp( parentName , "" ) )
   {
      return 0 ;
   }
   vtkMRMLNode* node = NULL ;
   //     vtkMRMLNode* node = scene->GetNodeByID( input->GetParentName() ) ;//doesn't work for unknown reasons
   //so we replaced it with a manual work around
   for( int i = 0 ; i < m_Scene->GetNumberOfNodes() ; i++ )
   {
      vtkMRMLNode* tempnode = m_Scene->GetNthNode( i) ;
      if( !strcmp( tempnode->GetName() , parentName ) )
      {
         node = tempnode ;
      }
   }
   if( !node )
   {
      std::cerr << "Parent node does not exist. Set parent before child" << std::endl ;
      return 1 ;
   }
   child->SetAndObserveTransformNodeID( node->GetID() ) ;
   return 0 ;
}
