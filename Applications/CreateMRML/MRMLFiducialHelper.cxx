#include "MRMLFiducialHelper.h"

MRMLFiducialHelper::MRMLFiducialHelper()
{
  m_X = 0.0 ;
  m_Y = 0.0 ;
  m_Z = 0.0 ;
  m_oW = 0.0 ;
  m_oX = 0.0 ;
  m_oY = 0.0 ;
  m_oZ = 0.0 ;
  m_Selected = true ;
  m_Visible = true ;
  m_TextScale = 6.0 ;
  m_RGB.resize( 3 , .5 ) ;
}

void  MRMLFiducialHelper::SetTextScale( double scale )
{
  m_TextScale = scale ;
}

std::string MRMLFiducialHelper::GetType()
{
  return "FiducialList" ;
}

double  MRMLFiducialHelper::GetTextScale()
{
 return m_TextScale ; 
}

void MRMLFiducialHelper::Print()
{
  MRMLColorableHelper::Print() ;
  std::cout << "Id: "<< m_Id << std::endl ;
  std::cout << "Label: " << m_Label << std::endl ;
  std::cout <<  "Position: " << m_X << " " <<  m_Y << " " << m_Z << std::endl ;
  std::cout <<  "Orientation: " << m_oW << " " << m_oX << " " <<  m_oY << " " <<  m_oZ << std::endl ;
  std::cout <<  "Is Selected: " << m_Selected << std::endl ;
  std::cout << "Is Visible: " << m_Visible << std::endl ;
}

void MRMLFiducialHelper::SetId( std::string id)
{
  m_Id = id ;
}

std::string MRMLFiducialHelper::GetId()
{
  return m_Id ;
}

void MRMLFiducialHelper::SetLabelText( std::string label )
{
  m_Label = label ;
}

std::string MRMLFiducialHelper::GetLabel()
{
  return m_Label ;
}


void MRMLFiducialHelper::SetPosition( float x , float y , float z )
{
  m_X = x ;
  m_Y = y ;
  m_Z = z ;
}


std::vector< float > MRMLFiducialHelper::GetPosition()
{
  std::vector< float > xyz( 3 ) ;
  xyz[ 0 ] = m_X ;
  xyz[ 1 ] = m_Y ;
  xyz[ 2 ] = m_Z ;
  return xyz ;
}

void MRMLFiducialHelper::SetOrientation( float w , float x , float y , float z )
{
  m_oW = z ;
  m_oX = x ;
  m_oY = y ;
  m_oZ = z ;
}


std::vector< float > MRMLFiducialHelper::GetOrientation()
{
  std::vector< float > wxyz( 4 ) ;
  wxyz[ 0 ] = m_oW ;
  wxyz[ 1 ] = m_oX ;
  wxyz[ 2 ] = m_oY ;
  wxyz[ 3 ] = m_oZ ;
  return wxyz ;
}

void MRMLFiducialHelper::Selected( bool selected )
{
  m_Selected = selected ;
}

bool MRMLFiducialHelper::IsSelected()
{
  return m_Selected ;
}

void MRMLFiducialHelper::Visibility( bool visible )
{
  m_Visible = visible ;
}

bool MRMLFiducialHelper::IsVisible()
{
  return m_Visible ;
}

int MRMLFiducialHelper::SetSelectedColor( float R , float G , float B )
{
   int error1 , error2 , error3 ;
   CORRECTVALUEMACRO( R , error1 ) ;
   CORRECTVALUEMACRO( G , error2 ) ;
   CORRECTVALUEMACRO( R , error3 ) ;
   if( error1 || error2 || error3 )
   {
      return 1 ;
   }
   m_RGB[ 0 ] = R ;
   m_RGB[ 1 ] = G ;
   m_RGB[ 2 ] = B ;
   return 0 ;
}
  
std::vector< float > MRMLFiducialHelper::GetSelectedColor()
{
   return m_RGB ;
}

