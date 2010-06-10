// AttributeMeshSpatialObject.txx: implementation of the AttributeMeshSpatialObject class.
//
//////////////////////////////////////////////////////////////////////

#include "AttributeMeshSpatialObject.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
template < class TMesh, class TAttributeType>
AttributeMeshSpatialObject<TMesh, TAttributeType>::AttributeMeshSpatialObject() 
{
    m_AttributeDimension = 1 ;
    this->SetTypeName("AttributeMeshSpatialObject");
    m_Attributes.clear () ;
    m_Attributes.resize(1);
    m_Attributes[0].clear();
    //cout << "constructor done" << endl ;
}

template < class TMesh, class TAttributeType>
AttributeMeshSpatialObject<TMesh, TAttributeType>::~AttributeMeshSpatialObject()
{

}

template < class TMesh, class TAttributeType>
void AttributeMeshSpatialObject<TMesh, TAttributeType>::ReadAttributes ( std::string attrFilename ) 
{
    int i, j ;
    std::ifstream attrFile ;

    attrFile.open ( attrFilename.c_str() ) ;

    int nPts = this->m_Mesh->GetNumberOfPoints () ;
    
    m_Attributes.clear () ;
    m_Attributes.resize(m_AttributeDimension);
    for (i = 0 ; i < m_AttributeDimension ; i++)
       m_Attributes[i].resize (nPts) ;
    
    //cout << m_AttributeDimension << endl ;
    TAttributeType temp ;
    for ( i = 0 ; i < nPts ; i++ )
    {
        for ( j = 0 ; j < m_AttributeDimension ; j++ )
        {
            attrFile >> temp ;
            m_Attributes[j][i] = temp ;
        }
    }
    attrFile.close () ;
}

template < class TMesh, class TAttributeType>
typename AttributeMeshSpatialObject<TMesh, TAttributeType>::TAttributeListType AttributeMeshSpatialObject<TMesh, TAttributeType>::GetAttribute (unsigned int vertIndex) 
{
  TAttributeListType attr;
  
  if ( vertIndex < m_Attributes[0].size() )
  {
    attr.resize (m_AttributeDimension) ;
    //attr.clear () ;

    for (int i = 0 ; i < m_AttributeDimension ; i++ )
      attr[i] = m_Attributes[i][vertIndex] ;
    //cout << m_AttributeDimension << endl ;
  }
  else
  {
    //cout << "oops" ;
    attr.resize (0) ;
  }
  return attr ;
}

template < class TMesh, class TAttributeType>
void AttributeMeshSpatialObject<TMesh, TAttributeType>::writeToFile(const char *filename) 
{
  std::ofstream outputFile ;
  outputFile.open ( filename ) ;

  int nPts = this->m_Mesh->GetNumberOfPoints () ;
  int i, j ;

  typename TMesh::PointType point ;
    
  outputFile << "ObjectType = Mesh Attribute" << endl ;
  outputFile << "NDims = 3" << endl ;
  outputFile << "AttrDims = " << m_AttributeDimension << endl ;
  outputFile << "ID = 0" << endl ;
  outputFile << "TransformMatrix = 1 0 0 0 1 0 0 0 1" << endl ;
  outputFile << "Offset = 0 0 0" << endl ;
  outputFile << "CenterOfRotation = 0 0 0" << endl ;
  outputFile << "ElementSpacing = 1 1 1" << endl ;
  outputFile << "PointType = MET_FLOAT" << endl ;
  outputFile << "PointDataType = MET_FLOAT" << endl ;
  outputFile << "CellDataType = MET_FLOAT" << endl ;
  outputFile << "NCellTypes = 1" << endl ;
  outputFile << "PointDim = ID x y ..." << endl ;
  outputFile << "NPoints = " << nPts << endl ;
  outputFile << "Points = " << endl ;

  for ( i = 0 ; i < nPts ; i++ )
  {
    if ( this->m_Mesh->GetPoint (i, &point) ) 
    {
      outputFile << i << " " ;
      outputFile << point[0] << " " ;
      outputFile << point[1] << " " ;
      outputFile << point[2] << " " ;
      for ( j = 0 ; j < m_AttributeDimension ; j++ )
      {
          outputFile << m_Attributes[j][i] << " " ;
      }
      outputFile << endl ;
    }
  }

  int nTris = this->m_Mesh->GetNumberOfCells () ;
  typename TMesh::CellAutoPointer cell ;

  outputFile << endl ;
  outputFile << "CellType = TRI" << endl ;
  outputFile << "NCells = " << nTris << endl ;
  outputFile << "Cells = " << endl ;
    
  for ( i = 0 ; i < nTris ; i++ )
  {
    if ( this->m_Mesh->GetCell ( i, cell ) )
    {
      typename TMesh::CellType::PointIdConstIterator points = cell->GetPointIds() ;
      outputFile << i << " " ;
      outputFile << points[0] << " " ;
      outputFile << points[1] << " " ;
      outputFile << points[2] << endl ;
    }
  }

  outputFile.close () ;
}

template < class TMesh, class TAttributeType>
void AttributeMeshSpatialObject<TMesh, TAttributeType>::loadFromFile(const char *filename) 
{
  std::ifstream inputFile ;
  inputFile.open ( filename ) ;

  int i, j ;
  int nPts, nTris ;
  int dummy ;
  char line[200] ;
  bool skip ;

  // skip the beginning 
  skip = true ;
  while ( skip )
  {
    inputFile.getline ( line, 200 ) ;
    skip = strncmp ( line, "AttrDims = ", 11 ) ;
  }
  // get the number of points
  m_AttributeDimension = atoi ( line+11 ) ;
  //cout << m_AttributeDimension << endl ;

  // skip more
  skip = true ;
  while ( skip )
  {
    inputFile.getline ( line, 200 ) ;
    skip = strncmp ( line, "NPoints = ", 10 ) ;
  }
  // get the number of points
  nPts = atoi ( line+10 ) ;
  //cout << nPts << endl ;

  // skip more
  skip = true ;
  while ( skip )
  {
    inputFile.getline ( line, 200 ) ;
    skip = strncmp ( line, "Points = ", 9 ) ;
  }
    
  // Add Points
  TAttributeListType attribute ;
  attribute.clear();
  attribute.resize(m_AttributeDimension);
  m_Attributes.resize (m_AttributeDimension) ;
  for ( i = 0 ; i < m_AttributeDimension ; i++ )
    m_Attributes[i].resize(nPts) ;

  for ( i = 0 ; i < nPts ; i++ )
  {
    typename MeshType::PointType pt;
    inputFile >> dummy >> pt[0] >> pt[1] >> pt[2] ;
    
    for ( j = 0 ; j < m_AttributeDimension ; j++ )
    {
      inputFile >> attribute[j] ;
      m_Attributes[j][i] = attribute[j] ;
    }
    this->m_Mesh->SetPoint(dummy,pt);
    
  }

  // skip the beginning 
  skip = true ;
  while ( skip )
  {
    inputFile.getline ( line, 200 ) ;
    skip = strncmp ( line, "NCells = ", 9 ) ;
  }
  // get the number of points
  nTris = atoi ( line+9 ) ;
    
  // skip more
  skip = true ;
  while ( skip )
  {
    inputFile.getline ( line, 200 ) ;
    skip = strncmp ( line, "Cells = ", 8 ) ;
  }
   
  // Add Cells
  typedef typename MeshType::CellType CellType;
  typedef typename CellType::CellAutoPointer CellAutoPointer;
  this->m_Mesh->SetCellsAllocationMethod ( MeshType::CellsAllocatedDynamicallyCellByCell );
  
  typedef typename MeshType::CellType  CellInterfaceType;
  typedef itk::TriangleCell<CellInterfaceType> TriangleCellType;
    
  int v0, v1, v2 ;

  for ( i = 0 ; i < nTris ; i++ )
  {
    CellAutoPointer cell;
    cell.TakeOwnership ( new TriangleCellType ); 
        
    inputFile >> dummy >> v0 >> v1 >> v2 ;
    cell->SetPointId(0, v0) ;
    cell->SetPointId(1, v1) ;
    cell->SetPointId(2, v2) ;
      
    this->m_Mesh->SetCell(i,cell);
  } 

  inputFile.close () ;
}

template < class TMesh, class TAttributeType>
void AttributeMeshSpatialObject<TMesh, TAttributeType>::SetAttributeDimension (int d)
{
  if ( d > 0 )
     m_AttributeDimension = d ;
}

template < class TMesh, class TAttributeType>
int AttributeMeshSpatialObject<TMesh, TAttributeType>::GetAttributeDimension ()
{
  return m_AttributeDimension ;
}
