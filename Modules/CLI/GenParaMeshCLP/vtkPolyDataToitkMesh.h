#ifndef __vtkPolyDataToitkMesh_h__
#define __vtkPolyDataToitkMesh_h__

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"


/** 
  \class vtkPolyDataToitkMesh
  \brief 
    \warning
  \sa 
  */

class vtkPolyDataToitkMesh
{

 public:

  vtkPolyDataToitkMesh( void );
  virtual ~vtkPolyDataToitkMesh( void );

  typedef itk::DefaultDynamicMeshTraits<double, 3, 3,double,double> TriangleMeshTraits;
  typedef itk::Mesh<double,3, TriangleMeshTraits> TriangleMeshType;

  /**
  The SetInput method provides pointer to the vtkPolyData
  */
  void SetInput( vtkPolyData * polydata);
  TriangleMeshType * GetOutput();
  void ConvertvtkToitk();

  TriangleMeshType::Pointer  m_itkMesh;

  vtkPolyData           * m_PolyData;

  
};

#endif
