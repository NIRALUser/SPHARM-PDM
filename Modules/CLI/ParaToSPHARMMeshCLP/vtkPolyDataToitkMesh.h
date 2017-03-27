#ifndef __vtkPolyDataToitkMesh_h__
#define __vtkPolyDataToitkMesh_h__

#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkPolyData.h"
#include <vtkVersion.h>
#include "itkDefaultDynamicMeshTraits.h"
#include "itkMesh.h"
#include "itkTriangleCell.h"
#include <vtkSmartPointer.h>

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

  typedef itk::DefaultDynamicMeshTraits<float, 3, 3, double> TriangleMeshTraits;
  typedef itk::Mesh<float, 3, TriangleMeshTraits>            TriangleMeshType;

  /**
  The SetInput method provides pointer to the vtkPolyData
  */
  void SetInput( vtkSmartPointer<vtkPolyData> polydata);

  TriangleMeshType * GetOutput();

  void ConvertvtkToitk();

  TriangleMeshType::Pointer m_itkMesh;

  vtkSmartPointer<vtkPolyData> m_PolyData;

};

#endif
