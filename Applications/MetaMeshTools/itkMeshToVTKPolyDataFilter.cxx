/*=========================================================================

MN, 09/21/2005, marc@bwh.harvard.edu

=========================================================================*/
#ifndef _itkMeshToVTKPolyDataFilter_cxx
#define _itkMeshToVTKPolyDataFilter_cxx

#include "itkMeshToVTKPolyDataFilter.h"

namespace itk
{



/**
 * Constructor
 */
MeshToVTKPolyDataFilter
::MeshToVTKPolyDataFilter()
{

  m_polyData = vtkPolyData::New();

}


/**
 * Destructor
 */
MeshToVTKPolyDataFilter
::~MeshToVTKPolyDataFilter()
{
  if( m_polyData )
    {
    m_polyData->Delete();
    m_polyData = 0;
    }
}



/**
 * Set mesh as input 
 */
void
MeshToVTKPolyDataFilter
::SetInput( MeshToVTKPolyDataFilter::MeshType * inputMesh )
{
  m_mesh = inputMesh;

  // connect the itk to the vtk pipeline
  /*m_polyData->SetUpdateInformationCallback( m_mesh->GetUpdateInformationCallback());
  m_polyData->SetPipelineModifiedCallback( m_mesh->GetPipelineModifiedCallback());
  m_polyData->SetWholeExtentCallback( m_mesh->GetWholeExtentCallback());
  m_polyData->SetSpacingCallback( m_mesh->GetSpacingCallback());
  m_polyData->SetOriginCallback( m_mesh->GetOriginCallback());
  m_polyData->SetScalarTypeCallback( m_mesh->GetScalarTypeCallback());
  m_polyData->SetNumberOfComponentsCallback( m_mesh->GetNumberOfComponentsCallback());
  m_polyData->SetPropagateUpdateExtentCallback( m_mesh->GetPropagateUpdateExtentCallback());
  m_polyData->SetUpdateDataCallback( m_mesh->GetUpdateDataCallback());
  m_polyData->SetDataExtentCallback( m_mesh->GetDataExtentCallback());
  m_polyData->SetBufferPointerCallback( m_mesh->GetBufferPointerCallback());
  m_polyData->SetCallbackUserData( m_mesh->GetCallbackUserData());*/

}

/**
 * Get polyData as output
 */
const vtkPolyData *
MeshToVTKPolyDataFilter
::GetOutput() const
{
  return m_polyData;
}

/**
 * Do the conversion
 */

void
MeshToVTKPolyDataFilter
::Update()
{
  GenerateData();
}

void
MeshToVTKPolyDataFilter
::GenerateData()
{

  //std::cout << "Generating data." << std::endl;
  //std::cout << "May need to clear out the polyData object before putting data in it. So that the function can be called twice without appending data." << std::endl;

  m_polyData->Reset();

  // Get the number of points in the mesh
  int numPoints = m_mesh->GetNumberOfPoints();
  if(numPoints == 0)
  {
    m_mesh->Print(std::cerr);
  }

  // Create the vtkPoints object and set the number of points
  vtkPoints* vpoints = vtkPoints::New();
  vpoints->SetNumberOfPoints(numPoints);

  // iterate over all the points in the itk mesh filling in
  // the vtkPoints object as we go
  MeshType::PointsContainer::Pointer points = m_mesh->GetPoints();
  for(MeshType::PointsContainer::Iterator i = points->Begin(); i != points->End(); ++i)
  {
    // Get the point index from the point container iterator
    int idx = i->Index();
    // Set the vtk point at the index with the the coord array from itk
    // itk returns a const pointer, but vtk is not const correct, so
    // we have to use a const cast to get rid of the const
    //vpoints->SetPoint(idx, const_cast<FLOAT*>(i->Value().GetDataPointer()));
    vpoints->SetPoint(idx, const_cast<vtkFloatingPointType*>(i->Value().GetDataPointer()));
  }
  // Set the points on the vtk grid
  m_polyData->SetPoints(vpoints);

  // Now create the cells using the MulitVisitor
  // 1. Create a MultiVisitor
  MeshType::CellType::MultiVisitor::Pointer mv = 
    MeshType::CellType::MultiVisitor::New();

  // 2. Create a triangle and quadrilateral visitor
  TriangleVisitor::Pointer tv = TriangleVisitor::New();
  //QuadrilateralVisitor::Pointer qv =  QuadrilateralVisitor::New();

  // 3. Set up the visitors
  int vtkCellCount = 0;       // running counter for current cell being inserted into vtk
  int numCells = m_mesh->GetNumberOfCells();
  int *types = new int[numCells];  // type array for vtk
  // create vtk cells and estimate the size
  vtkCellArray* cells = vtkCellArray::New();
  //cells->EstimateSize(numCells, 4);
  cells->EstimateSize(numCells, 3);

  // Set the TypeArray CellCount and CellArray for both visitors
  tv->SetTypeArray(types);
  tv->SetCellCounter(&vtkCellCount);
  tv->SetCellArray(cells);

  //qv->SetTypeArray(types);
  //qv->SetCellCounter(&vtkCellCount);
  //qv->SetCellArray(cells);

  // add the visitors to the multivisitor
  mv->AddVisitor(tv);
  //mv->AddVisitor(qv);
  // Now ask the mesh to accept the multivisitor which
  // will Call Visit for each cell in the mesh that matches the
  // cell types of the visitors added to the MultiVisitor
  m_mesh->Accept(mv);

  int iNumberOfTriangles = cells->GetNumberOfCells();
  //std::cout << "number of triangles = " << iNumberOfTriangles << std::endl;

  // print out the triangles

  /*

  vtkIdType npts;
  vtkIdType *pts;
  for ( int iI=0; iI<iNumberOfTriangles; iI++ ) {
    cells->GetCell( iI, npts, pts);



    {
      std::cout << iI << ": ";
      // output the triangle
      double dPt[3];
      for ( int iJ=0; iJ<3; iJ++ ) {
	vtkIdType ptIndx = pts[iJ];
	vpoints->GetPoint(pts[ptIndx], dPt);

	// now print it out
	std::cout << "( " << dPt[0] << " , " << dPt[1] << " , " << dPt[2] << "); ";

      }

      if ( npts!=3 ) {
	std::cerr << " Npts = " << npts << std::endl;
      } 



      std::cout << std::endl;

    }
  }
  
  */

  /*std::cout << "Number of mesh cells " << numCells << std::endl;
  std::cout << "Number of points " << numPoints << std::endl;
  std::cout << "vtkCellCount = " << vtkCellCount << std::endl;*/



  // Now set the cells on the vtk polydata
  m_polyData->SetPolys(cells);

  // Clean up vtk objects (no vtkSmartPointer ... )
  cells->Delete();
  vpoints->Delete();

}




} // end namespace itk

#endif

