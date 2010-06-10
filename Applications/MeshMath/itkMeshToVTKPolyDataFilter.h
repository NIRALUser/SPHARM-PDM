/*=========================================================================

MN, 09/21/05, marc@bwh.harvard.edu

=========================================================================*/
#ifndef __itkMeshToVTKPolyDataFilter_h
#define __itkMeshToVTKPolyDataFilter_h

#include "itkMesh.h"
#include "itkLineCell.h"
#include "itkTriangleCell.h"
#include "vtkPolyData.h"
#include "itkQuadrilateralCell.h"
#include "itkCellInterface.h"
#include "itkDefaultDynamicMeshTraits.h"
#include "itkPolygonCell.h"
#include "itkTetrahedronCell.h"
#include "vtkCellArray.h"

namespace itk
{
  
/** \class MeshToVTKPolyDataFilter
 * \brief Converts VTK ployData into ITK mesh data and plugs a 
 *  vtk data pipeline to an ITK datapipeline.   
 *
 * \ingroup   ImageFilters     
 */
class ITK_EXPORT MeshToVTKPolyDataFilter : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MeshToVTKPolyDataFilter       Self;
  typedef ProcessObject             Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(MeshToVTKPolyDataFilter, ProcessObject);

  /** Some typedefs. */

  const static int PointDimension = 3;
  const static int MaxCellDimension = 3;

  typedef itk::DefaultDynamicMeshTraits< float, PointDimension, MaxCellDimension, double > MeshTraits;
  typedef itk::Mesh< float, PointDimension, MeshTraits > MeshType;


  /*!
  \brief implementation for VTK cell visitors
  \todo documentation
  */
  class VistVTKCellsClass
  {
    vtkCellArray* m_Cells;
    int* m_LastCell;
    int* m_TypeArray;
  public:
    // typedef the itk cells we are interested in
    typedef  itk::CellInterface<  MeshType::PixelType, 
                                          MeshType::CellTraits >  CellInterfaceType;
    typedef itk::LineCell<CellInterfaceType>      floatLineCell;
    typedef itk::PolygonCell<CellInterfaceType>      floatPolygonCell;
    typedef itk::TriangleCell<CellInterfaceType>      floatTriangleCell;
    typedef itk::QuadrilateralCell<CellInterfaceType> floatQuadrilateralCell;
    typedef itk::TetrahedronCell<CellInterfaceType> floatTetrahedronCell;

    /*! Set the vtkCellArray that will be constructed
    */
    void SetCellArray(vtkCellArray* a) 
    {
      m_Cells = a;
    }

    /*! 
    Set the cell counter pointer
    */
    void SetCellCounter(int* i)
    {
      m_LastCell = i;
    }

    /*!
    Set the type array for storing the vtk cell types
    */
    void SetTypeArray(int* i)
    {
      m_TypeArray = i;
    }

    /*!
    Visit a line and create the VTK_LINE cell   
    */
    void Visit(unsigned long , floatLineCell* t)
    {

      vtkIdType tmppts[2];

      typedef floatLineCell::PointIdIterator floatLinePointIdIterator;
      floatLinePointIdIterator pointIditer = t->PointIdsBegin();
      floatLinePointIdIterator pointIdend = t->PointIdsEnd();

      int iI=0;

      while( pointIditer != pointIdend )
	  {
	    tmppts[iI] = *pointIditer;
	    iI++;
	    ++pointIditer;
	  }

      m_Cells->InsertNextCell(2,  tmppts );
      m_TypeArray[*m_LastCell] = VTK_LINE;
      (*m_LastCell)++;
    }

    /*!
    Visit a line and create the VTK_POLYGON cell   
    */
    void Visit(unsigned long , floatPolygonCell* t)
    {
 
      unsigned long num = t->GetNumberOfVertices();
      if (num > 4) {

	vtkIdType *tmppts = new vtkIdType[num];

	typedef floatPolygonCell::PointIdIterator floatPolygonPointIdIterator;
	floatPolygonPointIdIterator pointIditer = t->PointIdsBegin();
	floatPolygonPointIdIterator pointIdend = t->PointIdsEnd();

	int iI=0;

	while( pointIditer != pointIdend )
	  {
	    tmppts[iI] = *pointIditer;
	    iI++;
	    ++pointIditer;
	  }
	

        m_Cells->InsertNextCell(num, tmppts );
        m_TypeArray[*m_LastCell] = VTK_POLYGON;
        (*m_LastCell)++;

	delete [] tmppts;

      }
    }

    /*!
    Visit a triangle and create the VTK_TRIANGLE cell   
    */
    void Visit(unsigned long , floatTriangleCell* t)
    {

      vtkIdType tmppts[3];

      typedef floatTriangleCell::PointIdIterator floatTrianglePointIdIterator;
      floatTrianglePointIdIterator pointIditer = t->PointIdsBegin();
      floatTrianglePointIdIterator pointIdend = t->PointIdsEnd();

      int iI=0;

      while( pointIditer != pointIdend )
	{
	  tmppts[iI] = *pointIditer;
	  iI++;
	  ++pointIditer;
	}


      m_Cells->InsertNextCell(3,  tmppts );

      m_TypeArray[*m_LastCell] = VTK_TRIANGLE;
      (*m_LastCell)++;
    }

    /*! 
    Visit a triangle and create the VTK_QUAD cell 
    */
    void Visit(unsigned long , floatQuadrilateralCell* t)
    {

      vtkIdType tmppts[4];

      typedef floatQuadrilateralCell::PointIdIterator floatQuadPointIdIterator;
      floatQuadPointIdIterator pointIditer = t->PointIdsBegin();
      floatQuadPointIdIterator pointIdend = t->PointIdsEnd();

      int iI=0;

      while( pointIditer != pointIdend )
	{
	  tmppts[iI] = *pointIditer;
	  iI++;
	  ++pointIditer;
	}


      m_Cells->InsertNextCell(4,  tmppts );
      m_TypeArray[*m_LastCell] = VTK_QUAD;
      (*m_LastCell)++;
    }
    
    void Visit(unsigned long , floatTetrahedronCell* t)
    {

      vtkIdType tmppts[4];

      typedef floatTetrahedronCell::PointIdIterator floatTetraPointIdIterator;
      floatTetraPointIdIterator pointIditer = t->PointIdsBegin();
      floatTetraPointIdIterator pointIdend = t->PointIdsEnd();

      int iI=0;

      while( pointIditer != pointIdend )
	{
	  tmppts[iI] = *pointIditer;
	  iI++;
	  ++pointIditer;
	}

      m_Cells->InsertNextCell(4,  tmppts );
      m_TypeArray[*m_LastCell] = VTK_TETRA;
      (*m_LastCell)++;
    }

  };

  typedef itk::CellInterfaceVisitorImplementation<MeshType::PixelType, MeshType::CellTraits,
    itk::QuadrilateralCell< itk::CellInterface<MeshType::PixelType, MeshType::CellTraits > >, 
    VistVTKCellsClass> QuadrilateralVisitor;

  typedef itk::CellInterfaceVisitorImplementation<MeshType::PixelType,
    MeshType::CellTraits,
    itk::TriangleCell<itk::CellInterface<MeshType::PixelType, MeshType::CellTraits > >, 
    VistVTKCellsClass> TriangleVisitor;

  /** Get the output in the form of vtkPolyData **/
  const vtkPolyData*  GetOutput() const;

  /** Set the input in the form of an itk mesh */
  void SetInput( MeshType * );

  /** This call delegate the update to the importer */
  void Update();
  void GenerateData();

protected:
  MeshToVTKPolyDataFilter(); 
  virtual ~MeshToVTKPolyDataFilter(); 

private:
  MeshToVTKPolyDataFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  MeshType::Pointer m_mesh;
  vtkPolyData* m_polyData;

};

} // end namespace itk

#endif



