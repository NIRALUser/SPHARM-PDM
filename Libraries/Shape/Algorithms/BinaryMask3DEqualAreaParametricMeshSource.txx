/*=========================================================================
 Binary Mask to Equal Area Parametrized Surface Mesh
  by Martin Styner
=========================================================================*/
#ifndef _BinaryMask3DEqualAreaParametricMeshSource_txx
#define _BinaryMask3DEqualAreaParametricMeshSource_txx

#include "BinaryMask3DEqualAreaParametricMeshSource.h"
#include "EqualAreaParametricMeshNewtonIterator.h"
#include "itkNumericTraits.h"
#include "TSurfaceNet.h"
#include "TVoxelVolume.h"

#include <itkRegularSphereMeshSource.h>
#include <itkCastImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <math.h>

namespace itk
{

template <class TInputImage>
BinaryMask3DEqualAreaParametricMeshSource<TInputImage>
::BinaryMask3DEqualAreaParametricMeshSource()
{
  // Modify superclass default values, can be overridden by subclasses
  this->SetNumberOfRequiredInputs(1);
  this->SetNumberOfIndexedOutputs(2);
  OutputMeshPointer output
    = dynamic_cast<OutputMeshType *>(this->MakeOutput(0).GetPointer() );

  this->SetNumberOfIndexedOutputs(2);
  this->ProcessObject::SetNthOutput(0, output.GetPointer() );
  this->ProcessObject::SetNthOutput(1, output.GetPointer() );

  m_NumberOfIterations = 500;
  m_ObjectValue = NumericTraits<InputPixelType>::One;
  m_EulerNum = 0;

  m_InitParametricMeshSet = false;
  m_InitParametricMesh = NULL;
}

template <class TInputImage>
BinaryMask3DEqualAreaParametricMeshSource<TInputImage>
::~BinaryMask3DEqualAreaParametricMeshSource()
{
}

template <class TInputImage>
typename BinaryMask3DEqualAreaParametricMeshSource<TInputImage>::OutputMeshType
* BinaryMask3DEqualAreaParametricMeshSource<TInputImage>::GetSurfaceMesh()
  {
  return dynamic_cast<OutputMeshType *>(this->ProcessObject::GetOutput(0) );
  }

template <class TInputImage>
typename BinaryMask3DEqualAreaParametricMeshSource<TInputImage>::OutputMeshType
* BinaryMask3DEqualAreaParametricMeshSource<TInputImage>::GetParametrizationMesh()
  {
  return dynamic_cast<OutputMeshType *>(this->ProcessObject::GetOutput(1) );
  }

template <class TInputImage>
void
BinaryMask3DEqualAreaParametricMeshSource<TInputImage>
::SetInput(const InputImageType* image)
{
  this->ProcessObject::SetNthInput(0,
                                   const_cast<InputImageType *>( image ) );
}

/** Generate the data */
template <class TInputImage>
void
BinaryMask3DEqualAreaParametricMeshSource<TInputImage>
::GenerateData()
{

  if( this->GetNumberOfInputs() < 1 )
    {
    throw BinaryMask3DEqualAreaParametricMeshSourceException(
            __FILE__, __LINE__, "BinaryMask3DEqualAreaParametricMeshSource : Binary image mask not set");
    }

  InputImageConstPointer image =
    static_cast<const InputImageType *>(this->ProcessObject::GetInput(0) );

  // std::cout<<"Extracting Surface net" << std::endl;

  typedef  MinimumMaximumImageCalculator<InputImageType> MinMaxCalcType;
  typename MinMaxCalcType::Pointer minMaxCalc = MinMaxCalcType::New();
  minMaxCalc->SetImage(image);
  minMaxCalc->Compute();
  if( minMaxCalc->GetMinimum() >= minMaxCalc->GetMaximum() )
    {
    throw BinaryMask3DEqualAreaParametricMeshSourceException(__FILE__, __LINE__,
                                                             "BinaryMask3DEqualAreaParametricMeshSource : Image empty");
    }
  // std::cout << "Min/max: " << minMaxCalc->GetMinimum() <<", "
  //        << minMaxCalc->GetMaximum() << std::endl;

  // Casting to unsigned char
  typedef  CastImageFilter<InputImageType, CharImageType> CastFilterType;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(image);
  castFilter->Update();
  typename CharImageType::Pointer castImage = castFilter->GetOutput();

  // std::cout<<"Setting Label " << m_ObjectValue << " to 255, other to 0" << std::endl;
  CharImageIterator iter(castImage, castImage->GetLargestPossibleRegion() );
  while( !iter.IsAtEnd() )
    {
    if( iter.Get() == m_ObjectValue )
      {
      iter.Set(255);
      }
    else
      {
      iter.Set(0);
      }
    ++iter;
    }

  CharImageRegion imageRegion = castImage->GetLargestPossibleRegion();
  int             dim[3];
  dim[0] = imageRegion.GetSize(0);
  dim[1] = imageRegion.GetSize(1);
  dim[2] = imageRegion.GetSize(2);
  CharImageIndex nullIndex;
  nullIndex[0] = 0;
  nullIndex[1] = 0;
  nullIndex[2] = 0;

  const unsigned char *data = (const unsigned char *) &( (*castImage)[nullIndex]);

  TVoxelVolume vox(dim[0], dim[1], dim[2], data);
  TSurfaceNet  surfaceNet(vox, 127);

  IteratorSurfaceNet net;
  // map SurfaceNet to IteratorSurfaceNet
  // allocated vert and face data
  net.nvert = surfaceNet.size();
  net.face = new int[net.nvert * 16]; // upper board, worst case 50% empty
  net.vert = new IteratorSurfaceVertex[net.nvert];
  for( int i = 0; i < net.nvert; i++ )
    {
    TVertex vert = surfaceNet.GetVertex(i);
    (net.vert[i]).x = vert.wh_[0];
    (net.vert[i]).y = vert.wh_[1];
    (net.vert[i]).z = vert.wh_[2];
    (net.vert[i]).count = vert.count_;
    for( int j = 0; j < 14; j++ )
      {
      (net.vert[i]).neighb[j] = vert.neighb_[j];
      }
    }
  int * curface = net.face;
  for( int i = 0; i < net.nvert; i++ )
    {
    for( int j = 0; j < net.vert[i].count - 1; j += 2 )
      {
      int second = (j + 2) % net.vert[i].count;
      if( net.vert[i].neighb[j] > i &&             // 1st corner creates a square
          net.vert[i].neighb[j + 1] > i &&
          net.vert[i].neighb[second] > i )
        {
        *curface++ = i;
        *curface++ = (net.vert[i]).neighb[j];
        *curface++ = (net.vert[i]).neighb[j + 1];
        *curface++ = (net.vert[i]).neighb[second];
        }
      }
    }
  net.nface = (curface - net.face) / 4;

  m_EulerNum = net.nvert - net.nface;
  if( m_EulerNum != 2 )
    {
    char message[1000];
    sprintf(message, "Warning: Euler equation not satisfied. Euler Number : %d - %d = %d",
            net.nvert, net.nface, m_EulerNum);
    throw BinaryMask3DEqualAreaParametricMeshSourceException(__FILE__, __LINE__, message);
    }
  else
    {
    std::cout << "Euler Number ok = " << m_EulerNum << std::endl;
    }

  OutputMeshPointer surfaceMesh = OutputMeshType::New();

  // std::cout<<"Extracting Surface net done " << std::endl;

  // std::cout<<"Computing Initial Parametrization" << std::endl;

  double minErr     = 1e-8;
  double error = 1000000.0;
  char   state[100];

  EqualAreaParametricMeshParameter par;
  // These parameters are exmpirically determined, see
  // C. Brechbühler,  "Constrained Optimization",  Swiss Federal Institute of Technology, Communication Technology
  // Laboratory,
  // Image Science, BIWI-TR-166, 1995 ftp://ftp.vision.ee.ethz.ch/publications   echreports/eth_biwi_00034.ps.gz
  par.max_active = 1000;
  par.print_itn  = -2;
  par.delta      = 3e-7;
  par.constr_tol = 1e-3;
  par.line_tol   = 1e-5;
  par.ineq_low   = 1e-7;
  par.ineq_init  = 1e-2;
  par.ineq_final = 1e-6;
  par.ineq_slack = 2.0;
  par.newton_tol = 1e-4;
  par.rho_init   = 1;
  par.c0rho      = 1;
  par.c1rho      = 0.25;
  par.c2rho      = 0.012;
  par.rho_limit  = 3e-2;
  par.step_small = 0.5;
  par.step_large = 1.0;

  EqualAreaParametricMeshNewtonIterator * optim = NULL;
  double *                                xvec = new double[3 * net.nvert];
  if( m_InitParametricMeshSet )
    {
    // map InitParametricMesh
    PointsContainerPointer points = PointsContainer::New();
    points = m_InitParametricMesh->GetPoints();
    for( int i = 0; i < net.nvert; i++ )
      {
      PointType curPoint = points->GetElement(i);
      xvec[3 * i + 0] = curPoint[0]; xvec[3 * i + 1] = curPoint[1]; xvec[3 * i + 2] = curPoint[2];
      }
    optim = new EqualAreaParametricMeshNewtonIterator(net, par, xvec);
    }
  else
    {
    optim = new EqualAreaParametricMeshNewtonIterator(net, par);
    }
  optim->get_solution(xvec);  // this is the initial parametrization
  // map xvec(initial parametrization) to InitParametricMesh
    {
    m_InitParametricMesh = OutputMeshType::New();
    PointsContainerPointer points = PointsContainer::New();
    for( int i = 0; i < net.nvert; i++ )
      {
      double curVertex[3];
      curVertex[0] = xvec[3 * i + 0]; curVertex[1] = xvec[3 * i + 1]; curVertex[2] = xvec[3 * i + 2];
      points->InsertElement(i, PointType(curVertex) );
      }
    m_InitParametricMesh->SetPoints(points);

    // same connectivity as incoming mesh
    m_InitParametricMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
    for( int i = 0; i < net.nface; i++ )
      {
        { // first triangle
        CellType::CellAutoPointer cellpointer;
        cellpointer.TakeOwnership(new TriangleType);

        unsigned long triPoints[3];
        triPoints[0] = net.face[4 * i + 0];
        triPoints[1] = net.face[4 * i + 1];
        triPoints[2] = net.face[4 * i + 2];
        cellpointer->SetPointIds(triPoints);

        m_InitParametricMesh->SetCell( 2 * i + 0, cellpointer);
        }
        { // second triangle
        CellType::CellAutoPointer cellpointer;
        cellpointer.TakeOwnership(new TriangleType);

        unsigned long triPoints[3];
        triPoints[0] = net.face[4 * i + 2];
        triPoints[1] = net.face[4 * i + 3];
        triPoints[2] = net.face[4 * i + 0];
        cellpointer->SetPointIds(triPoints);

        m_InitParametricMesh->SetCell( 2 * i + 1, cellpointer);
        }
      }
    }

  // std::cout<<"Computing Parametrization Optimization" << std::endl;
  // FINALLY! run the optimization
  unsigned int cnt = 0;
  while( cnt<m_NumberOfIterations && error> minErr )
    {
    cnt++;
    error = optim->iterate();
    }

  if( error <= minErr )
    {
    sprintf(state, "Opt.o.k. it=%d", cnt);
    }
  else
    {
    sprintf(state, "Opt.err. %.2e", error);
    }

  optim->get_solution(xvec);

  // std::cout<<"Computing Parametrization done" << std::endl;

  // Get Results
  OutputMeshPointer paraMesh = OutputMeshType::New();

  // map xvec to ITK Mesh for paraMesh
    {
    PointsContainerPointer points = PointsContainer::New();
    for( int i = 0; i < net.nvert; i++ )
      {
      double curVertex[3];
      curVertex[0] = xvec[3 * i + 0];
      curVertex[1] = xvec[3 * i + 1];
      curVertex[2] = xvec[3 * i + 2];
      // check whether curVertex has Nan's if so raise exceptiont
      if( vnl_math_isnan(curVertex[0]) || vnl_math_isnan(curVertex[1]) || vnl_math_isnan(curVertex[2]) )
        {
        throw BinaryMask3DEqualAreaParametricMeshSourceException(
                __FILE__, __LINE__, "Numerical error in the parameterization, contains NaN");
        }

      points->InsertElement(i, PointType(curVertex) );
      }
    paraMesh->SetPoints(points);

    // same connectivity as incoming mesh
    paraMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
    for( int i = 0; i < net.nface; i++ )
      {
        { // first triangle
        CellType::CellAutoPointer cellpointer;
        cellpointer.TakeOwnership(new TriangleType);

        unsigned long triPoints[3];
        triPoints[0] = net.face[4 * i + 0];
        triPoints[1] = net.face[4 * i + 1];
        triPoints[2] = net.face[4 * i + 2];
        cellpointer->SetPointIds(triPoints);

        paraMesh->SetCell( 2 * i + 0, cellpointer);
        }
        { // second triangle
        CellType::CellAutoPointer cellpointer;
        cellpointer.TakeOwnership(new TriangleType);

        unsigned long triPoints[3];
        triPoints[0] = net.face[4 * i + 2];
        triPoints[1] = net.face[4 * i + 3];
        triPoints[2] = net.face[4 * i + 0];
        cellpointer->SetPointIds(triPoints);

        paraMesh->SetCell( 2 * i + 1, cellpointer);
        }
      }
    }

  // map IteratorSurfaceNet to ITK Mesh for surfaceMesh
    {
    InputSpacingType spacing = image->GetSpacing();

    PointsContainerPointer points = PointsContainer::New();
    for( int i = 0; i < net.nvert; i++ )
      {
      double curVertex[3];
      // scale from index space to mm space

      // TODO: use image transform to get to coordinate
      PointType     phyPoint;
      PointType     origin;
      itk::Index<3> index;

      origin = image->GetOrigin();
      // The mesh needs a flip regardless of orientation in X and Y
      index[0] = (net.vert[i]).x; index[1] = (net.vert[i]).y; index[2] = (net.vert[i]).z;
      image->TransformIndexToPhysicalPoint(index, phyPoint);

      // The mesh needs a translation if origin is not 0 in X and Y
      curVertex[0] = -phyPoint[0];
      curVertex[1] = -phyPoint[1];
      curVertex[2] = phyPoint[2];

      points->InsertElement(i, PointType(curVertex) );
      }

    surfaceMesh->SetPoints(points);

    surfaceMesh->SetCellsAllocationMethod( OutputMeshType::CellsAllocatedDynamicallyCellByCell );
    for( int i = 0; i < net.nface; i++ )
      {
        { // first triangle
        CellType::CellAutoPointer cellpointer;
        cellpointer.TakeOwnership(new TriangleType);

        unsigned long triPoints[3];
        triPoints[0] = net.face[4 * i + 0];
        triPoints[1] = net.face[4 * i + 1];
        triPoints[2] = net.face[4 * i + 2];
        cellpointer->SetPointIds(triPoints);

        surfaceMesh->SetCell( 2 * i + 0, cellpointer);
        }
        { // second triangle
        CellType::CellAutoPointer cellpointer;
        cellpointer.TakeOwnership(new TriangleType);

        unsigned long triPoints[3];
        triPoints[0] = net.face[4 * i + 2];
        triPoints[1] = net.face[4 * i + 3];
        triPoints[2] = net.face[4 * i + 0];
        cellpointer->SetPointIds(triPoints);

        surfaceMesh->SetCell( 2 * i + 1, cellpointer);
        }
      }
    }

  this->ProcessObject::SetNthOutput(0, surfaceMesh.GetPointer() );
  this->ProcessObject::SetNthOutput(1, paraMesh.GetPointer() );

  std::cout << "Computing done" << std::endl;

  this->GetOutput()->SetBufferedRegion( this->GetOutput()->GetRequestedRegion() );

  delete net.face;
  delete net.vert;
}

/** PrintSelf */
template <class TInputImage>
void
BinaryMask3DEqualAreaParametricMeshSource<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os, indent);

  os << indent
     << "ObjectValue: "
     << static_cast<NumericTraits<unsigned char>::PrintType>(m_ObjectValue)
     << std::endl;

  os << indent
     << "NumberOfInterations: "
     << static_cast<NumericTraits<unsigned int>::PrintType>(m_NumberOfIterations)
     << std::endl;

}

}

#endif
