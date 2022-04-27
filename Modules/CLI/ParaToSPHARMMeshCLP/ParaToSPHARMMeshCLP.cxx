/*
 * compute the SPHARM coefficients and the associated Mesh
 *
 * author:  Martin Styner
 *
 */

#ifdef _MSC_VER
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4503 )
#pragma warning ( disable: 4284 )
#endif

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>

#include <itkMeshSpatialObject.h>
#include <itkMesh.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkTriangleCell.h>
#include <itkDefaultDynamicMeshTraits.h>
#include <itkProcessObject.h>

#include "ParametricMeshToSPHARMSpatialObjectFilter.h"
#include "SphericalHarmonicSpatialObject.h"
#include "SphericalHarmonicMeshSource.h"
#include "SphericalHarmonicMedialAxisMeshSource.h"
#include "SphericalHarmonicCoefficientFileWriter.h"
#include "SphericalHarmonicCoefficientFileReader.h"
#include "itkMesh3DProcrustesAlignFilter.h"

#include "vtkPolyDataToitkMesh.h"
#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkDoubleArray.h"
#include "vtkIterativeClosestPointTransform.h"
#include "vtkLandmarkTransform.h"
#include "vtkPointLocator.h"
#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <itkAffineTransform.h>
#include <itkTransformFileWriter.h>

#include "argio.h"

#include "ParaToSPHARMMeshCLPCLP.h"

using namespace std;
vtkSmartPointer<vtkPolyData> ReadPolyData(std::string filePath);

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

#define dimension 3

  typedef itk::DefaultDynamicMeshTraits<float, dimension, dimension, double> MeshTrait;
  typedef itk::Mesh<float, dimension, MeshTrait>                             MeshType;

  /** Hold on to the type information specified by the template parameters. */
  typedef  MeshType::Pointer    MeshPointer;



  /** Some convenient typedefs. */
  typedef  MeshType::Pointer                MeshPointer;
  typedef  MeshType::PointsContainerPointer PointsContainerPointer;
  typedef  MeshType::CellsContainerPointer  CellsContainerPointer;
  typedef  MeshType::PointType              PointType;
  typedef itk::MeshSpatialObject<MeshType>                     MeshSpatialObjectType;
  typedef neurolib::ParametricMeshToSPHARMSpatialObjectFilter  SPHARMFilterType;
  typedef itk::Mesh3DProcrustesAlignFilter<MeshType, MeshType> ProcrustesFilterType;
  typedef itk::AffineTransform<double, dimension>              TransformType;

  std::string  outFileName("dummy");
  const char * base_string;
  base_string = outbase.c_str();

  vtkSmartPointer<vtkPolyData>      convSurfaceMesh = vtkSmartPointer<vtkPolyData>::New();

  MeshType::Pointer surfaceMesh = MeshType::New();
  MeshType::Pointer paraMesh = MeshType::New();
  MeshPointer       regParaTemplateMesh;
  MeshPointer       newParaMesh = MeshType::New();
//   typedef itk::SpatialObjectReader<3,float,MeshTrait> ReaderType;
//   ReaderType::Pointer readerSH = ReaderType::New();
  neurolib::SphericalHarmonicSpatialObject::Pointer flipTemplateSO = neurolib::SphericalHarmonicSpatialObject::New();

  try
    {
    convSurfaceMesh = ReadPolyData( inSurfFile.c_str() );

    vtkPolyDataToitkMesh VTKITKConverter;
    VTKITKConverter.SetInput( convSurfaceMesh );
    surfaceMesh = VTKITKConverter.GetOutput();
//     cout<<"test3"<<endl;
// delete (VTKITKConverter);
//     cout<<"test9"<<endl;
    }
  catch( itk::ExceptionObject &ex )
    {
    std::cout << ex.GetDescription() << std::endl;
    return 1;
    }

  try
    {
    vtkSmartPointer<vtkPolyData> convParaMesh = vtkSmartPointer<vtkPolyData>::New();
    convParaMesh = ReadPolyData( inParaFile.c_str() );

    vtkPolyDataToitkMesh VTKITKConverter2;
    VTKITKConverter2.SetInput( convParaMesh );
    paraMesh = VTKITKConverter2.GetOutput();

    // delete (VTKITKConverter2);
    }
  catch( itk::ExceptionObject &ex )
    {
    std::cout << ex.GetDescription() << std::endl;
    return 1;
    }

  if( flipTemplateFileOn )
    {
    if( debug )
      {
      std::cout << flipTemplateFile.c_str() << endl;
      }

    try
      {
      typedef neurolib::SphericalHarmonicCoefficientFileReader CoefReaderType;
      CoefReaderType::Pointer                                coefreader = CoefReaderType::New();
      neurolib::SphericalHarmonicSpatialObject::CoefListType coeflist;

      coefreader->SetFileName(flipTemplateFile.c_str() );
      coefreader->Update();
      coefreader->GetOutput(coeflist);

      flipTemplateSO->SetCoefs(coeflist);
      }
    catch( itk::ExceptionObject &ex )
      {
      std::cout << ex.GetDescription() << std::endl;
      return 1;
      }
    }
  if( debug )
    {
    std::cout << "Preprocessing done" << std::endl;
    }
  try
    {
    if( regParaTemplateFileOn )
      {
      try
        {
        if( debug )
          {
          std::cout << "Reading in regTemplateMesh" << std::endl;
          }

        vtkSmartPointer<vtkPolyData> convRegParaTemplateMesh = vtkSmartPointer<vtkPolyData>::New();
        convRegParaTemplateMesh = ReadPolyData( regParaTemplateFile.c_str() );

        vtkPolyDataToitkMesh VTKITKConverter4;
        VTKITKConverter4.SetInput(convRegParaTemplateMesh);
        regParaTemplateMesh = VTKITKConverter4.GetOutput();
        regParaTemplateMesh->Update();

        // delete (VTKITKConverter4);

        // 2. use ICP VTK filter -> T
        vtkSmartPointer<vtkIterativeClosestPointTransform> ICPTransform = vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
        vtkSmartPointer<vtkLandmarkTransform>              LT = ICPTransform->GetLandmarkTransform();
        LT->SetModeToRigidBody();
        ICPTransform->SetSource(convRegParaTemplateMesh);
        ICPTransform->SetTarget(convSurfaceMesh);
        ICPTransform->SetMaximumNumberOfIterations(50);

        ICPTransform->Update();

        // 3. get points on surface closest to T*selected points from regTemplate
        double transPoints[3][3];
        if( regParaPointsOn )
          {
          if( debug )
            {
            std::cout << "Reading in regParaPoints" << std::endl;
            }
          std::string   line;
          long int      index;
          char          colon;
          double        points[3][3];
          std::ifstream ptReader;
          ptReader.open(regParaPoints.c_str(), std::ios::in);
          int    pointnum;
          int    i;
          double x, y, z;
          if( ptReader.fail() )
            {
            std::cout << "Error reading regParaPoints file" << std::endl;
            return 1;
            }
          if( ptReader.is_open() )
            {
            std::getline(ptReader, line);
            pointnum = 0;
            while( !ptReader.eof() && pointnum < 3 )
              {
              ptReader >> index >> colon >> x >> y >> z;
              points[pointnum][0] = x;
              points[pointnum][1] = y;
              points[pointnum][2] = z;
              pointnum++;
              }

            if( pointnum < 3 )
              {
              std::cout << "Not enough points" << std::endl;
              return 1;
              }
            ptReader.close();
            }
          // Transform points
          for( i = 0; i < 3; i++ )
            {
            ICPTransform->TransformPoint(points[i], transPoints[i]);
            }
          }
        else
          {
          std::cout << "Need regParaPoints along with regParaTemplate" << std::endl;
          return 1;
          }
        // 4. compute and apply rotation of parameter
        if( debug )
          {
          std::cout << "Calculating parameter rotation" << std::endl;
          }
        vtkSmartPointer<vtkPointLocator> ptLoc = vtkSmartPointer<vtkPointLocator>::New();
        ptLoc->SetDataSet(convSurfaceMesh);
        // Find the id of the transformed points
        long int id[3];
        id[0] = ptLoc->FindClosestPoint(transPoints[0]);
        id[1] = ptLoc->FindClosestPoint(transPoints[1]);
        id[2] = ptLoc->FindClosestPoint(transPoints[2]);

        // Get the parametrization for the points
        PointsContainerPointer paraPts = paraMesh->GetPoints();
        PointType              pts[3] = {paraPts->GetElement(id[0]), paraPts->GetElement(id[1]), paraPts->GetElement(
                                           id[2])};

        // Calculate the theta and phi for the points
        double thetaPts[3], phiPts[3];
        for( int i = 0; i < 3; i++ )
          {
          thetaPts[i] = acos(pts[i][2]);           // 0 .. M_PI
          phiPts[i] = atan2(pts[i][1], pts[i][0]); // -M_PI ... M_PI
          if( phiPts[i] < 0 )
            {
            phiPts[i] = phiPts[i] + 2 * M_PI;                 // phiPts: 0...2*M_PI
            }
          }

        /*Compute the rotation so that
        a) The north pole is aligned to pts[0]
        (theta,phi) = (thetaPts[0],phiPts[0]) maps to (0,0)
        Rotation matrix given by Ry(-thetaPts[0])*Rz(-phiPts[0])
        where Rz(phi) and Ry(theta) are rotations about the z and y
        axes by angles phi and theta respectively.
        Foley, vanDam, et. al. 2nd ed. in C, page 215.

        b) phi for pts[1] after alignment is 0 or Pi
        c) The longitude for pts[2] is between 0 and Pi.

        */
        TransformType::Pointer    rotMatPhi = TransformType::New();
        TransformType::Pointer    rotMatTheta = TransformType::New();
        TransformType::Pointer    rotMat = TransformType::New();
        TransformType::MatrixType matrixPhi, matrixTheta, matrix;

        double phiRot = -phiPts[0];

        matrixPhi[0][0] = cos(phiRot);
        matrixPhi[0][1] = -sin(phiRot);
        matrixPhi[0][2] = 0;
        matrixPhi[1][0] = sin(phiRot);
        matrixPhi[1][1] = cos(phiRot);
        matrixPhi[1][2] = 0;
        matrixPhi[2][0] = 0;
        matrixPhi[2][1] = 0;
        matrixPhi[2][2] = 1;

        rotMatPhi->SetMatrix(matrixPhi);

        double thetaRot = -thetaPts[0];

        matrixTheta[0][0] = cos(thetaRot);
        matrixTheta[0][1] = 0;
        matrixTheta[0][2] = sin(thetaRot);
        matrixTheta[1][0] = 0;
        matrixTheta[1][1] = 1;
        matrixTheta[1][2] = 0;
        matrixTheta[2][0] = -sin(thetaRot);
        matrixTheta[2][1] = 0;
        matrixTheta[2][2] = cos(thetaRot);

        rotMatTheta->SetMatrix(matrixTheta);

        rotMatPhi->Compose(rotMatTheta);
        rotMat = rotMatPhi;
        matrix = rotMat->GetMatrix();

        // Calculate where pts map to
        PointType rotPts[3];
        rotPts[0] =  rotMat->TransformPoint(pts[0]);
        rotPts[1] =  rotMat->TransformPoint(pts[1]);
        rotPts[2] =  rotMat->TransformPoint(pts[2]);

        double rotPhi[3];
        for( int i = 0; i < 3; i++ )
          {
          if( rotPts[i][2] > 1 )
            {
            rotPts[i][2] = 1;
            }
          if( rotPts[i][2] < -1 )
            {
            rotPts[i][2] = -1;
            }
          rotPhi[i] = atan2(rotPts[i][1], rotPts[i][0]);
          if( rotPhi[i] < 0 )
            {
            rotPhi[i] = rotPhi[i] + 2 * M_PI;
            }
          }

        int    phiFlag = 0;
        double testPhi3 = rotPhi[2] - rotPhi[1];
        if( testPhi3 < 0 )
          {
          testPhi3 = testPhi3 + 2 * M_PI;
          }
        if( testPhi3 > M_PI )
          {
          phiFlag = 1;
          }

        PointsContainerPointer paraPoints = paraMesh->GetPoints();
        int                    nvert = paraPoints->Size();
        for( int i = 0; i < nvert; i++ )
          {
          PointType curPoint =  paraPoints->GetElement(i);
          PointType rotPoint = rotMat->TransformPoint(curPoint);
          if( rotPoint[2] > 1 )
            {
            rotPoint[2] = 1;
            }
          if( rotPoint[2] < -1 )
            {
            rotPoint[2] = -1;
            }
          double theta = acos(rotPoint[2]);
          double phi = atan2(rotPoint[1], rotPoint[0]);
          if( phi < 0 )
            {
            phi = phi + 2 * M_PI;
            }
          double phiNew = phi - rotPhi[1] + (phiFlag * M_PI);
          if( phiNew < 0 )
            {
            phiNew =  phiNew + 2 * M_PI;
            }
          if( phiNew >= 2 * M_PI )
            {
            phiNew = phiNew - 2 * M_PI;
            }

          PointType newPoint;
          newPoint[0] = sin(theta) * cos(phiNew);
          newPoint[1] = sin(theta) * sin(phiNew);
          newPoint[2] = rotPoint[2];
          paraPoints->SetElement(i, newPoint);
          }

        PointType testPts[3] = {paraPts->GetElement(id[0]), paraPts->GetElement(id[1]), paraPts->GetElement(id[2])};
        double    testTheta[3], testPhi[3];

        if( debug )
          {
          std::cout << "Final param.s of template points: " << std::endl;
          }
        for( int i = 0; i < 3; i++ )
          {
          testTheta[i] = acos(testPts[i][2]);  // 0 .. M_PI
          testPhi[i] = atan2(testPts[i][1], testPts[i][0]);
          if( testPhi[i] < 0 )
            {
            testPhi[i] = testPhi[i] + 2 * M_PI;
            }
          if( debug )
            {
            std::cout << "ThetaPhi[" << i << "] = (" << testTheta[i] << "," << testPhi[i] << "); ";
            }
          }
        if( debug )
          {
          std::cout << std::endl;
          }

        // 5. use that parametrization
        CellsContainerPointer paraCells = paraMesh->GetCells();
        newParaMesh->SetPoints(paraPoints);
        newParaMesh->SetCells(paraCells);
        MeshSpatialObjectType::Pointer meshSO = MeshSpatialObjectType::New();
        meshSO->SetMesh(newParaMesh);

        // if( writeFlipRegPara)
        //         {
        //           if(debug) std::cout << "saving regPara parametrization" << std::endl;
        //           MeshWriterType::Pointer metaWriter =  MeshWriterType::New();
        //           metaWriter->SetInput(meshSO);
        //           outFileName.erase();
        //           outFileName.append(base_string);
        //           outFileName.append("_regPara_para.meta");
        //           metaWriter->SetFileName(outFileName.c_str());
        //           metaWriter->Update();
        //         }

        // 7. in if regTemplate, add section for if regParaTemplate and use S2 for procrustes alignment
        }
      catch( itk::ExceptionObject &ex )
        {
        std::cout << ex.GetDescription() << "while doing regParaTemplateFile" << std::endl;
        return 1;
        }
      }

    SPHARMFilterType::Pointer spharmFilter = SPHARMFilterType::New();
    spharmFilter->SetInputSurfaceMesh(surfaceMesh);

    spharmFilter->SetInputParametrizationMesh(paraMesh);

    if( regParaTemplateFileOn )
      {
      spharmFilter->SetInputParametrizationMesh(newParaMesh);
      }
    else
      {
      spharmFilter->SetInputParametrizationMesh(paraMesh);
      }

    spharmFilter->SetDegree(spharmDegree);
    if( NoParaAlignFlag || regParaTemplateFileOn )
      {
      spharmFilter->SetParaEllipseAlignment(false);
      }
    else
      {
      spharmFilter->SetParaEllipseAlignment(true);
      }
    if( flipTemplateFileOn && !regParaTemplateFileOn )
      {
      spharmFilter->SetFlipTemplate(flipTemplateSO);
      }
    if( finalFlipIndex > 0 && finalFlipIndex <= 7 )
      {
      spharmFilter->SetParametrizationFlipIndex(finalFlipIndex);
      }
    spharmFilter->GenerateData();
    // spharmFilter->Update();

    if( debug )
      {
      std::cout << "saving par aligned data" << std::endl;
      }

    neurolib::SphericalHarmonicSpatialObject::Pointer      spharmSO = spharmFilter->GetOutput();
    neurolib::SphericalHarmonicSpatialObject::CoefListType coeflist;
    spharmSO->GetCoefs(coeflist);

    // save associated surface
    neurolib::SphericalHarmonicMeshSource::Pointer meshsrc = neurolib::SphericalHarmonicMeshSource::New();
    meshsrc->SetCoefs(coeflist);
    meshsrc->SetLevel(subdivLevel);
    meshsrc->Update();

    MeshType* meshSH;
    meshSH = meshsrc->GetOutput();
    vtkSmartPointer<vtkPolyDataWriter> vtkwriter = vtkSmartPointer<vtkPolyDataWriter>::New();

	 if(debug)
		 std::cout<<"saving medial axis mesh source"<<std::endl;

    // save associated surface
	 neurolib::SphericalHarmonicMedialAxisMeshSource::Pointer medialmeshsrc = neurolib::SphericalHarmonicMedialAxisMeshSource::New();
	 medialmeshsrc->SetCoefs(coeflist);
	 medialmeshsrc->SetThetaPhiIteration(thetaIteration,phiIteration);
	 medialmeshsrc->Update();

	 if(medialMesh)
	 {
		MeshType* medialmeshSH;
		medialmeshSH = medialmeshsrc->GetOutput();
		
        itkMeshTovtkPolyData ITKVTKConverter5;
        ITKVTKConverter5.SetInput( medialmeshSH );
		
        vtkSmartPointer<vtkPolyData> polydataAtt = ITKVTKConverter5.GetOutput();
		// START Add scalars for visualization
		double* Radius;
		Radius=medialmeshsrc->GetRadius();

        vtkSmartPointer<vtkFloatArray> scalars_radius = vtkSmartPointer<vtkFloatArray>::New();
		scalars_radius->SetNumberOfComponents(1);
		scalars_radius->SetName("radius");

		int medialMesh_npoints = thetaIteration * phiIteration;
		for(int i=0; i<thetaIteration; i++)
		  for (int j=0; j<phiIteration; j++)
		     scalars_radius->InsertNextValue(Radius[i]);

		polydataAtt->GetPointData()->AddArray(scalars_radius);   

		std::ofstream outfile_r;
		outFileName.erase();
		outFileName.append(base_string);
		outFileName.append("_medialMeshRadius.txt");
		outfile_r.open(outFileName.c_str(), ios::out);
		// print the header
		outfile_r << "NUMBER_OF_POINTS=" << medialMesh_npoints << std::endl;
		outfile_r << "DIMENSION=" << 1 << std::endl;
		outfile_r << "TYPE=Scalar" << std::endl;
		for(int i=0; i<thetaIteration; i++)
		  for (int j=0; j<phiIteration; j++)
		    outfile_r<<Radius[i] << std::endl;
		outfile_r.close();

		double* Area;
		Area=medialmeshsrc->GetArea();

        vtkSmartPointer<vtkFloatArray> scalars_area = vtkSmartPointer<vtkFloatArray>::New();
		scalars_area->SetNumberOfComponents(1);
		scalars_area->SetName("area");

		for(int i=0; i<thetaIteration; i++)
		  for (int j=0; j<phiIteration; j++)
		     scalars_area->InsertNextValue(Area[i]);

		polydataAtt->GetPointData()->AddArray(scalars_area);   

		std::ofstream outfile_a;
		outFileName.erase();
		outFileName.append(base_string);
		outFileName.append("_medialMeshArea.txt");
		outfile_a.open(outFileName.c_str(), ios::out);
		// print the header
		outfile_a << "NUMBER_OF_POINTS=" << medialMesh_npoints << std::endl;
		outfile_a << "DIMENSION=" << 1 << std::endl;
		outfile_a << "TYPE=Scalar" << std::endl;
		for(int i=0; i<thetaIteration; i++)
		  for (int j=0; j<phiIteration; j++)
		    outfile_a<<Area[i] << std::endl;
		outfile_a.close();

		double* partialArea;
		partialArea=medialmeshsrc->GetPartialArea();

        vtkSmartPointer<vtkFloatArray> scalars_partialArea = vtkSmartPointer<vtkFloatArray>::New();
		scalars_partialArea->SetNumberOfComponents(1);
		scalars_partialArea->SetName("partial_area");

		for(int i=0; i<medialMesh_npoints; i++)
		     scalars_partialArea->InsertNextValue(partialArea[i]);

		polydataAtt->GetPointData()->AddArray(scalars_partialArea);   

		std::ofstream outfile_pa;
		outFileName.erase();
		outFileName.append(base_string);
		outFileName.append("_medialMeshPartialArea.txt");
		outfile_pa.open(outFileName.c_str(), ios::out);
		// print the header
		outfile_pa << "NUMBER_OF_POINTS=" << medialMesh_npoints << std::endl;
		outfile_pa << "DIMENSION=" << 1 << std::endl;
		outfile_pa << "TYPE=Scalar" << std::endl;
		for(int i=0; i<medialMesh_npoints; i++)
		    outfile_pa<<partialArea[i] << std::endl;
		outfile_pa.close();

		double* partialRadius;
		partialRadius=medialmeshsrc->GetPartialRadius();

        vtkSmartPointer<vtkFloatArray> scalars_partialRadius = vtkSmartPointer<vtkFloatArray>::New();
		scalars_partialRadius->SetNumberOfComponents(1);
		scalars_partialRadius->SetName("partial_radius");

		for(int i=0; i<medialMesh_npoints; i++)
		     scalars_partialRadius->InsertNextValue(partialRadius[i]);

		polydataAtt->GetPointData()->AddArray(scalars_partialRadius);   

		std::ofstream outfile_pr;
		outFileName.erase();
		outFileName.append(base_string);
		outFileName.append("_medialMeshPartialRadius.txt");
		outfile_pr.open(outFileName.c_str(), ios::out);
		// print the header
		outfile_pr << "NUMBER_OF_POINTS=" << medialMesh_npoints << std::endl;
		outfile_pr << "DIMENSION=" << 1 << std::endl;
		outfile_pr << "TYPE=Scalar" << std::endl;
		for(int i=0; i<medialMesh_npoints; i++)
		    outfile_pr<<partialRadius[i] << std::endl;
		outfile_pr.close();

		// END Add scalars for visualization

		vtkwriter->SetInputData(polydataAtt);
		outFileName.erase();
		outFileName.append(base_string);
        outFileName.append("_SPHARMMedialMesh.vtk");

		vtkwriter->SetFileName(outFileName.c_str() );
		vtkwriter->Write();
	 }
	 
     vtkSmartPointer<vtkPolyData> medialAxis = vtkSmartPointer<vtkPolyData>::New();
	 medialAxis = medialmeshsrc->GetOutputMedialAxis();

	 
	 double* Theta;
	 double* Radius;
	 double* Area;
	 Theta=medialmeshsrc->GetTheta();
	 Radius=medialmeshsrc->GetRadius();
	 Area=medialmeshsrc->GetArea();
	 
	 outFileName.erase();
	 outFileName.append(base_string);
     outFileName.append("_MedialAxisScalars.csv");
	 
	 std::ofstream ScalarFile(outFileName.c_str());
	 ScalarFile<<"Theta,Radius,Area"<<std::endl;
     vtkSmartPointer<vtkDoubleArray> theta = vtkSmartPointer<vtkDoubleArray>::New();
     theta->SetName("_theta");
     vtkSmartPointer<vtkDoubleArray> radius = vtkSmartPointer<vtkDoubleArray>::New();
     radius->SetName("_radius");
     vtkSmartPointer<vtkDoubleArray> area = vtkSmartPointer<vtkDoubleArray>::New();
     area->SetName("_area");

	 for(int i=0; i<thetaIteration; i++)
         {
             ScalarFile<<Theta[i]<<","<<Radius[i]<<","<<Area[i]<<std::endl;
             theta->InsertNextValue(Theta[i]);
             radius->InsertNextValue(Radius[i]);
             area->InsertNextValue(Area[i]);
             medialAxis->GetPointData()->AddArray(theta);
             medialAxis->GetPointData()->AddArray(radius);
             medialAxis->GetPointData()->AddArray(area);
         }
         ScalarFile.close();

         vtkwriter->SetInputData(medialAxis );
         outFileName.erase();
         outFileName.append(base_string);
         outFileName.append("_SPHARMMedialAxis.vtk");

         vtkwriter->SetFileName(outFileName.c_str() );
         vtkwriter->Write();

         //
         itkMeshTovtkPolyData ITKVTKConverter;
         ITKVTKConverter.SetInput( meshSH );

         if( writePara )
         {
           if( debug )
           {
             std::cout << "writing para mesh data" << std::endl;
           }
           MeshType *             _paraMesh = meshsrc->GetOutputParaMesh();
           PointsContainerPointer paraPoints = _paraMesh->GetPoints();

           outFileName.erase();
           outFileName.append(base_string);
           outFileName.append("_para.vtk");

           itkMeshTovtkPolyData ITKVTKConverter4;
           ITKVTKConverter4.SetInput(_paraMesh);
           vtkwriter->SetInputData(ITKVTKConverter4.GetOutput() );
           vtkwriter->SetFileName(outFileName.c_str() );
           vtkwriter->Write();

           // DIMENSION = 1 ; NUMBER_OF_POINTS ; TYPE = Scalar
           {
               // write phi
               outFileName.erase();
               outFileName.append(base_string);
               outFileName.append("_paraPhi.txt");
               std::ofstream efile(outFileName.c_str(), std::ios::out);
               if( !efile )
               {
                   std::cerr << "Error: open of file \"" << outFileName  << "\" failed." << std::endl;
                   exit(-1);
               }
               efile.precision(10);
               int nvert = paraPoints->Size();
               efile << "NUMBER_OF_POINTS = " << nvert << std::endl;
               efile << "DIMENSION = 1" << std::endl;
               efile << "TYPE = Scalar" << std::endl;
               vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
               array->SetName("_paraPhi");

               for( int i = 0; i < nvert; i++ )
               {
                   PointType curPoint =  paraPoints->GetElement(i);
                   double    phi = atan2(curPoint[1], curPoint[0]) + M_PI; // 0 .. 2 * M_PI
                   efile << phi << endl;
                   // if (i != nvert - 1) efile << " ";
                   array->InsertNextValue(phi);
               }
               efile.close();
               ITKVTKConverter.GetOutput()->GetPointData()->AddArray(array);
           }
           {
               // write phi half
               outFileName.erase();
               outFileName.append(base_string);
               outFileName.append("_paraPhiHalf.txt");
               std::ofstream efile(outFileName.c_str(), std::ios::out);
               if( !efile )
               {
                   std::cerr << "Error: open of file \"" << outFileName  << "\" failed." << std::endl;
                   exit(-1);
               }
               efile.precision(10);

               int nvert = paraPoints->Size();
               efile << "NUMBER_OF_POINTS = " << nvert << std::endl;
               efile << "DIMENSION = 1" << std::endl;
               efile << "TYPE = Scalar" << std::endl;
               vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
               array->SetName("_paraPhiHalf");
               for( int i = 0; i < nvert; i++ )
               {
                   PointType curPoint =  paraPoints->GetElement(i);
                   double    phi = atan2(curPoint[1], curPoint[0]) + M_PI;
                   if( phi > M_PI )
                   {
                       phi = 2 * M_PI - phi;       // 0 .. M_PI ..0
                   }
                   efile << phi << endl;
                   // if (i != nvert - 1) efile << " ";
                   array->InsertNextValue(phi);
               }
               efile.close();
               ITKVTKConverter.GetOutput()->GetPointData()->AddArray(array);
           }
           {
               // write theta
               outFileName.erase();
               outFileName.append(base_string);
               outFileName.append("_paraTheta.txt");
               std::ofstream efile(outFileName.c_str(), std::ios::out);
               if( !efile )
                 {
                 std::cerr << "Error: open of file \"" << outFileName << "\" failed." << std::endl;
                 exit(-1);
                 }
               efile.precision(10);
               int nvert = paraPoints->Size();
               efile << "NUMBER_OF_POINTS = " << nvert << std::endl;
               efile << "DIMENSION = 1" << std::endl;
               efile << "TYPE = Scalar" << std::endl;
               vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
               array->SetName("_paraTheta");
               for( int i = 0; i < nvert; i++ )
               {
                   PointType curPoint =  paraPoints->GetElement(i);
                   double    curTheta = atan(curPoint[2] / sqrt(curPoint[0] * curPoint[0] + curPoint[1] * curPoint[1]) ) + M_PI_2;
                   //0 .. M_PI
                   efile << curTheta << endl;
                   // if (i != nvert - 1) efile << " ";
                   array->InsertNextValue(curTheta);
               }
               efile.close();
               ITKVTKConverter.GetOutput()->GetPointData()->AddArray(array);
           }
           {
               // write theta/M_PI * (phi/M_PI/2 + 1)
               outFileName.erase();
               outFileName.append(base_string);
               outFileName.append("_paraMix.txt");
               std::ofstream efile(outFileName.c_str(), std::ios::out);
               if( !efile )
               {
                   std::cerr << "Error: open of file \"" << outFileName << "\" failed." << std::endl;
                   exit(-1);
               }
               efile.precision(10);

               int nvert = paraPoints->Size();
               efile << "NUMBER_OF_POINTS = " << nvert << std::endl;
               efile << "DIMENSION = 1" << std::endl;
               efile << "TYPE = Scalar" << std::endl;
               vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
               array->SetName("_paraMix");
               for( int i = 0; i < nvert; i++ )
               {
                   PointType curPoint =  paraPoints->GetElement(i);
                   double    phi = atan2(curPoint[1], curPoint[0]) + M_PI;
                   if( phi > M_PI )
                   {
                       phi = 2 * M_PI - phi;            // 0 .. M_PI ..0
                   }
                   double curTheta = atan(curPoint[2] / sqrt(curPoint[0] * curPoint[0] + curPoint[1] * curPoint[1]) ) + M_PI_2;
                   // 0 .. M_PI
                   efile << curTheta / M_PI * ( phi / M_PI + 1) << endl;
                   // if (i != nvert - 1) efile << " ";
                   array->InsertNextValue(curTheta/ M_PI * ( phi / M_PI + 1));
               }
               efile.close();
               ITKVTKConverter.GetOutput()->GetPointData()->AddArray(array);
           }
         }
         outFileName.erase();
         outFileName.append(base_string);
         outFileName.append("_SPHARM.vtk");
         vtkwriter->SetFileName(outFileName.c_str() );

         vtkwriter->SetInputData(ITKVTKConverter.GetOutput());
         vtkwriter->SetFileName(outFileName.c_str());
         vtkwriter->Write();

         if( debug )
           {
           std::cout << "saving par aligned coefs" << std::endl;
           }

         // save coefficients too
         outFileName.erase();
         outFileName.append(base_string);
         outFileName.append("_SPHARM.coef");
         typedef neurolib::SphericalHarmonicCoefficientFileWriter CoefWriterType;
         CoefWriterType::Pointer coefwriter = CoefWriterType::New();

         coefwriter->SetFileName(outFileName.c_str() );
         coefwriter->SetInput(coeflist);
         coefwriter->Update();

    if( !NoParaAlignFlag )
      {
      if( debug )
        {
        std::cout << "generating ellipse aligned data" << std::endl;
        }
      // get EllipseAlign
      neurolib::SphericalHarmonicSpatialObject::CoefListType * ellipsoidAlign;
      ellipsoidAlign = spharmFilter->GetEllipseAlignCoef();

      if( debug )
        {
        std::cout << "writing ellipse aligned data" << std::endl;
        }
      neurolib::SphericalHarmonicMeshSource::Pointer meshEllisrc = neurolib::SphericalHarmonicMeshSource::New();
      meshEllisrc->SetCoefs(*ellipsoidAlign);
      meshEllisrc->SetLevel(subdivLevel);
      meshEllisrc->Update();
      MeshType * ellipseMesh = meshEllisrc->GetOutput();

      itkMeshTovtkPolyData ITKVTKConverter2;
      ITKVTKConverter2.SetInput(ellipseMesh);

      if( writePara )
      {
          MeshType *             _paraMesh = meshsrc->GetOutputParaMesh();
          PointsContainerPointer paraPoints = _paraMesh->GetPoints();
          // DIMENSION = 1 ; NUMBER_OF_POINTS ; TYPE = Scalar
            {
              // write phi
              int nvert = paraPoints->Size();
              vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
              array->SetName("_paraPhi");

              for( int i = 0; i < nvert; i++ )
              {
                  PointType curPoint =  paraPoints->GetElement(i);
                  double    phi = atan2(curPoint[1], curPoint[0]) + M_PI; // 0 .. 2 * M_PI
                  array->InsertNextValue(phi);
              }
              ITKVTKConverter2.GetOutput()->GetPointData()->AddArray(array);
            }
            {
              // write phi half
                int nvert = paraPoints->Size();
                vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
                array->SetName("_paraPhiHalf");
                for( int i = 0; i < nvert; i++ )
                {
                  PointType curPoint =  paraPoints->GetElement(i);
                  double    phi = atan2(curPoint[1], curPoint[0]) + M_PI;
                  if( phi > M_PI )
                  {
                    phi = 2 * M_PI - phi;       // 0 .. M_PI ..0
                  }
                  array->InsertNextValue(phi);
                }
                ITKVTKConverter2.GetOutput()->GetPointData()->AddArray(array);
            }
            {
                // write theta
                int nvert = paraPoints->Size();
                vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
                array->SetName("_paraTheta");
                for( int i = 0; i < nvert; i++ )
                {
                  PointType curPoint =  paraPoints->GetElement(i);
                  double    curTheta = atan(curPoint[2] / sqrt(curPoint[0] * curPoint[0] + curPoint[1] * curPoint[1]) ) + M_PI_2;
                  //0 .. M_PI
                  array->InsertNextValue(curTheta);
                }
                ITKVTKConverter2.GetOutput()->GetPointData()->AddArray(array);
            }
            {
                // write theta/M_PI * (phi/M_PI/2 + 1)
                int nvert = paraPoints->Size();
                vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
                array->SetName("_paraMix");
                for( int i = 0; i < nvert; i++ )
                {
                  PointType curPoint =  paraPoints->GetElement(i);
                  double    phi = atan2(curPoint[1], curPoint[0]) + M_PI;
                  if( phi > M_PI )
                    {
                    phi = 2 * M_PI - phi;            // 0 .. M_PI ..0
                    }
                  double curTheta = atan(curPoint[2] / sqrt(curPoint[0] * curPoint[0] + curPoint[1] * curPoint[1]) ) + M_PI_2;
                  // 0 .. M_PI
                  array->InsertNextValue(curTheta/ M_PI * ( phi / M_PI + 1));
                }
                ITKVTKConverter2.GetOutput()->GetPointData()->AddArray(array);
            }
      }

      vtkwriter->SetInputData(ITKVTKConverter2.GetOutput() );
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_SPHARM_ellalign.vtk");
      vtkwriter->SetFileName(outFileName.c_str() );
      vtkwriter->Write();

      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_SPHARM_ellalign.coef");
      coefwriter->SetFileName(outFileName.c_str() );
      coefwriter->SetInput(*ellipsoidAlign);
      coefwriter->Update();
      }

    if( regTemplateFileOn || regParaTemplateFileOn )
      {
        std::cout<<"regTemplateFileOn: "<<regTemplateFileOn<<std::endl;
        std::cout<<"regParaTemplateFileOn: "<<regParaTemplateFileOn<<std::endl;

      if( debug )
        {
        std::cout << "writing proc aligned data" << std::endl;
        }
      // load RegTemplateMesh
      MeshType::Pointer RegTemplateMesh = MeshType::New();

      if( regTemplateFileOn )
        {
        vtkSmartPointer<vtkPolyData> convRegTemplateMesh = vtkSmartPointer<vtkPolyData>::New();
        convRegTemplateMesh = ReadPolyData( regTemplateFile.c_str() );

        vtkPolyDataToitkMesh VTKITKConverter3;
        VTKITKConverter3.SetInput( convRegTemplateMesh );
        RegTemplateMesh = VTKITKConverter3.GetOutput();

        }

      ProcrustesFilterType::Pointer procrustesFilter = ProcrustesFilterType::New();
      procrustesFilter->SetNumberOfInputs(2);
      if( regTemplateFileOn )
        {
        procrustesFilter->SetInput(0, RegTemplateMesh);
        }

      if( regParaTemplateFileOn )
        {
        procrustesFilter->SetInput(0, regParaTemplateMesh);
        }
      procrustesFilter->SetInput(1, meshSH);
      procrustesFilter->SetUseInitialAverageOff();
      procrustesFilter->SetUseNormalizationOff();
      procrustesFilter->SetUseScalingOff();
      procrustesFilter->SetUseSingleIterationOn();
      procrustesFilter->Update();

      if( procrustesTransformationOutputOn )
        {

        // there are two transforms here
        // but because single iteration is on
        // and useInitialAverage is off
        // the transformation of the second one
        // will be the transformation that maps back to the first one
        // this is what we are interested in here

        typedef itk::TransformFileWriter TransformFileWriterType;
        TransformFileWriterType::Pointer transformFileWriter = TransformFileWriterType::New();
        transformFileWriter->SetInput( procrustesFilter->GetTransform( 1 ) );
        transformFileWriter->SetFileName( procrustesTransformationFile );

        try
          {
          transformFileWriter->Update();
          }
        catch( itk::ExceptionObject & e )
          {
          std::cerr << "Non-fatal exception while trying to write the Procrustes affine transform" << std::endl;
          std::cerr << e << std::endl;
          }

        }

      MeshType::Pointer RegisteredMesh = procrustesFilter->GetOutput(1);

      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_SPHARM_procalign.vtk");
      std::cout << "procalign" << std::endl;
      itkMeshTovtkPolyData ITKVTKConverter3;
      ITKVTKConverter3.SetInput(RegisteredMesh);
      vtkwriter->SetInputData(ITKVTKConverter3.GetOutput() );
      vtkwriter->SetFileName(outFileName.c_str() );
      vtkwriter->Write();

      }
    }

  catch( itk::ExceptionObject &e )
    {
    e.Print(std::cout);
    return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;

}

vtkSmartPointer<vtkPolyData> ReadPolyData(std::string filePath)
{
    size_t found = filePath.rfind(".vtp");
    if (found != std::string::npos)
    {
        vtkSmartPointer<vtkXMLPolyDataReader> VTPreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        VTPreader->SetFileName(filePath.c_str());
        VTPreader->Update();
        return VTPreader->GetOutput();
    }
    found = filePath.rfind(".vtk");
    if ( found != std::string::npos )
    {
        vtkSmartPointer<vtkPolyDataReader> VTKreader = vtkSmartPointer<vtkPolyDataReader>::New();
        VTKreader->SetFileName(filePath.c_str());
        VTKreader->Update();
        return VTKreader->GetOutput();
    }
    return NULL;
}
