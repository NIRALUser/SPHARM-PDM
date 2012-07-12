/*=========================================================================

  Author: Christine Xu

=========================================================================*/

#include <itkTriangleCell.h>

#include "SphericalHarmonicMedialAxisMeshSource.h"

#include <math.h>
#include <stdio.h>

#include <iostream>

namespace neurolib
{

#ifndef M_PI
#define M_PI 3.1415926535897932
#endif
#ifndef M_PI_2
#define M_PI_2 1.5707963267948966
#endif

typedef itk::TriangleCell<SphericalHarmonicMedialAxisMeshSource::CellType> TriangleType;
typedef itk::LineCell<SphericalHarmonicMedialAxisMeshSource::CellType> LineType;

SphericalHarmonicMedialAxisMeshSource::SphericalHarmonicMedialAxisMeshSource()
{
	m_Dimension = 3;
	
	m_Degree = 0;
	m_FromL = 0;
	m_ToL = 0;

	//	GenerateData();
}

SphericalHarmonicMedialAxisMeshSource::~SphericalHarmonicMedialAxisMeshSource()
{
}

void SphericalHarmonicMedialAxisMeshSource::SetDegree(unsigned int d)
{
	if( m_Coefs.empty() || d > floor(sqrt( (double) m_Coefs.size() ) ) - 1 )
	{
		throw SphericalHarmonicPolynomialException(
				__FILE__, __LINE__, "The maximum degree of the coefficient exceeds the size of coefficient list.");
	}
	m_Degree = d;
	m_ToL = d;
}

void SphericalHarmonicMedialAxisMeshSource::SetCoefs(SphericalHarmonicMedialAxisMeshSource::CoefListType& coeflist)
{
	m_Coefs = coeflist;
	m_Degree =  (int) floor(sqrt( (double) coeflist.size() ) ) - 1;
	m_ToL  = m_Degree;
	// std::cout<<"m_Degree = "<<m_Degree <<std::endl;
}

void SphericalHarmonicMedialAxisMeshSource::SetThetaPhiIteration(int theta, int phi)
{
	m_theta=theta;
	m_phi=phi;
	//	std::cout << "here" << std::endl;

	///////////
}

void SphericalHarmonicMedialAxisMeshSource::GenerateData()
{

  // std::cout << "IN GENERATE DATA AAAAAAA " << std::endl;
        
	if( m_Coefs.empty() )
	{
		throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, "Coefficients mustn't be empty.");
	}
	if( m_FromL > m_ToL )
	{
		throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
				"The starting degree should be smaller or equal to the ending degree.");
	}
	if( m_ToL > m_Degree )
	{
		throw SphericalHarmonicPolynomialException(__FILE__, __LINE__,
				"The evalueated degree mustn't exceed the size of the coefficients.");
	}
	//	std::cout << "plop" << std::endl;
	int n_vert=m_phi * m_theta;
// 	int n_triag;

  // Allocate datas
// 	Point3* all_vert = new Point3[m_phi * m_theta];
// 	Point3* all_triangs = new Point3[n_triag * 3]; // all possible vertices in triangs
// 	int *   triangs = new int[3 * n_triag];
// 	Point3* mesh = new Point3[m_phi * m_theta];
	m_PhiTable = new double[m_phi];
	m_ThetaTable = new double[m_theta];

	//	std::cout << m_theta << " " << m_phi << std::endl;

	m_ThetaPhiTable = new double[m_phi*m_theta*2];
	Point3* vertex = new Point3[n_vert];
	Point3* medialVertex = new Point3[m_theta];
	m_radius = new double[m_theta];
	double* surface = new double[m_theta];
	double volume=0;

	setPhiThetaTable();

	SphericalHarmonicPolynomial<3> SPHARM;
	try
	{
		SPHARM.SetCoefs(m_Coefs);
		SPHARM.SetDegree(m_Degree);
    // Calculate mesh
		for(int i=0; i<m_theta; i++)
		{
			Point3 p1;
			p1[0]=0;
			p1[1]=0;
			p1[2]=0;
			for(int j=0; j<m_phi; j++)
			{
				SPHARM.Evaluate(m_FromL, m_ToL, m_ThetaPhiTable[(i*m_phi+j)*2+1], m_ThetaPhiTable[(i*m_phi+j)*2], vertex[i*m_phi + j]);
				
				//Averaging each points
				p1[0]+=vertex[i*m_phi+j][0];
				p1[1]+=vertex[i*m_phi+j][1];
				p1[2]+=vertex[i*m_phi+j][2];
			}
			medialVertex[i][0]=p1[0]/m_phi;
			medialVertex[i][1]=p1[1]/m_phi;
			medialVertex[i][2]=p1[2]/m_phi;
		}
	}
	catch( SphericalHarmonicPolynomialException ex )
	{
		throw SphericalHarmonicPolynomialException(__FILE__, __LINE__, ex.GetDescription() );
	}
	
	for(int i=0; i<m_theta; i++)
	{
		m_radius[i]=0;
		surface[i]=0;
		for(int j=0; j<m_phi; j++)
		{
			double distance=0,triangle=0,x,y,z,temp;
			double* p1=new double[3];
			double* p2=new double[3];
			p1[0]=vertex[i*m_phi+j][0];
			p1[1]=vertex[i*m_phi+j][1];
			p1[2]=vertex[i*m_phi+j][2];
			p2[0]=vertex[i*m_phi+(j+1)%m_phi][0];
			p2[1]=vertex[i*m_phi+(j+1)%m_phi][1];
			p2[2]=vertex[i*m_phi+(j+1)%m_phi][2];
			
			x=p1[0]-medialVertex[i][0];
			y=p1[1]-medialVertex[i][1];
			z=p1[2]-medialVertex[i][2];
			distance=sqrt(x*x+y*y+z*z);
			m_radius[i]+=distance;
			
			temp=medialVertex[i][0]*p2[1]-medialVertex[i][0]*p1[1]+p1[0]*medialVertex[i][1]-p1[0]*p2[1]+p2[0]*p1[1]-p2[0]*medialVertex[i][1];
			temp=temp*temp;
			triangle+=temp;
			temp=medialVertex[i][1]*p2[2]-medialVertex[i][1]*p1[2]+p1[1]*medialVertex[i][2]-p1[1]*p2[2]+p2[1]*p1[2]-p2[1]*medialVertex[i][2];
			temp=temp*temp;
			triangle+=temp;
			temp=medialVertex[i][2]*p2[0]-medialVertex[i][2]*p1[0]+p1[2]*medialVertex[i][0]-p1[2]*p2[0]+p2[2]*p1[0]-p2[2]*medialVertex[i][0];
			temp=temp*temp;
			triangle+=temp;
			triangle=sqrt(triangle)/2;
			
			surface[i]+=triangle;
		}
		m_radius[i]/=m_phi;
	}
	
	//Volume
	for(int i=0; i<m_theta-1; i++)
	{
		double h, x, y, z;
		x=medialVertex[i+1][0]-medialVertex[i][0];
		y=medialVertex[i+1][1]-medialVertex[i][1];
		z=medialVertex[i+1][2]-medialVertex[i][2];
		h=sqrt(x*x+y*y+z*z);
		volume+=(surface[i]+surface[i+1])*h/2;
	}
	
	MeshType::Pointer outputMesh = this->GetOutput();
	m_outputMedialAxis = vtkSmartPointer<vtkPolyData>::New();

	PointsContainerPointer points = PointsContainer::New();
	for( int i = 0; i < n_vert; i++ )
		points->InsertElement(i, PointType(vertex[i]) );
	outputMesh->SetPoints(points);
	
	vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> Lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyLine> Line = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkFloatArray> Radius = vtkSmartPointer<vtkFloatArray>::New();
	Radius->SetNumberOfComponents(m_theta);
	Radius->SetName("Radius");
	vtkSmartPointer<vtkFloatArray> RadiusScaled = vtkSmartPointer<vtkFloatArray>::New();
	RadiusScaled->SetNumberOfComponents(m_theta);
	RadiusScaled->SetName("Radius_Scaled");
	vtkSmartPointer<vtkFloatArray> Theta = vtkSmartPointer<vtkFloatArray>::New();
	Theta->SetNumberOfComponents(m_theta);
	Theta->SetName("Theta");
	vtkSmartPointer<vtkFloatArray> ThetaScaled = vtkSmartPointer<vtkFloatArray>::New();
	ThetaScaled->SetNumberOfComponents(m_theta);
	ThetaScaled->SetName("Theta_Scaled");
	vtkSmartPointer<vtkFloatArray> Surface = vtkSmartPointer<vtkFloatArray>::New();
	Surface->SetNumberOfComponents(m_theta);
	Surface->SetName("Surface");
	vtkSmartPointer<vtkFloatArray> SurfaceScaled = vtkSmartPointer<vtkFloatArray>::New();
	SurfaceScaled->SetNumberOfComponents(m_theta);
	SurfaceScaled->SetName("Surface_Scaled");
	vtkSmartPointer<vtkFloatArray> Volume = vtkSmartPointer<vtkFloatArray>::New();
	Volume->SetNumberOfComponents(1);
	Volume->SetName("Volume");
	
	
	double* ThetaBounds=new double[2];
	GetBounds(m_ThetaTable,ThetaBounds,m_theta);
	double* RadiusBounds=new double[2];
	GetBounds(m_radius,RadiusBounds,m_theta);
	double* SurfaceBounds=new double[2];
	GetBounds(surface,SurfaceBounds,m_theta);
	Line->GetPointIds()->SetNumberOfIds(m_theta);
	
	for(int i=0; i<m_theta; i++)
	{
		Points->InsertNextPoint(medialVertex[i]);
		Line->GetPointIds()->SetId(i,i);
		Theta->InsertComponent(0,i,m_ThetaTable[i]);
		ThetaScaled->InsertComponent(0,i,100*(m_ThetaTable[i]-ThetaBounds[0])/(ThetaBounds[1]-ThetaBounds[0]));
		Radius->InsertComponent(0,i,m_radius[i]);
		RadiusScaled->InsertComponent(0,i,100*(m_radius[i]-RadiusBounds[0])/(RadiusBounds[1]-RadiusBounds[0]));
		Surface->InsertComponent(0,i,surface[i]);
		SurfaceScaled->InsertComponent(0,i,100*(surface[i]-SurfaceBounds[0])/(SurfaceBounds[1]-SurfaceBounds[0]));
	}
	
	Volume->InsertComponent(0,0,volume);
	Lines->InsertNextCell(Line);
	m_outputMedialAxis->SetPoints(Points);
	m_outputMedialAxis->SetLines(Lines);
	m_outputMedialAxis->GetPointData()->AddArray(Theta);
	m_outputMedialAxis->GetPointData()->AddArray(ThetaScaled);
	m_outputMedialAxis->GetPointData()->AddArray(Radius);
	m_outputMedialAxis->GetPointData()->AddArray(RadiusScaled);
	m_outputMedialAxis->GetPointData()->AddArray(Surface);
	m_outputMedialAxis->GetPointData()->AddArray(SurfaceScaled);
	m_outputMedialAxis->GetPointData()->AddArray(Volume);

  /**
	 * Specify the method used for allocating cells
	*/
	outputMesh->SetCellsAllocationMethod( MeshType::CellsAllocatedDynamicallyCellByCell );
	
	
	for(int i=0; i<m_theta-1; i++ )
	{
		for(int j=0; j<m_phi-1; j++)
		{
			{
				CellType::CellAutoPointer cellpointer;
				cellpointer.TakeOwnership(new TriangleType);
			/**
				* Assign the points to the tetrahedron through their identifiers.
			*/
				unsigned long triPoints[3];
				
				triPoints[0] = i*m_phi+j;
				triPoints[2] = (i+1)*m_phi+j+1;
				triPoints[1] = (i+1)*m_phi+j;
				cellpointer->SetPointIds(triPoints);
				outputMesh->SetCell(2*(i*(m_phi-1)+j), cellpointer);
				
				CellType::CellAutoPointer cellpointer2;
				cellpointer2.TakeOwnership(new TriangleType);
				
				triPoints[0] = i*m_phi+j;
				triPoints[2] = (i)*m_phi+j+1;
				triPoints[1] = (i+1)*m_phi+j+1;
				cellpointer2->SetPointIds(triPoints);
				outputMesh->SetCell(2*(i*(m_phi-1)+j)+1, cellpointer2);
			}
		}

	}
	
	delete[] surface;
	delete[] vertex;
	delete[] medialVertex;

	//	std::cout << "Number of elements in Radius " <<  (sizeof(m_radius) / sizeof(m_radius[0])) << std::endl;
	//std::cout << "Number of elements in Theta " <<  (sizeof(m_ThetaTable) / sizeof(m_ThetaTable[0])) << std::endl;
}

void SphericalHarmonicMedialAxisMeshSource::setPhiThetaTable()
{
	double step=2*M_PI/(m_phi-1), currentValue=0;
	for(int i=0; i<m_phi; i++)
	{
		m_PhiTable[i]=currentValue;
		currentValue+=step;
	}
	
	step=M_PI/(m_theta-1);
	currentValue=0;
	for(int i=0; i<m_theta; i++)
	{
		m_ThetaTable[i]=currentValue;
		currentValue+=step;
	}

	for(int i=0; i<m_theta; i++)
	{
		for(int j=0; j<m_phi; j++)
		{
			m_ThetaPhiTable[(i*m_phi+j)*2]=acos(-sin(m_ThetaTable[i])*cos(m_PhiTable[j]));
			if(sin(m_ThetaTable[i])*sin(m_ThetaTable[i])*cos(m_PhiTable[j])*cos(m_PhiTable[j])==1 || cos(m_ThetaTable[i])/sqrt(1-sin(m_ThetaTable[i])*sin(m_ThetaTable[i])*cos(m_PhiTable[j])*cos(m_PhiTable[j]))>=1)
				m_ThetaPhiTable[(i*m_phi+j)*2+1]=0;
			else if(cos(m_ThetaTable[i])/sqrt(1-sin(m_ThetaTable[i])*sin(m_ThetaTable[i])*cos(m_PhiTable[j])*cos(m_PhiTable[j]))<=-1)
				m_ThetaPhiTable[(i*m_phi+j)*2+1]=M_PI;
			else
			{
				if(sin(m_ThetaTable[i])*sin(m_PhiTable[j])<0)
					m_ThetaPhiTable[(i*m_phi+j)*2+1]=2*M_PI-acos(cos(m_ThetaTable[i])/sqrt(1-sin(m_ThetaTable[i])*sin(m_ThetaTable[i])*cos(m_PhiTable[j])*cos(m_PhiTable[j])));
				else
					m_ThetaPhiTable[(i*m_phi+j)*2+1]=acos(cos(m_ThetaTable[i])/sqrt(1-sin(m_ThetaTable[i])*sin(m_ThetaTable[i])*cos(m_PhiTable[j])*cos(m_PhiTable[j])));
			}
		}
	}
}

void SphericalHarmonicMedialAxisMeshSource::GetBounds(double Vector[], double Bounds[], int Size)
{
	double min=100000,max=-1000000;
	for(int i=0; i<Size; i++)
	{
		if(Vector[i]<min)
			min=Vector[i];
		if(Vector[i]>max)
			max=Vector[i];
	}
	Bounds[0]=min;
	Bounds[1]=max;
}

// void SphericalHarmonicMeshSource::set_up_icosahedron_triangs(Point3* all_vert,
// 		Point3* all_triangs,
// 		int subdiv,
// 		int n_vert,
// 		int m_phi,
// 		int m_theta,
// 		double *icos,
// 		int n_triangs,
// 		int *triangs)
// {
// 	int    i, n, m, k, numtriags;
// 	double x1, x2, y1, y2, z1, z2, x3, y3, z3;
// 	double dx12, dy12, dz12, dx23, dy23, dz23;
// 	double length;
// 
// 	double epsilon = 0.00001; // machine epsilon??
// 
// 	memcpy(all_vert, vert, 12 * sizeof(Point3) );
// 
//   // std::cout<<"after memcpy"<<std::endl;
// 
// 	k = 12;
// 	for( i = 0; i < 30; i++ )
// 	{
// 		x1 = vert[edge[i][0]][0];
// 		y1 = vert[edge[i][0]][1];
// 		z1 = vert[edge[i][0]][2];
// 		x2 = vert[edge[i][1]][0];
// 		y2 = vert[edge[i][1]][1];
// 		z2 = vert[edge[i][1]][2];
// 		dx12 = (x2 - x1) / subdiv;
// 		dy12 = (y2 - y1) / subdiv;
// 		dz12 = (z2 - z1) / subdiv;
// 		for( n = 1; n < subdiv; n++ )
// 		{
// 			all_vert[k][0] = x1 + n * dx12;
// 			all_vert[k][1] = y1 + n * dy12;
// 			all_vert[k][2] = z1 + n * dz12;
// 			length = sqrt( (double) all_vert[k][0] * all_vert[k][0]
// 					+ all_vert[k][1] * all_vert[k][1]
// 					+ all_vert[k][2] * all_vert[k][2]);
// 			all_vert[k][0] /= length;
// 			all_vert[k][1] /= length;
// 			all_vert[k][2] /= length;
// 			k++;
// 		}
// 	}
// 
// 	if( subdiv > 2 )
// 	{
// 		for( i = 0; i < 20; i++ )
// 		{
// 			x1 = vert[triang[i][0]][0];
// 			y1 = vert[triang[i][0]][1];
// 			z1 = vert[triang[i][0]][2];
// 			x2 = vert[triang[i][1]][0];
// 			y2 = vert[triang[i][1]][1];
// 			z2 = vert[triang[i][1]][2];
// 			x3 = vert[triang[i][2]][0];
// 			y3 = vert[triang[i][2]][1];
// 			z3 = vert[triang[i][2]][2];
// 			dx12 = (x2 - x1) / subdiv;
// 			dy12 = (y2 - y1) / subdiv;
// 			dz12 = (z2 - z1) / subdiv;
// 			dx23 = (x3 - x2) / subdiv;
// 			dy23 = (y3 - y2) / subdiv;
// 			dz23 = (z3 - z2) / subdiv;
// 
// 			n = 1;
// 
// 			do
// 			{
// 				for( m = 1; m <= n; m++ )
// 				{
// 					all_vert[k][0] = x1 + (n + 1) * dx12 + m * dx23;
// 					all_vert[k][1] = y1 + (n + 1) * dy12 + m * dy23;
// 					all_vert[k][2] = z1 + (n + 1) * dz12 + m * dz23;
// 					length = sqrt( (double) all_vert[k][0] * all_vert[k][0]
// 							+ all_vert[k][1] * all_vert[k][1]
// 							+ all_vert[k][2] * all_vert[k][2]);
// 					all_vert[k][0] /= length;
// 					all_vert[k][1] /= length;
// 					all_vert[k][2] /= length;
// 					k++;
// 				}
// 				n++;
// 			}
// 			while( n <= (subdiv - 2) );
// 		}
// 	}
// 	numtriags = 0;
// 
//   // std::cout<<"before get triangulation"<<std::endl;
//   // std::cout<<n_triangs<<std::endl;
// 
//   // get triangulation
// 	if( subdiv > 1 )
// 	{
// 		for( i = 0; i < 20; i++ )
// 		{
// 			x1 = vert[triang[i][0]][0];
// 			y1 = vert[triang[i][0]][1];
// 			z1 = vert[triang[i][0]][2];
// 			x2 = vert[triang[i][1]][0];
// 			y2 = vert[triang[i][1]][1];
// 			z2 = vert[triang[i][1]][2];
// 			x3 = vert[triang[i][2]][0];
// 			y3 = vert[triang[i][2]][1];
// 			z3 = vert[triang[i][2]][2];
// 			dx12 = (x2 - x1) / subdiv;
// 			dy12 = (y2 - y1) / subdiv;
// 			dz12 = (z2 - z1) / subdiv;
// 			dx23 = (x3 - x2) / subdiv;
// 			dy23 = (y3 - y2) / subdiv;
// 			dz23 = (z3 - z2) / subdiv;
// 
// 			n = 1;
// 
// 			do
// 			{
// 				for( m = 1; m <= n; m++ )
// 				{
//           // Draw lower triangle
// 					all_triangs[numtriags][0] = x1 + n * dx12 + m * dx23;
// 					all_triangs[numtriags][1] = y1 + n * dy12 + m * dy23;
// 					all_triangs[numtriags][2] = z1 + n * dz12 + m * dz23;
// 					length = sqrt( (double) all_triangs[numtriags][0] * all_triangs[numtriags][0]
// 							+ all_triangs[numtriags][1] * all_triangs[numtriags][1]
// 							+ all_triangs[numtriags][2] * all_triangs[numtriags][2]);
// 					all_triangs[numtriags][0] /= length;
// 					all_triangs[numtriags][1] /= length;
// 					all_triangs[numtriags][2] /= length;
// 					numtriags++;
// 					all_triangs[numtriags][0] = x1 + (n - 1) * dx12 + (m - 1) * dx23;
// 					all_triangs[numtriags][1] = y1 + (n - 1) * dy12 + (m - 1) * dy23;
// 					all_triangs[numtriags][2] = z1 + (n - 1) * dz12 + (m - 1) * dz23;
// 					length = sqrt( (double) all_triangs[numtriags][0] * all_triangs[numtriags][0]
// 							+ all_triangs[numtriags][1] * all_triangs[numtriags][1]
// 							+ all_triangs[numtriags][2] * all_triangs[numtriags][2]);
// 					all_triangs[numtriags][0] /= length;
// 					all_triangs[numtriags][1] /= length;
// 					all_triangs[numtriags][2] /= length;
// 					numtriags++;
// 					all_triangs[numtriags][0] = x1 + n * dx12 + (m - 1) * dx23;
// 					all_triangs[numtriags][1] = y1 + n * dy12 + (m - 1) * dy23;
// 					all_triangs[numtriags][2] = z1 + n * dz12 + (m - 1) * dz23;
// 					length = sqrt( (double) all_triangs[numtriags][0] * all_triangs[numtriags][0]
// 							+ all_triangs[numtriags][1] * all_triangs[numtriags][1]
// 							+ all_triangs[numtriags][2] * all_triangs[numtriags][2]);
// 					all_triangs[numtriags][0] /= length;
// 					all_triangs[numtriags][1] /= length;
// 					all_triangs[numtriags][2] /= length;
// 					numtriags++;
// 					if( m != n )
// 					{
//             // Draw lower left triangle
// 						all_triangs[numtriags][0] = x1 + n * dx12 + m * dx23;
// 						all_triangs[numtriags][1] = y1 + n * dy12 + m * dy23;
// 						all_triangs[numtriags][2] = z1 + n * dz12 + m * dz23;
// 						length = sqrt( (double) all_triangs[numtriags][0] * all_triangs[numtriags][0]
// 								+ all_triangs[numtriags][1] * all_triangs[numtriags][1]
// 								+ all_triangs[numtriags][2] * all_triangs[numtriags][2]);
// 						all_triangs[numtriags][0] /= length;
// 						all_triangs[numtriags][1] /= length;
// 						all_triangs[numtriags][2] /= length;
// 						numtriags++;
// 						all_triangs[numtriags][0] = x1 + (n - 1) * dx12 + m * dx23;
// 						all_triangs[numtriags][1] = y1 + (n - 1) * dy12 + m * dy23;
// 						all_triangs[numtriags][2] = z1 + (n - 1) * dz12 + m * dz23;
// 						length = sqrt( (double) all_triangs[numtriags][0] * all_triangs[numtriags][0]
// 								+ all_triangs[numtriags][1] * all_triangs[numtriags][1]
// 								+ all_triangs[numtriags][2] * all_triangs[numtriags][2]);
// 						all_triangs[numtriags][0] /= length;
// 						all_triangs[numtriags][1] /= length;
// 						all_triangs[numtriags][2] /= length;
// 						numtriags++;
// 						all_triangs[numtriags][0] = x1 + (n - 1) * dx12 + (m - 1) * dx23;
// 						all_triangs[numtriags][1] = y1 + (n - 1) * dy12 + (m - 1) * dy23;
// 						all_triangs[numtriags][2] = z1 + (n - 1) * dz12 + (m - 1) * dz23;
// 						length = sqrt( (double) all_triangs[numtriags][0] * all_triangs[numtriags][0]
// 								+ all_triangs[numtriags][1] * all_triangs[numtriags][1]
// 								+ all_triangs[numtriags][2] * all_triangs[numtriags][2]);
// 						all_triangs[numtriags][0] /= length;
// 						all_triangs[numtriags][1] /= length;
// 						all_triangs[numtriags][2] /= length;
// 						numtriags++;
// 					}
// 				}
// 				n++;
// 			}
// 			while( n <= subdiv );
// 		}
// 	}
// 
//   // std::cout<<"before indexing of triangs"<<std::endl;
// 
//   // indexing of triangs
// 	if( subdiv == 1 )
// 	{
// 		memcpy(triangs, triang, 20 * 3 * sizeof(int) );
// 		numtriags = 20;
// 	}
// 	else
// 	{
//     // find for every point in triangle list the corresponding index in all_vert
//     // initialize
// 		for( i = 0; i < numtriags; i++ )
// 		{
// 			triangs[i] = -1;
// 		}
//     // find indexes
// 		for( i = 0; i < n_vert; i++ )
// 		{
// 			for( int j = 0; j < numtriags; j++ )
// 			{
// 				if( triangs[j] < 0 )
// 				{
// 					if( (fabs(all_vert[i][0] - all_triangs[j][0]) < epsilon) &&
// 										(fabs(all_vert[i][1] - all_triangs[j][1]) < epsilon) &&
// 										(fabs(all_vert[i][2] - all_triangs[j][2]) < epsilon ) )
// 					{
// 						triangs[j] = i;
// 					}
// 				}
// 			}
// 		}
//     // for(i=0; i<n_vert; i++)
//     //  std::cout<<triangs[3*i]<<","<<triangs[3*i+1]<<","<<triangs[3*i+2]<<std::endl;
// 		for( i = 0; i < numtriags; i++ )
// 		{
// 			if( triangs[i] == -1 )
// 			{
// 				std::cerr << " - " << i << " :" << all_triangs[i][0]
// 						<< "," << all_triangs[i][1] << "," << all_triangs[i][2] << std::endl;
// 			}
// 		}
// 
//     // numtriags is the number of vertices in triangles -> divide it by 3
// 		numtriags = numtriags / 3;
// 	}
//   // std::cout<<"before get phi   heta"<<std::endl;
//   // get phi   eta
// 	for( i = 0; i < n_vert; i++ )
// 	{
// 		icos[2 * i] = atan2(all_vert[i][1], all_vert[i][0]) + M_PI;
// 		icos[2 * i + 1] = atan(all_vert[i][2]
//                            / sqrt( (double) all_vert[i][0] * all_vert[i][0]
// 				+ all_vert[i][1] * all_vert[i][1]) ) + M_PI_2;
//     // std::cout<<icos[2*i]<<","<<icos[2*i+1]<<std::endl;
// 	}
// }

// void SphericalHarmonicMeshSource::interpol_vert(int m_phi,
// 		int m_theta,
// 		Point3 *mesh,
// 		int n_vertex,
// 		double *icos,
// 		Point3 *vertex)
// {
// 	int    i;
// 	double phi, theta;
// 	double x, y;
// 	int    xu, xd, yu, yd;
// 	double xi, yi, zi;
// 	double ksi, eta;
// 
// 	for( i = 0; i < n_vertex; i++ )
// 	{
// 		phi = icos[2 * i]; // phi and theta for every vertex
// 		theta = icos[2 * i + 1];
// 
// 		x = (phi - 1e-5) * ( (float)m_phi) / (2 * M_PI);
// 		xd = (int) x;
// 		xu = xd + 1;
// 		ksi = x - xd;
// 		if( xu >= m_phi )
// 		{
// 			xd = mod(xd, m_phi);
// 			xu = mod(xu, m_phi);
// 		}
// 
// 		y = theta * m_theta / M_PI - 0.5;
// 		yd = (int) y;
// 		yu = yd + 1;
// 		eta = y - yd;
// 		if( yu == m_theta )
// 		{
// 			yu = yu - 1;
// 			xd = mod(xd + m_theta, m_phi);
// 			xu = mod(xu + m_theta, m_phi);
// 		}
// 		if( yd == -1 )
// 		{
// 			yd = 0;
// 			xd = mod(xd + m_theta, m_phi);
// 			xu = mod(xu + m_theta, m_phi);
// 		}
// 
// 		interpol_2d(m_phi, m_theta, xd, xu, yd, yu, ksi, eta, mesh, &xi, &yi, &zi);  // vertex
// 		vertex[i][0] = xi;
// 		vertex[i][1] = yi;
// 		vertex[i][2] = zi;
// 
// 	}
// }
// 
// void SphericalHarmonicMeshSource::interpol_2d(int m_phi,
// 															 int m_theta,
// 															 int xd,
// 															 int xu,
// 															 int yd,
// 															 int yu,
// 															 double ksi,
// 															 double eta,
// 															 Point3 *mesh,
// 															 double *xi,
// 															 double *yi,
// 															 double *zi)
// {
// 	double f00, f10, f01, f11;
// 
// 	f00 = mesh[xd + m_phi * yd][0];
// 	f10 = mesh[xu + m_phi * yd][0];
// 	f01 = mesh[xd + m_phi * yu][0];
// 	f11 = mesh[xu + m_phi * yu][0];
// 
// 	*xi = f00 + (f10 - f00) * ksi + (f01 - f00) * eta
// 			+ (f11 + f00 - f10 - f01) * ksi * eta;
// 
// 	f00 = mesh[xd + m_phi * yd][1];
// 	f10 = mesh[xu + m_phi * yd][1];
// 	f01 = mesh[xd + m_phi * yu][1];
// 	f11 = mesh[xu + m_phi * yu][1];
// 
// 	*yi = f00 + (f10 - f00) * ksi + (f01 - f00) * eta
// 			+ (f11 + f00 - f10 - f01) * ksi * eta;
// 
// 	f00 = mesh[xd + m_phi * yd][2];
// 	f10 = mesh[xu + m_phi * yd][2];
// 	f01 = mesh[xd + m_phi * yu][2];
// 	f11 = mesh[xu + m_phi * yu][2];
// 
// 	*zi = f00 + (f10 - f00) * ksi + (f01 - f00) * eta
// 			+ (f11 + f00 - f10 - f01) * ksi * eta;
// }
// 
// int SphericalHarmonicMeshSource::mod(int x, int y)
// {
// 	if( (x % y) >= 0 )
// 	{
// 		return x % y;
// 	}
// 	else
// 	{
// 		return y + (x % y);
// 	}
// }

} // end namespace neurolib
