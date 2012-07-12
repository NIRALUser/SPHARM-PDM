#ifndef _namicSphericalHarmonicMedialAxisMeshSource_h
#define _namicSphericalHarmonicMedialAxisMeshSource_h

/** \class SphericalHarmonicMedialAxisMeshSource
 *
 *  \brief This class provides convenient class for generating 3D meshes represented by Spherical Harmonics descriptor and the medial axis associated.
 *
 *  \author Jean-Baptiste Berger
 */

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>

#include "SphericalHarmonicPolynomial.h"
#include "SphericalHarmonicSpatialObject.h"

#include "itkMeshSource.h"
#include "itkMesh.h"

namespace neurolib
{

	class SphericalHarmonicMedialAxisMeshSourceException : public itk::ExceptionObject
	{
		public:
			/** Run-time information. */
			itkTypeMacro( SphericalHarmonicMedialAxisMeshSourceException, ExceptionObject );

			/** Constructor. */
			SphericalHarmonicMedialAxisMeshSourceException(const char *file, unsigned int line,
                                       const char* message = "Error in generating meshes from Spherical Harmonics") :
					ExceptionObject(file, line)
					{
						SetDescription(message);
					}

					/** Constructor. */
					SphericalHarmonicMedialAxisMeshSourceException(const std::string & file, unsigned int line,
                                       const char* message = "Error in generating meshes from Spherical Harmonics") :
							ExceptionObject(file, line)
							{
								SetDescription(message);
							}

	};

	typedef SphericalHarmonicSpatialObject::MeshType SphericalHarmonicMedialAxisMeshSourceMeshType;
// default dimension for mesh is 3
	class SphericalHarmonicMedialAxisMeshSource : public itk::MeshSource<SphericalHarmonicMedialAxisMeshSourceMeshType>
	{
		public:
			typedef SphericalHarmonicMedialAxisMeshSource         Self;
			typedef SphericalHarmonicMedialAxisMeshSourceMeshType MeshType;
			typedef itk::MeshSource<MeshType>                     Superclass;
			typedef itk::SmartPointer<Self>                       Pointer;
			typedef itk::SmartPointer<const Self>                 ConstPointer;

			/** Method for creation through the object factory. */
			itkNewMacro(Self);

			/** Run-time type information (and related methods). */
			itkTypeMacro(SphericalHarmonicMedialAxisMeshSource, MeshSource);

			typedef SphericalHarmonicSpatialObject::ScalarType   ScalarType;
			typedef SphericalHarmonicSpatialObject::CoefType     CoefType;
			typedef SphericalHarmonicSpatialObject::CoefListType CoefListType;

			typedef double Point3[3];

			void SetCoefs(CoefListType& coeflist);

			typedef SphericalHarmonicSpatialObject::MeshTrait MeshTrait;
			typedef MeshType::PointsContainerPointer          PointsContainerPointer;
			typedef MeshType::PointsContainer                 PointsContainer;
			typedef MeshType::PointType                       PointType;
			typedef MeshType::CellsContainerPointer           CellsContainerPointer;
			typedef MeshType::CellsContainer                  CellsContainer;
			typedef MeshType::CellType                        CellType;

// 			itkGetConstReferenceMacro(Level, unsigned int);
// 			itkSetMacro(Level, unsigned int);

			itkGetConstReferenceMacro(Degree, unsigned int);
			void SetDegree(unsigned int d);

			itkGetConstReferenceMacro(FromL, unsigned int);
			itkSetMacro(FromL, unsigned int);

			itkGetConstReferenceMacro(ToL, unsigned int);
			itkSetMacro(ToL, unsigned int);

			itkGetConstReferenceMacro(Dimension, unsigned int);

			void SetThetaPhiIteration(int theta, int phi);
			
			double* GetTheta(){return m_ThetaTable;}
			double* GetRadius(){return m_radius;}
			
			vtkSmartPointer<vtkPolyData> GetOutputMedialAxis()
			{
				return m_outputMedialAxis;
			}

		protected:
			SphericalHarmonicMedialAxisMeshSource();
			~SphericalHarmonicMedialAxisMeshSource();
			
			void GenerateData();
			
			void setPhiThetaTable();
			void GetBounds(double[],double[],int);
// 
// 			void set_up_icosahedron_triangs(Point3* all_vert, Point3* all_triangs, int subdiv, int n_vert, int n_phi, int n_theta,
// 													  double *icos, int n_triangs,
// 													  int *triangs);
// 
// 			void interpol_vert(int n_phi, int n_theta, Point3 *mesh, int n_vertex, double *icos, Point3 *vertex);
// 
// 			int mod(int x, int y);
// 
// 			void interpol_2d(int n_phi, int n_theta, int xd, int xu, int yd, int yu, double ksi, double eta, Point3 *mesh,
// 								  double *xi, double *yi,
// 								  double *zi);

		private:
// 			MeshType::Pointer outputParaMesh;
// 			MeshType::Pointer m_outputMedialAxis;
			vtkSmartPointer<vtkPolyData> m_outputMedialAxis;

			unsigned int m_Dimension; // dimension of the output mesh

			unsigned int m_Degree; // degree of the coefficients
			CoefListType m_Coefs;
			int m_theta;
			int m_phi;

			unsigned int m_FromL;
			unsigned int m_ToL;
			
			double *m_PhiTable;
			double *m_ThetaTable;
			double *m_ThetaPhiTable;
			double *m_radius;

// 			double *m_icos;   // icosahedron
	};

}

#endif
