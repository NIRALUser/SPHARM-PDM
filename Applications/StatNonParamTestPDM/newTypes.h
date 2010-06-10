#ifndef NEWTYPES_H
#define NEWTYPES_H

#include "constants.h"

#include <itkMeshSpatialObject.h>
#include <itkMesh.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkDefaultDynamicMeshTraits.h>

typedef itk::DefaultDynamicMeshTraits< double , dimension, dimension, double, double, double > MeshTrait;
typedef itk::Mesh<double,dimension, MeshTrait> MeshType;

/** Hold on to the type information specified by the template parameters. */
typedef  MeshType::Pointer              MeshPointer;
typedef  MeshTrait::PointType           MeshPointType;
typedef  MeshTrait::PixelType           MeshPixelType; 
typedef  MeshType::Pointer              MeshPointer;
typedef  MeshType::CellTraits           CellTraits;
typedef  MeshType::PointsContainerPointer PointsContainerPointer;
typedef  MeshType::PointsContainer      PointsContainer;
typedef  MeshType::CellsContainerPointer CellsContainerPointer;
typedef  MeshType::CellsContainer       CellsContainer;
typedef  MeshType::PointType            PointType;
typedef  MeshType::CellType             CellType;
typedef  itk::TriangleCell<CellType>   TriangleType;

typedef itk::MeshSpatialObject<MeshType> MeshSpatialObjectType;
typedef itk::SpatialObjectWriter<dimension,double,MeshTrait> MeshWriterType;
typedef itk::SpatialObjectReader<dimension,double,MeshTrait> MeshReaderType;

typedef itk::DefaultDynamicMeshTraits < double, dimension, dimension, double, double > MeshTraitsType ; 
typedef itk::Mesh < double, dimension, MeshTraitsType > itkMeshType ;
typedef itk::MeshSpatialObject < itkMeshType > itkMeshSOType ;
typedef itk::MetaMeshConverter < dimension, double, MeshTraitsType > MeshConverterType ;

#endif
