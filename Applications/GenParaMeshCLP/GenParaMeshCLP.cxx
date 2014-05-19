/*
 * compute the spherical parametrization of a binary mask
 *
 * author:  Martin Styner
 *
 */

#include <fstream>
#include <string>
#include <string.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkMeshSpatialObject.h>
#include <itkMesh.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkImageFileReader.h>

#include "BinaryMask3DEqualAreaParametricMeshSource.h"

#include "GenParaMeshCLPCLP.h"
#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataWriter.h"

using namespace std;
void WriteEulerFile( std::string outEulerName, int Eulernum);

// static int debug = 0;

int main( int argc, char * argv[] )
{

  PARSE_ARGS;
  typedef itk::Image<unsigned short, 3>                             ImageType;
  typedef itk::ImageFileReader<ImageType>                           VolumeReaderType;
  typedef itk::BinaryMask3DEqualAreaParametricMeshSource<ImageType> MeshSourceType;
  typedef MeshSourceType::OutputMeshType                            OutputMeshType;
  typedef itk::MeshSpatialObject<OutputMeshType>                    MeshSpatialObjectType;
  //   typedef itk::SpatialObjectWriter<3,float,OutputMeshType::MeshTraits> MeshWriterType;
  typedef itk::ImageFileReader<ImageType> VolumeReaderType;

  MeshSourceType::Pointer meshsrc = MeshSourceType::New();

  ImageType::Pointer image;
  bool               initParaFile;

  if( initParaFileName == "NULL" )
    {
    initParaFile = false;
    }
  else
    {
    initParaFile = true;
    }

  if( debug )
    {
    std::cout << "Reading Image: " << std::endl;
    }
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(infile);
  imageReader->Update();
  image = imageReader->GetOutput();
  OutputMeshType::Pointer mesh;
  OutputMeshType::Pointer parmesh;

  ofstream log;

  if( logFile )
    {
    log.open(outLogName.c_str(), ios::out | ios::app);
    }

  // if necessary read parmesh
  typedef itk::SpatialObjectReader<3, float, OutputMeshType::MeshTraits> ReaderType;
  ReaderType::Pointer readerSH = ReaderType::New();
  if( initParaFile )
    {
    try
      {
      readerSH->SetFileName(initParaFileName);
      readerSH->Update();
      }
    catch( itk::ExceptionObject ex )
      {
      std::cout << ex.GetDescription() << std::endl;
      return EXIT_FAILURE ;
      }
    ReaderType::SceneType::Pointer          scene1 = readerSH->GetScene();
    ReaderType::SceneType::ObjectListType * objList =  scene1->GetObjects(1, NULL);
    // TODO: plugin name if multiple object are present
    ReaderType::SceneType::ObjectListType::iterator it = objList->begin();
    itk::SpatialObject<3> *                         curObj = *it;
    MeshSpatialObjectType::Pointer                  paraSOMesh = dynamic_cast<MeshSpatialObjectType *>(curObj);
    parmesh = paraSOMesh->GetMesh();
    }

  try
    {
    if( debug )
      {
      std::cout << "Creating Para Surface Mesh: " << std::endl;
      }
    //MeshSourceType::Pointer meshsrc = MeshSourceType::New();
    meshsrc->SetInput(image);
    meshsrc->SetNumberOfIterations(numIterations);
    meshsrc->SetObjectValue(label);
    if( initParaFile )
      {
      meshsrc->SetInitParametricMesh(parmesh);
      }
    meshsrc->Update();
    // Output Mesh
    mesh = meshsrc->GetSurfaceMesh();
    // Create the mesh Spatial Object

   if (meshsrc->GetEulerNum() == 2 )
   {

    vtkPolyDataWriter *writer;
    writer = vtkPolyDataWriter::New();

      {
      // convert surfaces to VTK
      itkMeshTovtkPolyData * ITKVTKConverter = new itkMeshTovtkPolyData;
      ITKVTKConverter->SetInput( mesh);
      vtkPolyData *SurfMesh;
      SurfMesh = ITKVTKConverter->GetOutput();

      // Writing the file
      #if VTK_MAJOR_VERSION > 5
      writer->SetInputData(SurfMesh);
      #else
      writer->SetInput(SurfMesh);
      #endif
      writer->SetFileName(outSurfName.c_str() );

      writer->Write();

      }

    if( EulerFile )
      {
      WriteEulerFile(outEulerName, meshsrc->GetEulerNum() );
      }

    // Output Mesh
    parmesh = meshsrc->GetParametrizationMesh();
    // Create the mesh Spatial Object
      {
      // convert surfaces to VTK
      itkMeshTovtkPolyData * ITKVTKConverter2 = new itkMeshTovtkPolyData;
      ITKVTKConverter2->SetInput( parmesh );
      vtkPolyData *ParaMesh;
      ParaMesh = ITKVTKConverter2->GetOutput();
      // delete (ITKVTKConverter2);

      // Writing the file
      #if VTK_MAJOR_VERSION > 5
      writer->SetInputData(ParaMesh);
      #else
      writer->SetInput(ParaMesh);
      #endif
      writer->SetFileName(outParaName.c_str() );
      writer->Write();
      if( logFile )
        {
        log << "Computed " << infile << std::endl;
        }
      }
   }

    }
  catch( itk::ExceptionObject e )
    {
    if( EulerFile )
      {
      WriteEulerFile(outEulerName, meshsrc->GetEulerNum() );
      }

    e.Print(std::cout);
    if( logFile )
      {
      log << "Failed " << infile << " " << e.what();
      }
    return EXIT_FAILURE ;
    }
  return EXIT_SUCCESS ;
}

void WriteEulerFile( std::string outEulerName, int Eulernum)
{
  std::ofstream file( outEulerName.c_str() );

  if( file )
    {
    file << Eulernum << std::endl;
    file.close();
    }

}
