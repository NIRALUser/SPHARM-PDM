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
#include "vtkPolyDataToitkMesh.h"

#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>


using namespace std;
void WriteEulerFile( std::string outEulerName, int Eulernum);
vtkSmartPointer<vtkPolyData> ReadPolyData(std::string filePath);

// static int debug = 0;

int main( int argc, char * argv[] )
{

  PARSE_ARGS;
  typedef itk::Image<unsigned short, 3>                             ImageType;
  typedef itk::ImageFileReader<ImageType>                           VolumeReaderType;
  typedef itk::BinaryMask3DEqualAreaParametricMeshSource<ImageType> MeshSourceType;
  typedef MeshSourceType::OutputMeshType                            OutputMeshType;
  typedef itk::MeshSpatialObject<OutputMeshType>                    MeshSpatialObjectType;
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



  vtkSmartPointer<vtkPolyData>      vtkmesh = vtkSmartPointer<vtkPolyData>::New();
  vtkPolyDataToitkMesh VTKITKConverter;

  // if necessary read parmesh
  typedef itk::SpatialObjectReader<3, double, OutputMeshType::MeshTraits> ReaderType;
  ReaderType::Pointer readerSH = ReaderType::New();
  if( initParaFile )
    {
      if (strstr(initParaFileName.c_str(),".meta"))
	{
	try
	  {
	    if( debug )
	      {
	      std::cout << "Reading " << initParaFileName << std::endl;
	      }
	    readerSH->SetFileName(initParaFileName);
	    readerSH->Update();
      ReaderType::GroupType::Pointer group1 = readerSH->GetGroup();
	    ReaderType::GroupType::ObjectListType * childList = group1->GetChildren(1, NULL);
	    // TODO: plugin name if multiple object are present
	    ReaderType::GroupType::ObjectListType::iterator it = childList->begin();
	    itk::SpatialObject<3> *                         curObj = *it;
	    MeshSpatialObjectType::Pointer                  paraSOMesh = dynamic_cast<MeshSpatialObjectType *>(curObj);
	    parmesh = paraSOMesh->GetMesh();
	    }
	catch( itk::ExceptionObject ex )
	    {
	    std::cout << ex.GetDescription() << std::endl;
	    return EXIT_FAILURE ;
	    }
        } 


      if (strstr(initParaFileName.c_str(),".vtk"))
        {
	  // convert surfaces to meta
	  try
	    {
	      if( debug )
		{
		  std::cout << "Reading vtk " << initParaFileName << std::endl;
		}
	      vtkmesh = ReadPolyData( initParaFileName );
	      
	      VTKITKConverter.SetInput( vtkmesh);
	      parmesh = VTKITKConverter.GetOutput();
	      std::cout << "Converting vtk done"  << std::endl;
	    }
	  catch( itk::ExceptionObject ex )
	    {
	      std::cout << ex.GetDescription() << std::endl;
	      return EXIT_FAILURE ;
	    }
	}
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

    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();

      {
      // convert surfaces to VTK
      itkMeshTovtkPolyData ITKVTKConverter;
      ITKVTKConverter.SetInput( mesh);
      vtkSmartPointer<vtkPolyData> SurfMesh = vtkSmartPointer<vtkPolyData>::New();
      SurfMesh = ITKVTKConverter.GetOutput();

      // Writing the file
      writer->SetInputData(SurfMesh);
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
      itkMeshTovtkPolyData ITKVTKConverter2;
      ITKVTKConverter2.SetInput( parmesh );
      vtkSmartPointer<vtkPolyData> ParaMesh = vtkSmartPointer<vtkPolyData>::New();
      ParaMesh = ITKVTKConverter2.GetOutput();
      // delete (ITKVTKConverter2);

      // Writing the file
      writer->SetInputData(ParaMesh);
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
