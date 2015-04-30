#include "argio.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataToitkMesh.h"
#include "itkMesh3DProcrustesAlignFilter.h"

#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkDataSet.h"
#include "vtkPointData.h"
// bp2009
#include "itkSimplexMesh.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkFeatureEdges.h"
#include "vtkStripper.h"
#include "vtkTriangleFilter.h"
#include "vtkDelaunay3D.h"
#include "vtkDelaunay2D.h"
#include "vtkAppendPolyData.h"
#include "itkRegularSphereMeshSource.h"
#include "itkBinaryMask3DMeshSource.h"
#include "itkMeshSource.h"
#include "itkSimplexMeshToTriangleMeshFilter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkFloatArray.h"
#include "vtkColorTransferFunction.h"
#include "vtkLookupTable.h"
#include <vtkCurvatures.h>
#include <vtkSmartPointer.h>
#include "vtkDoubleArray.h"
#include <vtkCell.h>
#include <vtkPoints.h>
#include <vtkVersion.h>

#include <vtkTable.h>
#include <vtkDelimitedTextReader.h>

#if !defined(M_PI)
#define M_PI 3.14159265358979323846264338327950288   /* pi */
#endif
// bp2009
// bp2010
#include "vtkPointLocator.h"
// bp2010

#include <itkDefaultDynamicMeshTraits.h>
#include <itkMetaMeshConverter.h>
#include <itkVariableLengthVector.h>
#include <itkMetaArrayReader.h>
#include <itkMetaArrayWriter.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkSceneSpatialObject.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkContinuousIndex.h>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include <itkImageFileWriter.h>
#include <itkDanielssonDistanceMapImageFilter.h>
#include "itkBinaryThresholdImageFilter.h"
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkPlane.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDecimatePro.h>
#include <sstream>

//#include "MeshMathImpl.cxx"

// itk typedefs
typedef itk::DefaultDynamicMeshTraits<float, 3, 3, float, float> MeshTraitsType;
typedef itk::Mesh<float, 3, MeshTraitsType>                      MeshType;
typedef itk::MeshSpatialObject<MeshType>                         MeshSOType;
typedef itk::MetaMeshConverter<3, float, MeshTraitsType>         MeshConverterType;
typedef MeshTraitsType::PointType                                PointType;

typedef MeshTraitsType::CellType CellType;

typedef itk::SpatialObjectWriter<3, float, MeshTraitsType> MeshWriterType;
typedef itk::SceneSpatialObject<3>                         SceneSpatialObjectType;

typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> TriangleMeshTraits;
typedef itk::Mesh<double, 3, TriangleMeshTraits>                    TriangleMeshType;
typedef TriangleMeshTraits::PointType                               PointTriangleType;
typedef itk::MetaMeshConverter<3, double, TriangleMeshTraits>       TriangleMeshConverterType;
typedef itk::MeshSpatialObject<TriangleMeshType>                    TriangleMeshSOType;
typedef itk::SpatialObjectWriter<3, double, TriangleMeshTraits>     TriangleMeshWriterType;

typedef float                                ValueType;
typedef itk::VariableLengthVector<ValueType> MeasurementVectorType;
typedef itk::VariableLengthVector<float>     VectorType;
typedef itk::MetaArrayReader                 VectorReaderType;
typedef itk::MetaArrayWriter                 VectorWriterType;

typedef itk::Mesh3DProcrustesAlignFilter<MeshType, MeshType> ProcrustesFilterType;
typedef itk::SpatialObjectReader<3, float, MeshTraitsType>   ReaderType;

typedef MeshType::CellType                      CellType;
typedef itk::LineCell<CellType>                 LineType;
typedef itk::TriangleCell<CellType>             TriangleType;
typedef MeshType::CellsContainer::ConstIterator CellIterator;
typedef  MeshType::Pointer                      MeshPointer;

typedef float                                                           PixelType;
typedef itk::Image<PixelType, 3>                                        ImageType;
typedef itk::LinearInterpolateImageFunction<ImageType, double>          LinearInterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NearestNeighborInterpolatorType;
typedef itk::ImageFileReader<ImageType>                                 ImageReaderType;
typedef itk::ImageFileWriter<ImageType>                                 ImageWriterType;

typedef itk::Image<short, 3>                                                  ImageTypeShort;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageTypeShort>            Threshold;
typedef itk::DanielssonDistanceMapImageFilter<ImageTypeShort, ImageTypeShort> Distance;

// bp2009
typedef double                                                                  PixelDoubleType;
typedef itk::Image<PixelDoubleType, 3>                                          ImageDoubleType;
typedef itk::DefaultStaticMeshTraits<double, 3, 3, double, double, double>      SimplexMeshTraits;
typedef itk::SimplexMesh<double, 3, SimplexMeshTraits>                          SimplexMeshType;
typedef itk::RegularSphereMeshSource<TriangleMeshType>                          SphereMeshSourceType;
typedef itk::MeshSource<TriangleMeshType>                                       MeshSourceType;
typedef SphereMeshSourceType::PointType                                         SpherePointType;
typedef itk::SimplexMeshToTriangleMeshFilter<SimplexMeshType, TriangleMeshType> SimplexToTriangleType;
typedef float                                                                   LUTValueType;
// bp2009

static int GetIntFromString(const char *inputLine)
{
  std::string Line = inputLine;
  int         loc = Line.find("=", 1);
  std::string NbPoint = Line.erase(0, loc + 1);
  int         Val = atoi(NbPoint.c_str() );
  return Val;
}
int roundd (double d)
{
 return static_cast<int>(d + 0.5);

}
/* qsort float comparision function */
static int oldStyleCompare( const void* a_, const void* b_ )
{
  const float a = *(float *)a_;
  const float b = *(float *)b_;

  if( a < b )
    {
    return -1;
    }
  if( a > b )
    {
    return 1;
    }
  return 0;
}

int compare(const void * a, const void * b)
{
  return *(int *)a - *(int *)b;
}

int main(int argc, const char * *argv)
{
  // make sure the arguments are valid
  if( argc < 3 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help") )
    {
    std::cout << "Usage:" << std::endl;
    std::cout << argv[0] << " inputmesh/inputarray OutputFileName [options]" << std::endl;
    std::cout << " -subtract <meshfile>      Subtract mesh from inputmesh, write a KWMeshVisu readable text file"  << std::endl;
    std::cout << " -magnitude                Magnitude of the input metaArray file (mvh/mva) and writes a KWMeshVisu readable file"    << std::endl;
    std::cout << " -scaleMVA <double>        Scales the input metaArray file (mvh/mva) and writes a KWMeshVisu readable file"  << std::endl;
    std::cout << " -scaleMesh <DimX> <DimY> <DimZ>       Scales the input mesh file" << std::endl;
    std::cout << " -avgMesh <Meshfile1> <Meshfile2> ..."  << std::endl;
    std::cout << "   Compute the average mesh from inputmesh file1, file2..." << std::endl;
    std::cout << " -ave <Vectorfile1> <Vectorfile2> ..."  << std::endl;
    std::cout << "   Compute the average vector field from file1, file2... generated with -substract"
          << std::endl;
    std::cout << " -normave  <Vectorfile1> <Vectorfile2> ..." << std::endl;
    std::cout << "   Works as the \"-ave\" option, but the average vector are projected on the normal at each point "    << std::endl;
    std::cout << " -InvVect <VectorFile>     Invert all the vectors created with the -substract option and write a KWMeshVisu readable file"   << std::endl;
    std::cout << " -magdir <VectorFile>  Compute the signed magnitude of each of the vector from the vector field.(+ if in the normal direction, - otherwise)"   << std::endl;
    std::cout << " -magNormDir <VectorFile> Compute the signed magnitude of the normal projection of the vector field"   << std::endl;
    std::cout << " -applyVec <VectorFile> Deforme the mesh according to the vector field specified as input" << std::endl;
    std::cout << " -meshValues                Find the points and cells in a mesh. The outputfile is a textfile with the values"<< std::endl;
    std::cout << " -avgGaussMesh <Meshfile1> <Meshfile2> ... -gaussMeshPara <mean>,<stdev>,<val1>,<val2>,... "<< std::endl;
    std::cout << "   Compute the gaussian average for mesh files." << std::endl;
    std::cout << "   The first parameter is the average, then the standard deviation of the Gaussian model and the rest are the values associated with the files"   << std::endl;
    std::cout << " -avgGaussKWM <txtfile1> <txtfile2>... -gaussKWMPara <mean>,<stdev>,<val1>,<val2>,... "
      << std::endl;
    std::cout << "  Compute the gaussian average for KWMeshVisu files." << std::endl;
    std::cout
    <<
    "  The first parameter is the average, then the standard deviation of the Gaussian model and the rest are the values associated with the files"
    << std::endl;
    std::cout
    <<
    " -alignMesh <Meshfile1> <Meshfile2>... Align all of the meshes to the inputmesh (== MeshFile0) using Procrustes alignment"
    << endl;
    std::cout << " -BadTriangle <thresh value> [-correctMesh correctFilename] "  << std::endl;
    std::cout
    <<
    "  Find the bad triangles in a Mesh. The <thresh value> is the value of the threshFactor to calculate the standard deviation for the bad triangles. The output is a KWMeshVisu text file with the values of the average of the triangles of the mesh "
    << std::endl;
    std::cout << "  To have a new Mesh with the correct triangles -correctMesh " <<  std::endl;
    std::cout << " -extraction extractFilename [-extractClosest] [-nn]" << std::endl;
    std::cout
    <<
    "  To extract an attribute.The Input is the Mesh, the extractFilename is the attribute image and the Output is a KWMeshVisu text file with the attribute extraction"
    << std::endl;
    std::cout << "  [-extractClosest]: extract closest attribute" << std::endl;
    std::cout << "  [-nn]: nearest neighbor interpolation (default: linear)" << std::endl;
    std::cout << " -value <file1> <file2>... " << std::endl;
    std::cout
    << "  Extract the 5th column from a textfile and write a KWMeshVisu file with the values obtained"
    <<  std::endl;
    std::cout << " -subKWM <textname>       Difference between 2 KWMeshVisu files" << std::endl;
    std::cout
    <<
    " -MaxColor <textfile>...  Compare each point in every files, find a max for every points, keep 5% near the max, the other values will be 0"
    << std::endl;
    std::cout
    <<
    " -dist_absolute <textfile>,<textfile>...  -result_absolute <textfile>,<textfile>... Absolute distance map between KWMeshVisu files"
    << std::endl;
    std::cout
    <<
    " -dist_relative <textfile>,<textfile>...  -result_relative <textfile>,<textfile>... Relative distance map between KWMeshVisu files (values between -1 & 1)"
    << std::endl;
    std::cout << " -label <textfile>      Separate every labels, find the mean..." << std::endl;
    std::cout
    << " -color -val <number_of_label>,<value_label>... -oldval <number_of_old_label>,<old_value_label>..."
    << std::endl;
    std::cout << "  To change the value of labels to see the evolution with KWMeshVisu. " <<  std::endl;
    std::cout << "  Value_label is when the label grow up. " <<  std::endl;
    std::cout << "  Old_value_label is for the label wich has already grown up." <<  std::endl;
    std::cout
    << " -first <textfile>...   Convert a column file into a line file with a comma between each value"
    << std::endl;
    // cchou MC2Origin
    std::cout << " -MC2Origin       Translate the Center of Mass to the Origin" << std::endl;
    // bp2009 StatsKWM
    std::cout << " -avgOneKWM       Computes the avg of an input KWMeshVisu readable file" << std::endl;
    std::cout << " -medianOneKWM    Computes the min of an input KWMeshVisu readable file" << std::endl;
    std::cout << " -minOneKWM       Computes the min of an input KWMeshVisu readable file" << std::endl;
    std::cout << " -maxOneKWM       Computes the max of an input KWMeshVisu readable file" << std::endl;
    std::cout << " -per1OneKWM      Computes the 1% percentile of an input KWMeshVisu readable file" << std::endl;
    std::cout << " -per99OneKWM     Computes the 99% percentile of an input KWMeshVisu readable file" << std::endl;
    // bp2009 StatsKWM
    // bp2009 FillHole
    std::cout << " -FillHole        Fills up a hole in a open mesh." << std::endl;
    std::cout << "      If more than one hole exists, this operation might have to be repeated."
      << std::endl;
    // bp2009 FillHole
    // bp2009 BordersOut
    std::cout << " -BordersOut      Outputs the borders of a mesh (if there)." << std::endl;
    // bp2009 BordersOut
    // bp2009 IsOpen
    std::cout << " -IsOpen          Gives back an integer defining whether the mesh is open or not" << std::endl;
    // bp2009 IsOpen
    // bp2009 CleanMesh
    std::cout
    << " -CleanMesh       Re-mesh the input mesh and gives back a new clean mesh without degenerated triangles"
    << std::endl;
    // bp2009 CleanMesh
    // bp2009 SmoothMesh
    std::cout << " -SmoothMesh iterations    Gives back a Laplacian smoothed surface" << std::endl;
    std::cout << "          Iterations defines how many times the Laplacian is applied."
      << std::endl;
    // bp2009 SmoothMesh
    // bp2009 FilterNormals
    std::cout
    <<
    " -FilterNormals direction <MeshFileVTK> <MeshFileVTKOut>  ... Changes homogeneously normals of the polygons in a mesh"
    << std::endl;
    std::cout
    <<
    "               direction= [1] normals outwards [-1] normals inwards"
    << std::endl;
    // bp2009 FilterNormals
    // bp2009 StatsROI
    std::cout << " -statsROI <txtROIFileIn>   Process a KWMeshVisu file, given a ROI Mask" << std::endl;
    std::cout << "       Outputs a new KWMeshVisu only with the info in the mask"
      << std::endl;
    // bp2009 StatsROI
    // bp2009 KWMtoPolyData
    std::cout
    <<
    " -KWMtoPolyData <txtFileIn> <nameScalarField>   Writes a KWM scalar field (1D) into a PolyData Field Data Scalar to visualize in Slicer3"
    << std::endl;
    // bp2009 KWMtoPolyData
    std::cout << " -significanceLevel <double> the min Pvalue for the Pval ColorMap " << std::endl;

    // bp2009 ProcessROI
    std::cout << " -processROI <txtROIFileIn> <MeshFileIn>  ... [TEMP - do not know where to put this]" << std::endl;
    std::cout << "   Gets stats for a distances ROI map" << std::endl;
    // bp2009 ProcessROI
    std::cout << " -surfaceArea <AttributeFile>   Computes surface area in a txt file" << std::endl;
    std::cout << " -lobarSurfaceArea <ParcellationAttributeFile>   Computes lobar surface area (output: csv file)"<<std::endl;
    std::cout << " -variance <AttributeFile2> <AttributeFile3>...   Compute variance across population" << std::endl;
    // bp2010 GetCurvatures
    std::cout
    <<
    " -GetCurvatures <txtFileOut_C> <txtFileOut_S> <txtFileOut_Gauss> <txtFileOut_mean>... Gets an assorted set of curvature measurements:"
    << std::endl;
    std::cout
    <<
    "   Koenderink curvature values for an input mesh (shape index = S, curvedness = C), Gaussian Curvature and Mean Curvature"
    << std::endl;
    // bp2010 GetCurvatures
    // bp2010 particleConsistency
    std::cout
    <<
    " -particleConsistency <vtkFileIn_1> <lptsFileIn_1> ... <vtkFileIn_n> <lptsFileIn_n> ... Generates new particle files where fliped particles does not appear"
    << std::endl;
    // bp2010 particleConsistency
    // bp2010 GetDirectionFeatures
    std::cout
    <<
    " -GetDirectionFeatures ... Generates three KWM loadable files that have directionality features with respect to X, Y and Z directions"
    << std::endl;
    // bp2010 GetDirectionFeatures
    std::cout
    <<
    " -closestPoint <InputAttributeFile1> <InputMesh2>   Computes interpolated attribute file (output file) for second mesh using closest point interpolation between two input meshes"
    << std::endl;
    std::cout << " -extractVertices <YLocationAttributeFile.txt> <ZLocationAttributeFile.txt>  Extract points and write 3 different files listing respectively X, Y and Z values"
    << std::endl;
    std::cout << " -mean <AttributeFile2> <AttributeFile3>... Compute mean scalar attribute file (assuming same number of vertices)" << std::endl;
    //cx2011 cart2bary
    std::cout << " -cart2bary <vtkPointFileIn> Project each point in <vtkPointFileIn> onto the <inputmesh> and then calculate the Barycentric coordinate of the projected point, output to <OutputFileName>" << std::endl;

    //bp2011
    std::cout << " -attSTD <attribute_file2> <attribute_file2> ...<attribute_fileN> ... Compute std for a collection of vector or scalar files" << std::endl;
    std::cout << " -VTKtolpts  ... Create a particle file out of a triangulated VTK mesh" << std::endl;
    //bp2012
    std::cout << " -relaxPolygons iterations... Mesh relaxation based in vtkWindowedSincPolyDataFilter" << std::endl;
    std::cout << " -decimateMesh <target reduction> ... Mesh decimation, reducing the number of points. Target reduction is a value from 0..1 that indicates in what % vertices should be reduced" << std::endl;
    //mst2013
    std::cout << " -listPointData    lists info about all the VTK point Data in the vtk file" << std::endl;
    std::cout << " -PointDataOp <name> <op> <val>   <name>: name of point data, " << std::endl
	      << "     <op>: [threshBelow | sub]  threshBelow: set to 0 all data below <value>, sub: subtract value"   << std::endl;
    //bp2013
    std::cout << " -translateMesh <tx_DimX> <tx_DimY> <tx_DimZ> ,  translates a mesh for a given amount " << std::endl;
    //bp2015
    std::cout << " -lookupPointData <csv_lookup> <scalar_field_name>,  substitutes scalar values in a given scalar field based on an input lookup " << std::endl;
    
    std::cout << " -verbose                   Verbose output" << std::endl;
    return 0;
    }

  // get the arguments
  char *      inputFilename = strdup(argv[1]);
  char *      outputFilename = strdup(argv[2]);
  std::string outputFilename2;
  int         DotLoc;

  bool colorOn =  ipExistsArgument(argv, "-color");
  // bool valOn=  ipExistsArgument(argv, "-val");
  float *num = new float[30];
  float *old_num = new float[30];
  // if(valOn)
  //  {
  //    char * tmp_str = ipGetStringArgument(argv, "-val", NULL);
  //  }
  // bool oldvalOn=  ipExistsArgument(argv, "-oldval");
  // if (oldvalOn)
  //  {
  //    char * tmp_str = ipGetStringArgument(argv, "-oldval", NULL);
  //  }

  const int maxNumFiles = 1000;
  int       nbfile = 0;
  char *    files[maxNumFiles];

  bool                     MaxColorOn = ipExistsArgument(argv, "-MaxColor");
  std::vector<std::string> MaxFiles;
  if( MaxColorOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-MaxColor", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      MaxFiles.push_back(files[i]);
      }
    }

  bool                     distAbsOn = ipExistsArgument(argv, "-dist_absolute");
  std::vector<std::string> distFiles;
  if( distAbsOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-dist_absolute", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      distFiles.push_back(files[i]);
      }
    }

  bool                     resultAbsOn = ipExistsArgument(argv, "-result_absolute");
  std::vector<std::string> resultFiles;
  if( resultAbsOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-result_absolute", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      resultFiles.push_back(files[i]);
      }
    }

  bool                     distRelOn = ipExistsArgument(argv, "-dist_relative");
  std::vector<std::string> distReFiles;
  if( distRelOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-dist_relative", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      distReFiles.push_back(files[i]);
      }
    }

  bool                     resultRelOn = ipExistsArgument(argv, "-result_relative");
  std::vector<std::string> restFiles;
  if( resultRelOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-result_relative", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      restFiles.push_back(files[i]);
      }
    }

  bool   subtractOn = ipExistsArgument(argv, "-subtract");
  char * subtractFile = ipGetStringArgument(argv, "-subtract", NULL);
  bool   magnitudeOn = ipExistsArgument(argv, "-magnitude");

  bool scaleMVAOn = ipExistsArgument(argv, "-scaleMVA");
  bool scaleOn = ipExistsArgument(argv, "-scaleMesh");
  // bool scale3DOn = ipExistsArgument(argv, "-scale3DMesh");
  double                   scaleFactor = 1.0;
  std::vector<std::string> vector_scaleFactor;
  int                      numscaleparams = nbfile + 3;
  float *                  ScaleParams = new float[numscaleparams];
  if( scaleMVAOn )
    {
    scaleFactor = ipGetDoubleArgument(argv, "-scaleMVA", 1.0);
    }

  else if( scaleOn )
    {
    char * tmp_str = ipGetStringArgument(argv, "-scaleMesh", NULL);
    int    nbfactor = ipExtractFloatTokens(ScaleParams, tmp_str, numscaleparams);
    if( nbfactor == 1 )
      {
      ScaleParams[1] = ScaleParams[0];
      ScaleParams[2] = ScaleParams[0];
      }
    }

  bool aveOn = ipExistsArgument(argv, "-ave");
  bool normaveOn = ipExistsArgument(argv, "-normave");
  bool averageOn = aveOn || normaveOn;
  if( aveOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-ave", files, maxNumFiles);
    }
  if( normaveOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-normave", files, maxNumFiles);
    }
  std::vector<std::string> AveFiles;
  if( averageOn )
    {
    for( int i = 0; i < nbfile; i++ )
      {
      AveFiles.push_back(files[i]);
      }
    }

  bool                     meshValueOn =  ipExistsArgument(argv, "-meshValues");
  bool                     alignMeshOn =  ipExistsArgument(argv, "-alignMesh");
  std::vector<std::string> AlignMeshFiles;
  std::vector<std::string> AlignMeshOutputFiles;
  if( alignMeshOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-alignMesh", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      AlignMeshFiles.push_back(files[i]);
      }
    AlignMeshOutputFiles.assign(AlignMeshFiles.begin(), AlignMeshFiles.end() );
    std::vector<std::string>::iterator inputFile = AlignMeshOutputFiles.begin();
    AlignMeshOutputFiles.insert(inputFile, inputFilename );
    for( int i = 0; i < nbfile + 1; i++ )
      {
      AlignMeshOutputFiles[i].erase( AlignMeshOutputFiles[i].size() - 5);
      AlignMeshOutputFiles[i].insert(AlignMeshOutputFiles[i].size(), "_align.meta");
      }
    }

  bool valueOn = ipExistsArgument(argv, "-value");
  // /////////////uncomment if you want to select with the first column (like only one label)
  // double Index = ipGetDoubleArgument(argv, "-value",5.0);

  bool   labelOn = ipExistsArgument(argv, "-label");
  char * labelFilename = ipGetStringArgument(argv, "-label", NULL);

  bool   subKWMOn = ipExistsArgument(argv, "-subKWM");
  char * subKWMFilename = ipGetStringArgument(argv, "-subKWM", NULL);

  bool   badtriangleMeshOn = ipExistsArgument(argv, "-BadTriangle");
  double threshFactor = ipGetDoubleArgument(argv, "-BadTriangle", 5.0);
  char * correctFilename = ipGetStringArgument(argv, "-BadTriangle", NULL);
  bool   correctbadtriangleMeshOn = ipExistsArgument(argv, "-correctMesh");

  bool   extractionOn = ipExistsArgument(argv, "-extraction");
  char * extractFilename = ipGetStringArgument(argv, "-extraction", NULL);
  bool   extractClosestOn = ipExistsArgument(argv, "-extractClosest");
  bool   nearestNeighborOn = ipExistsArgument(argv, "-nn");

  bool                     avgMeshOn = ipExistsArgument(argv, "-avgMesh");
  std::vector<std::string> AvgMeshFiles;
  if( avgMeshOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-avgMesh", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      AvgMeshFiles.push_back(files[i]);
      }
    }

  bool                     firstOn = ipExistsArgument(argv, "-first");
  std::vector<std::string> firstfiles;
  if( firstOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-first", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      firstfiles.push_back(files[i]);
      }
    }

  bool                     avgGaussMeshOn = ipExistsArgument(argv, "-avgGaussMesh");
  std::vector<std::string> AvgGaussMeshFiles;
  if( avgGaussMeshOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-avgGaussMesh", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      AvgGaussMeshFiles.push_back(files[i]);
      }
    }

  bool   gaussMeshParaOn =  ipExistsArgument(argv, "-gaussMeshPara");
  int    numparams = nbfile + 1 + 2;
  float *GaussParams = new float[numparams];
  float *Age = new float[nbfile + 1];
  if( gaussMeshParaOn )
    {
    char * tmp_str = ipGetStringArgument(argv, "-gaussMeshPara", NULL);
    int    nbage = ipExtractFloatTokens(GaussParams, tmp_str, numparams);
    if( numparams != nbage )
      {
      std::cerr << "we need " << numparams << " comma separated entries" << std::endl;
      exit(1);
      }
    else
      {
      nbage = nbage - 2;
      for( int j = 0; j < nbage; j++ )
        {
        Age[j] = GaussParams[j + 2];
        }
      }
    }

  bool                     avgGaussKWMOn =  ipExistsArgument(argv, "-avgGaussKWM");
  std::vector<std::string> avgGaussKWMFiles;
  if( avgGaussKWMOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-avgGaussKWM", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      avgGaussKWMFiles.push_back(files[i]);
      }
    }

  bool   gaussKWMParaOn =  ipExistsArgument(argv, "-gaussKWMPara");
  int    numpara = nbfile + 1 + 2;
  float *GaussKWMParams = new float[numpara];
  float *AgeKWM = new float[nbfile + 1];
  if( gaussKWMParaOn )
    {
    char * tmp_str = ipGetStringArgument(argv, "-gaussKWMPara", NULL);
    int    nbage = ipExtractFloatTokens(GaussKWMParams, tmp_str, numpara);
    if( numpara != nbage )
      {
      std::cerr << "we need " << numpara << " comma separated entries" << std::endl;
      exit(1);
      }
    else
      {
      nbage = nbage - 2;
      for( int j = 0; j < nbage; j++ )
        {
        AgeKWM[j] = GaussKWMParams[j + 2];
        }
      }

    }

  char *VectFile = NULL;
  bool  invvectOn = ipExistsArgument(argv, "-InvVect");
  if( invvectOn )
    {
    VectFile = ipGetStringArgument(argv, "-InvVect", NULL);
    }

  bool magdirOn = ipExistsArgument(argv, "-magdir");
  if( magdirOn )
    {
    VectFile = ipGetStringArgument(argv, "-magdir", NULL);
    }

  bool magNormDirOn = ipExistsArgument(argv, "-magNormDir");
  if( magNormDirOn )
    {
    VectFile = ipGetStringArgument(argv, "-magNormDir", NULL);
    }

  bool applyvecOn =  ipExistsArgument(argv, "-applyVec");
  if( applyvecOn )
    {
    VectFile = ipGetStringArgument(argv, "-applyVec", NULL);
    }

  // cchou: MC2OriginOn
  bool MC2OriginOn = ipExistsArgument(argv, "-MC2Origin");

  bool debug = ipExistsArgument(argv, "-verbose");

  // bp2009 testFill
  bool testFillOn = ipExistsArgument(argv, "-FillHole");
  bool IsOpenOn = ipExistsArgument(argv, "-IsOpen");
  bool BordersOutOn = ipExistsArgument(argv, "-BordersOut");
  // bp2009 testFill
  // bp2009 StatsKWM
  bool avgOneKWMOn = ipExistsArgument(argv, "-avgOneKWM");
  bool medianOneKWMOn = ipExistsArgument(argv, "-medianOneKWM");
  bool minOneKWMOn = ipExistsArgument(argv, "-minOneKWM");
  bool maxOneKWMOn = ipExistsArgument(argv, "-maxOneKWM");
  bool per1OneKWMOn = ipExistsArgument(argv, "-per1OneKWM");
  bool per99OneKWMOn = ipExistsArgument(argv, "-per99OneKWM");
  // bp2009 StatsKWM
  // bp2009 CleanMesh
  bool cleanMeshOn = ipExistsArgument(argv, "-CleanMesh");
  // bp2009 CleanMesh
  // bp2009 SmoothMesh
  bool smoothMeshOn = ipExistsArgument(argv, "-SmoothMesh");
  int  smoothIterationNb = ipGetIntArgument(argv, "-SmoothMesh", 1);
  // bp2009 SmoothMesh
  // bp2009 FilterNormals
  bool filterNormalsOn = ipExistsArgument(argv, "-FilterNormals");
  // bp2009 FilterNormals
  bool   statsROIOn = ipExistsArgument(argv, "-statsROI");
  char * ROIFileName = ipGetStringArgument(argv, "-statsROI", NULL);

  // bp2009 KWMtoPolyData
  bool PvalueColorMapOn = ipExistsArgument(argv, "-significanceLevel");
// double PvalueColorMapNb;
  std::string PvalueColorMapNb_string;
  double      PvalueColorMapNb(0.0);
  if( PvalueColorMapOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-significanceLevel", files, maxNumFiles);
    if( nbfile != 1 )
      {
      std::cerr << "Error: Incorrect number of arguments!" << std::endl;
      exit(1);
      }
    else
      {
      nbfile = nbfile + 1; PvalueColorMapNb_string.append(files[0]);
      // istringstream(PvalueColorMapNb_string) >> PvalueColorMapNb; //convert the string (the argument0 to a double
      PvalueColorMapNb = atof(PvalueColorMapNb_string.c_str() );

      }
    }
  //std::cout << PvalueColorMapNb << std::endl;

// bp2009 KWMtoPolyData
  bool KWMtoPolyDataOn = ipExistsArgument(argv, "-KWMtoPolyData");
  if( KWMtoPolyDataOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-KWMtoPolyData", files, maxNumFiles);
    if( nbfile != 2 )
      {
      std::cerr << "Error: Incorrect number of arguments!" << std::endl;
      exit(1);
      }
    }

  bool  surfaceAreaOn = ipExistsArgument(argv, "-surfaceArea");
  char *AttributeFileName = ipGetStringArgument(argv, "-surfaceArea", NULL);

  bool lobarSurfaceAreaOn = ipExistsArgument(argv, "-lobarSurfaceArea");
  AttributeFileName = ipGetStringArgument(argv, "-lobarSurfaceArea",NULL);

  std::vector<std::string> VarFiles;
  bool                     varianceOn = ipExistsArgument(argv, "-variance");
  if( varianceOn )
    {
    VarFiles.push_back(inputFilename);
    nbfile = ipGetStringMultipArgument(argv, "-variance", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      VarFiles.push_back(files[i]);
      }
    nbfile++;
    }

  // bp2010 ParticleConsistency
  bool                     particleOn = ipExistsArgument(argv, "-particleConsistency");
  std::vector<std::string> TestMeshFiles;
  if( particleOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-particleConsistency", files, maxNumFiles);
    for( int i = 0; i < nbfile; i++ )
      {
      TestMeshFiles.push_back(files[i]);
      }
    }
  // bp2010 ParticleConsistency
  // bp2009 GetCurvatures
  bool GetCurvaturesOn = ipExistsArgument(argv, "-GetCurvatures");
  // bp2009 GetCurvatures

  // bp2009 GetDirectionFeatures
  bool GetDirectionFeaturesOn = ipExistsArgument(argv, "-GetDirectionFeatures");
  // bp2009 GetDirectionFeatures

  std::vector<std::string> closestPointFiles;
  bool                     closestPointOn = ipExistsArgument(argv, "-closestPoint");
  if( closestPointOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-closestPoint", files, maxNumFiles);
    if( nbfile != 2 )
      {
      std::cerr << "Error: closestPoint option expects two files!\nPlease type:" << std::endl;
      std::cerr << argv[0] << "-help" << std::endl;
      exit(1);
      }
    else
      {
      closestPointFiles.push_back(files[0]);
      closestPointFiles.push_back(files[1]);
      }
    }

  std::vector<std::string> extractVerticesFiles;
  bool                     extractVerticesOn = ipExistsArgument(argv, "-extractVertices");
  if( extractVerticesOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-extractVertices", files, maxNumFiles);
    if( nbfile != 2 )
      {
      std::cerr << "Error: extractVertices option expects two files!\nPlease type:" << std::endl;
      std::cerr << argv[0] << "-help" << std::endl;
      exit(1);
      }
    else
      {
      extractVerticesFiles.push_back(files[0]);
      extractVerticesFiles.push_back(files[1]);
      }
    }

  std::vector<std::string> meanFiles;
  bool meanOn = ipExistsArgument(argv, "-mean");
  if (meanOn)
    {
      meanFiles.push_back(inputFilename);
      nbfile = ipGetStringMultipArgument(argv, "-mean", files, maxNumFiles);
      for(int i = 0 ; i < nbfile ; i++)
	meanFiles.push_back(files[i]);
      nbfile++;
    }

  //cx2011 cart2bary
  bool cart2baryOn = ipExistsArgument(argv, "-cart2bary");
  std::string vtkPointFile;
  if(cart2baryOn) 
    vtkPointFile = ipGetStringArgument(argv, "-cart2bary", NULL);
  //cx2011 cart2bary

  bool attSTDOn = ipExistsArgument(argv, "-attSTD");
  std::vector<std::string> VectorStdFiles;
  if (attSTDOn)
    {
      nbfile = ipGetStringMultipArgument(argv, "-attSTD", files, maxNumFiles); 
      for(int i = 0 ; i < nbfile ; i++) 
	VectorStdFiles.push_back(files[i]);
    }

  bool lptsOn = ipExistsArgument(argv, "-VTKtolpts");
  //bp2012 relaxPolygons
  bool relaxPolygonsOn = ipExistsArgument(argv, "-relaxPolygons");
  int relaxPolygonsNb = ipGetIntArgument(argv, "-relaxPolygons",1);
  //bp2012 relaxPolygons

 //bp2012 decimateMesh
  bool decimateMeshOn = ipExistsArgument(argv, "-decimateMesh");
  double decimateMeshTR = ipGetDoubleArgument(argv, "-decimateMesh",1.0);
  //bp2012 decimateMesh

  //ms2013  pointData
  bool listPointDataOn  = ipExistsArgument(argv,"-listPointData");
  //ms2013  pointData
  bool PointDataOpOn = ipExistsArgument(argv, "-PointDataOp");
  std::vector<std::string> pointDataOpCmd;
  if( PointDataOpOn )
    {
    nbfile = ipGetStringMultipArgument(argv, "-PointDataOp", files, 3);
    for( int i = 0; i < nbfile; i++ )
      {
      pointDataOpCmd.push_back(files[i]);
      }
    }

  //bp2013
  bool translateOn = ipExistsArgument(argv, "-translateMesh");
  std::vector<std::string> vector_translateFactor;
  int numtranslateparams = nbfile + 3;
  float *TranslateParams = new float[numscaleparams];
  if( translateOn )
    {
    char * tmp_str = ipGetStringArgument(argv, "-translateMesh", NULL);
    ipExtractFloatTokens(TranslateParams, tmp_str, numtranslateparams);
    }

  //bp2015
  bool lookupPointOn = ipExistsArgument(argv, "-lookupPointData");
  if( lookupPointOn )
      {
      nbfile = ipGetStringMultipArgument(argv, "-lookupPointData", files, maxNumFiles);
      if( nbfile != 2 )
        {
        std::cerr << "Error: Incorrect number of arguments!" << std::endl;
        exit(1);
        }
      }

  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  // PROCESSING STARTS HERE
  if( subtractOn )
    {
    // read in the input files
    MeshConverterType::Pointer converter = MeshConverterType::New();

    if( debug )
      {
      std::cout << "Reading input mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer inputMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer() );
    MeshType::Pointer   inputMesh = inputMeshSO->GetMesh();

    if( debug )
      {
      std::cout << "subtracting Mesh " << std::endl;
      }
    MeshSOType::Pointer sutractMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(subtractFile).GetPointer());
    MeshType::Pointer   subtractMesh = sutractMeshSO->GetMesh();

    // make sure the two meshes have the same number of verts
    if( inputMesh->GetNumberOfPoints() != subtractMesh->GetNumberOfPoints() )
      {
      std::cerr << "Meshes do not have same number of vertices: " << inputMesh->GetNumberOfPoints()
                << " vs " << subtractMesh->GetNumberOfPoints() << std::endl;
      exit(-1);
      }

    std::ofstream outfile;
    outfile.open( outputFilename );

    PointType *           point1 = new PointType;
    PointType *           point2 = new PointType;
    itk::Vector<float, 3> diff;

    // print the header
    outfile << "NUMBER_OF_POINTS=" << inputMesh->GetNumberOfPoints() << std::endl;
    outfile << "DIMENSION=" << 3 << std::endl;
    outfile << "TYPE=Vector" << std::endl;
    for( unsigned int i = 0; i < inputMesh->GetNumberOfPoints(); i++ )
      {
      if( inputMesh->GetPoint( i, point1 ) && subtractMesh->GetPoint(i, point2) )
        {
        diff = *point1 - *point2;
        outfile << diff[0] << " " << diff[1] << " " << diff[2] << std::endl;
        }
      else
        {
        return -1;
        }
      }
    outfile.close();

    outputFilename2 = std::string(outputFilename);
    DotLoc = outputFilename2.find(".txt", 0);
    if( DotLoc != -1 )
      {
      outputFilename2.insert(DotLoc, ".mag");
      }
    else
      {
      outputFilename2 = std::string(outputFilename) + std::string(".mag");
      }

    outfile.open( outputFilename2.c_str() );
    //  outputFilename => Filename
    // print the header
    outfile << "NUMBER_OF_POINTS=" << inputMesh->GetNumberOfPoints() << std::endl;
    outfile << "DIMENSION=" << 1 << std::endl;
    outfile << "TYPE=Scalar" << std::endl;
    for( unsigned int i = 0; i < inputMesh->GetNumberOfPoints(); i++ )
      {
      if( inputMesh->GetPoint( i, point1 ) && subtractMesh->GetPoint(i, point2) )
        {
        diff = *point1 - *point2;
        outfile << sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]) << std::endl;
        }
      else
        {
        return -1;
        }
      }
    outfile.close();

    delete ( point1 );
    delete ( point2 );
    }
  else if( scaleOn )
    {

    // know if meta or vtk file
    int    bemeta;
    char * pch;
    pch = strrchr(inputFilename, 'v');
    if( !pch || (strcmp("vtk", pch) != 0) )
      {
      bemeta = 1;
      }
    else
      {
      bemeta = 0;
      }

    if( bemeta == 1 )
      {
      MeshConverterType::Pointer converter = MeshConverterType::New();
      if( debug )
        {
        std::cout << "Reading input mesh " << inputFilename << std::endl;
        }
      MeshSOType::Pointer              inputMeshSO =
        dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
      MeshType::Pointer                inputMesh = inputMeshSO->GetMesh();
      MeshType::PointsContainerPointer points = inputMesh->GetPoints();

      MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
      for( unsigned int pointID = 0; pointID < points->Size(); pointID++ )
        {
        PointType curPoint =  points->GetElement(pointID);
        PointType vert;
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          vert[dim] = curPoint[dim] * ScaleParams[dim];
          }
        pointsTmp->InsertElement(pointID, vert);
        }

      inputMesh->SetPoints(pointsTmp);
      MeshWriterType::Pointer writer = MeshWriterType::New();
      writer->SetInput(inputMeshSO);
      writer->SetFileName(outputFilename);
      writer->Update();
      }

    else // =>vtk file
      {
      vtkPolyDataReader *meshin = vtkPolyDataReader::New();
      meshin->SetFileName(inputFilename);
      meshin->Update();
      double      x[3];
      vtkPoints * PointVTK = vtkPoints::New();
      PointVTK = meshin->GetOutput()->GetPoints();

      for( int PointId = 0; PointId < (meshin->GetOutput()->GetNumberOfPoints() ); PointId++ )
        {
        PointVTK->GetPoint(PointId, x);
        double vert[3];
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          vert[dim] = x[dim] * ScaleParams[dim];
          }
        PointVTK->SetPoint(PointId, vert);
        }
      vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
      #if VTK_MAJOR_VERSION > 5
      SurfaceWriter->SetInputData(meshin->GetOutput() );
      #else
      SurfaceWriter->SetInput(meshin->GetOutput() );
      #endif
      SurfaceWriter->SetFileName(outputFilename);
      SurfaceWriter->Update();
      }
    }

  else if( magnitudeOn )
    {
    if( debug )
      {
      std::cout << "Magnitude of Vector field" << std::endl;
      }
    MeasurementVectorType mv;
    // load the metaArray
    VectorReaderType::Pointer vectorReader = VectorReaderType::New();
    vectorReader->SetFileName(inputFilename);
    vectorReader->Update();
    vectorReader->GetOutput(MET_FLOAT, &mv, true);

    std::ofstream outfile;
    outfile.open( outputFilename );

    // print the header
    outfile << "NUMBER_OF_POINTS = " << mv.GetSize() / 3 << std::endl;
    outfile << "DIMENSION = " << 1 << std::endl;
    outfile << "TYPE=Scalar" << std::endl;
    for( unsigned int i = 0; i < mv.GetSize(); i += 3 )
      {
      double length = sqrt(mv[i] * mv[i] + mv[i + 1] * mv[i + 1] + mv[i + 2] * mv[i + 2]);
      if( scaleMVAOn )
        {
        outfile << length * scaleFactor  << std::endl;
        }
      else
        {
        outfile << length << std::endl;
        }
      }

    }
  else if( scaleMVAOn )
    {
    if( debug )
      {
      std::cout << "Scaling of Vector field" << std::endl;
      }
    MeasurementVectorType mv;
    // load the metaArray
    VectorReaderType::Pointer vectorReader = VectorReaderType::New();
    vectorReader->SetFileName(inputFilename);
    vectorReader->Update();
    vectorReader->GetOutput(MET_FLOAT, &mv, true);

    std::ofstream outfile;
    outfile.open( outputFilename );

    // print the header
    outfile << "NUMBER_OF_POINTS = " << mv.GetSize() / 3 << std::endl;
    outfile << "DIMENSION = " << 3 << std::endl;
    outfile << "TYPE = Vector" << std::endl;
    for( unsigned int i = 0; i < mv.GetSize(); i += 3 )
      {
      outfile << mv[i] * scaleFactor << " " << mv[i + 1] * scaleFactor << " "
              << mv[i + 2] * scaleFactor << std::endl;
      }

    }
  else if( meshValueOn )
    {

    // int num(0);
    // read input mesh
    MeshConverterType::Pointer converter = MeshConverterType::New();
    MeshSOType::Pointer              SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer                mesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer Points = mesh->GetPoints();
    int                              numberOfPoints = Points->Size();

    // open output file
    std::ofstream outfileNor;
    outfileNor.open( outputFilename );
    outfileNor << "Number of  Points = " << numberOfPoints << std::endl;
    // find points
    for( unsigned int pointID = 0; pointID < Points->Size(); pointID++ )
      {
      PointType curPoint =  Points->GetElement(pointID);
      outfileNor << "Points " << pointID << " " << curPoint << std::endl;
      }

    // find number&points of cells
    typedef CellType::PointIdIterator PointIdIterator;
    CellIterator cellIterator = mesh->GetCells()->Begin();
    CellIterator cellEnd      = mesh->GetCells()->End();
    PointType    curPoint;
    while( cellIterator != cellEnd )
      {
      CellType *                cell = cellIterator.Value();
      TriangleType *            line = dynamic_cast<TriangleType *>(cell);
      LineType::PointIdIterator pit = line->PointIdsBegin();
      int                       pIndex1, pIndex2, pIndex3;

      pIndex1 = *pit;
      ++pit;
      pIndex2 = *pit;
      ++pit;
      pIndex3 = *pit;

      outfileNor << "Cells  " << num << " " << pIndex1 << " " << pIndex2 << " " << pIndex3 << std::endl;

      ++cellIterator;
      ++num;
      }

    // close outputFile
    outfileNor.close();

    }
  else if( extractionOn )
    {

    PointType point;
    PointType pixel;

    // read input
    MeshConverterType::Pointer converter = MeshConverterType::New();
    MeshSOType::Pointer SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer                mesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer points = mesh->GetPoints();
    int                              numPoints = points->Size();
    double * *                       pointIndex = NULL;
    pointIndex = new double *[numPoints];
    for( int c = 0; c < numPoints; c++ )
      {
      pointIndex[c] = new double[3];
      for( int d = 0; d < 3; d++ )
        {
        pointIndex[c][d] = 0.0;
        }
      }
    std::vector<double> extractValue(numPoints,0.0);

    // load attribute image
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(extractFilename);
    reader->Update();

    ImageType::Pointer                 image = reader->GetOutput();
    Distance::VectorImageType::Pointer closestVectorImage;

    if( extractClosestOn )
      {
      Threshold::Pointer threshold = Threshold::New();
      threshold->SetInput(image);
      threshold->SetLowerThreshold(1);
      threshold->SetUpperThreshold(32000);

      Distance::Pointer distance = Distance::New();
      distance->UseImageSpacingOn();
      distance->SquaredDistanceOn();
      distance->SetInput(threshold->GetOutput() );
      distance->Update();
      closestVectorImage = distance->GetVectorDistanceMap();
      }

    // Origin & Spacing
    const ImageType::SpacingType spacing = image->GetSpacing();
    const ImageType::PointType   origin = image->GetOrigin();

    int Element[3] = {-1, -1, 1};

    if( !nearestNeighborOn )
      {
      LinearInterpolatorType::ContinuousIndexType point_index;
      LinearInterpolatorType::Pointer             interpolator = LinearInterpolatorType::New();
      interpolator->SetInputImage(image);
      // for all points in the mesh, extract image attribute
      for( int index = 0; index < numPoints; index++ )
        {
        PointType curPoint =  points->GetElement(index);
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          pixel[dim] = (Element[dim] * curPoint[dim] - origin[dim]) / spacing[dim];
          point_index[dim] = pixel[dim];
          pointIndex[index][dim] = point_index[dim];
          }
        extractValue[index] = interpolator->EvaluateAtContinuousIndex(point_index);

        if( extractValue[index] == 0 && extractClosestOn )
          {
          Distance::VectorImageType::PixelType vectorPixel;
          Distance::VectorImageType::IndexType vectorIndex;
          for( int indexDim = 0; indexDim < 3; indexDim++ )
            {
            vectorIndex[indexDim] =  (long int) roundd( point_index[indexDim] );
            }
          vectorPixel = closestVectorImage->GetPixel(vectorIndex);
          for( unsigned int dim = 0; dim < 3; dim++ )
            {
            point_index[dim] = vectorIndex[dim] + vectorPixel[dim];
            pointIndex[index][dim] = point_index[dim];
            }
          extractValue[index] = interpolator->EvaluateAtContinuousIndex(point_index);
          if( extractValue[index] == 0 )
            {
            std::cerr << "error at extracting closest point" << std::endl;
            exit(1);
            }
          }
        }
      }
    else
      {
      NearestNeighborInterpolatorType::ContinuousIndexType point_index;
      NearestNeighborInterpolatorType::Pointer             interpolator = NearestNeighborInterpolatorType::New();
      interpolator->SetInputImage(image);
      // for all points in the mesh, extract image attribute
      for( int index = 0; index < numPoints; index++ )
        {
        PointType curPoint =  points->GetElement(index);
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          pixel[dim] = (Element[dim] * curPoint[dim] - origin[dim]) / spacing[dim];
          point_index[dim] = pixel[dim];
          pointIndex[index][dim] = point_index[dim];
          }
        extractValue[index] = interpolator->EvaluateAtContinuousIndex(point_index);

        if( extractValue[index] == 0 && extractClosestOn )
          {
          Distance::VectorImageType::PixelType vectorPixel;
          Distance::VectorImageType::IndexType vectorIndex;
          for( int indexDim = 0; indexDim < 3; indexDim++ )
            {
            vectorIndex[indexDim] =  (long int) roundd( point_index[indexDim] );
            }
          vectorPixel = closestVectorImage->GetPixel(vectorIndex);
          for( unsigned int dim = 0; dim < 3; dim++ )
            {
            point_index[dim] = vectorIndex[dim] + vectorPixel[dim];
            pointIndex[index][dim] = point_index[dim];
            }
          extractValue[index] = interpolator->EvaluateAtContinuousIndex(point_index);
          if( extractValue[index] == 0 )
            {
            std::cerr << "error at extracting closest point" << std::endl;
            exit(1);
            }
          }
        }
      }

    // write ouput KWMeshVisu text file, attribute extraction
    std::ofstream outfileNor;
    outfileNor.open( outputFilename );
    outfileNor << "NUMBER_OF_POINTS=" << numPoints << std::endl;
    outfileNor << "DIMENSION=" << 1 << std::endl;
    outfileNor << "TYPE=Scalar" << std::endl;
    for( int i = 0; i < numPoints; i++ )
      {
      outfileNor << extractValue[i] << std::endl;
      }
    // close outputFile
    outfileNor.close();
    for( int c = 0; c < numPoints; c++ )
      {
      delete[] pointIndex[c];
      }
    delete[] pointIndex;

    }
  else if( avgMeshOn )
    {

    	if( debug )
    	{
		std::cout << "Reading mesh " << " " << argv[2] << std::endl;
    	}

	vtkPolyDataReader *meshin = vtkPolyDataReader::New();
	meshin->SetFileName(inputFilename);
	meshin->Update();
	vtkPolyData *polydata = meshin->GetOutput();
	vtkPoints *avgPoints = vtkPoints::New();
	avgPoints=meshin->GetOutput()->GetPoints();

	int numMeshes = AvgMeshFiles.size();
	for( int index = 0; index < numMeshes; index++ )
        {
	      try
	        {

	          if( debug )
	          {
		        std::cout << "Reading  mesh " << AvgMeshFiles[index] << std::endl;
		  }
		  vtkPolyDataReader *meshin2 = vtkPolyDataReader::New();
		  meshin2->SetFileName((char*)AvgMeshFiles[index].c_str());
		  meshin2->Update();

		  vtkPoints *meshPoints = vtkPoints::New();
		  meshPoints=meshin2->GetOutput()->GetPoints();
		  vtkPoints *tmpPoints = vtkPoints::New();
		  for( unsigned int pointID = 0; pointID < meshPoints->GetNumberOfPoints(); pointID++ )
                  {
		    double curPoint[3];
		    double avgPoint[3];
		    double tmpPoint[3];
                    meshPoints->GetPoint(pointID, curPoint);
 		    avgPoints->GetPoint(pointID, avgPoint);
	            for( unsigned int dim = 0; dim < 3; dim++ )
                    {
 		       tmpPoint[dim] = curPoint[dim] + avgPoint[dim];
		    }
		    tmpPoints->InsertPoint(pointID, tmpPoint);
		  }

		        avgPoints = tmpPoints;
       		}
	        catch( itk::ExceptionObject ex )
	        {
        		std::cerr << "Error reading meshfile:  " << AvgMeshFiles[index] << std::endl << "ITK error: "
                  << ex.GetDescription() << std::endl;
		        exit(-3);
        	}
		
	 }

	vtkPoints *tmpPoints = vtkPoints::New();
    	for( unsigned int pointID = 0; pointID < avgPoints->GetNumberOfPoints(); pointID++ )
  	{
	      double avgPoint[3];
	      double tmpPoint[3];
              avgPoints->GetPoint(pointID, avgPoint);
	      for( unsigned int dim = 0; dim < 3; dim++ )
	      {
	        tmpPoint[dim] = avgPoint[dim] / (numMeshes + 1);
              }
              tmpPoints->InsertPoint(pointID, tmpPoint);
      	}
	avgPoints = tmpPoints;

    polydata->SetPoints(avgPoints);
//    polydata->Update();//no more Update() function for vtkPolyData: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_Update
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(polydata);
    #else
    SurfaceWriter->SetInput(polydata);
    #endif
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();
    if (debug)
    {
       std::cout << "Writing new mesh " << outputFilename << std::endl;
    }



   /* MeshConverterType::Pointer converter = MeshConverterType::New();

    // read input
    if( debug )
      {
      std::cout << "Reading  mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer                surfaceMesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer avgPoints = surfaceMesh->GetPoints();

    int numMeshes = AvgMeshFiles.size();
    // for all meshes
    for( int index = 0; index < numMeshes; index++ )
      {
      try
        {
        if( debug )
          {
          std::cout << "Reading  mesh " << AvgMeshFiles[index] << std::endl;
          }
        MeshSOType::Pointer inputMeshSO =
          dynamic_cast<MeshSOType *>(converter->ReadMeta(AvgMeshFiles[index].c_str() ).GetPointer() );
        MeshType::Pointer                inputMesh = inputMeshSO->GetMesh();
        MeshType::PointsContainerPointer points = inputMesh->GetPoints();
        MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
        for( unsigned int pointID = 0; pointID < points->Size(); pointID++ )
          {
          PointType curPoint =  points->GetElement(pointID);
          PointType vert;
          PointType curPointAvg =  avgPoints->GetElement(pointID);
          for( unsigned int dim = 0; dim < 3; dim++ )
            {
            vert[dim] = curPoint[dim] + curPointAvg[dim];
            }
          pointsTmp->InsertElement(pointID, vert);
          }

        avgPoints = pointsTmp;
        }
      catch( itk::ExceptionObject ex )
        {
        std::cerr << "Error reading meshfile:  " << AvgMeshFiles[index] << std::endl << "ITK error: "
                  << ex.GetDescription() << std::endl;
        exit(-3);
        }
      }

    MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
    for( unsigned int pointID = 0; pointID < avgPoints->Size(); pointID++ )
      {
      PointType curPoint =  avgPoints->GetElement(pointID);
      PointType vert;
      for( unsigned int dim = 0; dim < 3; dim++ )
        {
        vert[dim] = curPoint[dim] / (numMeshes + 1);
        }
      pointsTmp->InsertElement(pointID, vert);
      }
    avgPoints = pointsTmp;

    surfaceMesh->SetPoints(avgPoints);
    SOMesh->SetMesh(surfaceMesh);
    MeshWriterType::Pointer writer = MeshWriterType::New();
    writer->SetInput(SOMesh);
    writer->SetFileName(outputFilename);
    writer->Update();
*/
    }
  else if( alignMeshOn )
    {

    MeshConverterType::Pointer converter = MeshConverterType::New();

    // read input
    MeshSOType::Pointer SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer() );
    MeshType::Pointer                surfaceMesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer Points = surfaceMesh->GetPoints();
    MeshType::Pointer                newPoints = MeshType::New();

    int numMeshes = AlignMeshFiles.size();

    // align every mesh to the first mesh
    ProcrustesFilterType::Pointer procrustesFilter = ProcrustesFilterType::New();
    procrustesFilter->SetNumberOfInputs(numMeshes + 1);
    procrustesFilter->SetUseInitialAverageOff();
    procrustesFilter->SetUseNormalizationOff();
    procrustesFilter->SetUseScalingOff();
    procrustesFilter->SetUseSingleIterationOn();
    procrustesFilter->SetAlignRotationOn();

    procrustesFilter->SetInput(0, surfaceMesh);
    // read meshes and use Mesh3DProcrustesAlignFilter to align
    for( int index = 0; index < numMeshes; index++ )
      {
      MeshSOType::Pointer inputMeshSO =
        dynamic_cast<MeshSOType *>(converter->ReadMeta(AlignMeshFiles[index].c_str() ).GetPointer() );
      MeshType::Pointer   inputMesh = inputMeshSO->GetMesh();
      procrustesFilter->SetInput(index + 1, inputMesh);
      }

    procrustesFilter->Update();

    // create the align MeshFile for the input
    MeshWriterType::Pointer writer = MeshWriterType::New();
    MeshType::Pointer       RegisteredMesh = procrustesFilter->GetOutput(0);
    SOMesh->SetMesh(RegisteredMesh);
    writer->SetInput(SOMesh);
    writer->SetFileName(AlignMeshOutputFiles[0].c_str() );
    writer->Update();
    // create align MeshFiles
    for( int index = 0; index < numMeshes; index++ )
      {
      MeshType::Pointer   _RegisteredMesh = procrustesFilter->GetOutput(index + 1);
      MeshSOType::Pointer inputMeshSO =
        dynamic_cast<MeshSOType *>(converter->ReadMeta(AlignMeshFiles[index].c_str() ).GetPointer() );
      inputMeshSO->SetMesh(_RegisteredMesh);
      writer->SetInput(inputMeshSO);
      writer->SetFileName(AlignMeshOutputFiles[index + 1].c_str() );
      writer->Update();
      }
    }
  else if( badtriangleMeshOn )
    {

    // read input Mesh
    MeshConverterType::Pointer converter = MeshConverterType::New();
    MeshSOType::Pointer              SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer                mesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer points = mesh->GetPoints();
    int                              numPoints = points->Size();
    double *                         avgTriangleArea = new double[numPoints];
    int *                            numTriangles = new int[numPoints];
    int *                            neighboor = new int[numPoints];
    for( int i = 0; i < numPoints; i++ )
      {
      avgTriangleArea[i] = 0.0;
      numTriangles[i] = 0;
      neighboor[i] = 0;
      }
    float  std = 0;
    double sum = 0;
    float  avgTriangle = 0;
    double max = 0, min = 0, dif = 0;
    double a = 0, b = 0, c = 0;
    bool   meshChanged = true;
    int    z = 0;

    while( meshChanged )
      {
      meshChanged = false;
      MeshType::PointType point;
      for( int dim = 0; dim < 3; dim++ )
        {
        point[dim] = 0;
        }

      // z!=0 load file (correctFilename) wich is already correct to improve the correction
      if( z != 0 )
        {
        MeshConverterType::Pointer       _converter = MeshConverterType::New();
        MeshSOType::Pointer              _SOMesh =
          dynamic_cast<MeshSOType *>(_converter->ReadMeta(correctFilename).GetPointer());
        MeshType::Pointer                _mesh = _SOMesh->GetMesh();
        MeshType::PointsContainerPointer _points = _mesh->GetPoints();
        }

      typedef CellType::PointIdIterator PointIdIterator;
      CellIterator cellIterator = mesh->GetCells()->Begin();
      CellIterator cellEnd      = mesh->GetCells()->End();
      // PointType    curPoint;

      int    _num = 0;
      // calculate average triangle for all points, max, min
      while( cellIterator != cellEnd )
        {
        CellType *                cell = cellIterator.Value();
        TriangleType *            line = dynamic_cast<TriangleType *>(cell);
        LineType::PointIdIterator pit = line->PointIdsBegin();
        MeshType::PointType       p1, p2, p3;
        int                       pIndex1, pIndex2, pIndex3;
        for( int dim = 0; dim < 3; dim++ )
          {
          pIndex1 = *pit;
          PointType _curPoint =  points->GetElement(pIndex1);
          p1[dim] = _curPoint[dim];
          }

        ++pit;
        for( int dim = 0; dim < 3; dim++ )
          {
          pIndex2 = *pit;
          PointType _curPoint =  points->GetElement(pIndex2);
          p2[dim] = _curPoint[dim];
          }

        ++pit;
        for( int dim = 0; dim < 3; dim++ )
          {
          pIndex3 = *pit;
          PointType _curPoint =  points->GetElement(pIndex3);
          p3[dim] = _curPoint[dim];
          }

        a = (p2[0] - p1[0]) * (p2[0] - p1[0]) + (p2[1] - p1[1]) * (p2[1] - p1[1]) + (p2[2] - p1[2]) * (p2[2] - p1[2]);
        c = (p3[0] - p1[0]) * (p3[0] - p1[0]) + (p3[1] - p1[1]) * (p3[1] - p1[1]) + (p3[2] - p1[2]) * (p3[2] - p1[2]);
        b = (p3[0] - p2[0]) * (p3[0] - p2[0]) + (p3[1] - p2[1]) * (p3[1] - p2[1]) + (p3[2] - p2[2]) * (p3[2] - p2[2]);

        double triangleArea = sqrt(2 * b * c + 2 * c * a + 2 * a * b - a * a - b * b - c * c) * 0.25;

        avgTriangleArea[pIndex1] += triangleArea;
        avgTriangleArea[pIndex2] += triangleArea;
        avgTriangleArea[pIndex3] += triangleArea;
        numTriangles[pIndex1]++;
        numTriangles[pIndex2]++;
        numTriangles[pIndex3]++;

        ++cellIterator;
        ++_num;
        }
      // calculate the average of the triangles
      for( int i = 0; i < numPoints; i++ )
        {
        avgTriangleArea[i] /= numTriangles[i];
        }
      // find max, min of all of the triangle area and calculate the total average triangle
      max = avgTriangleArea[0];
      min = avgTriangleArea[0];
      for( int i = 0; i < numPoints; i++ )
        {
        if( avgTriangleArea[i] > max )
          {
          max = avgTriangleArea[i];
          }
        if( avgTriangleArea[i] < min )
          {
          min = avgTriangleArea[i];
          }
        avgTriangle = avgTriangle + avgTriangleArea[i];
        }
      avgTriangle = avgTriangle / numPoints;
      // calculate the standard deviation
      for( int i = 0; i < numPoints; i++ )
        {
        dif = avgTriangleArea[i] - avgTriangle;
        sum = sum + dif * dif;
        }
      std = sqrt(sum / (numPoints - 1) );

      double areaThreshold = avgTriangle + threshFactor * std;
      int    n = 0;
      int    numNeighboor = 0;
      // correct points if it's necessary
      for( int i = 0; i < numPoints; i++ )
        {
        // if point i is a bad point
        if( avgTriangleArea[i]  >= areaThreshold )
          {
          for( int dim = 0; dim < 3; dim++ )
            {
            point[dim] = 0;
            }
          n = 0;
          numNeighboor = 0;
          for( int _c = 0; _c < numPoints; _c++ )
            {
            neighboor[_c] = 0;
            }
          PointType    curPoint =  points->GetElement(i);
          CellIterator _cellIterator = mesh->GetCells()->Begin();
          CellIterator _cellEnd      = mesh->GetCells()->End();
          bool         p2neighboor = 0, p3neighboor = 0;
          // find the neighboors and calculate the new point with the neighboors
          while( _cellIterator != _cellEnd )
            {
            CellType *                cell = _cellIterator.Value();
            TriangleType *            line = dynamic_cast<TriangleType *>(cell);
            LineType::PointIdIterator pit = line->PointIdsBegin();
            MeshType::PointType       p1, p2, p3, p;
            int                       pIndex1, pIndex2, pIndex3;
            // find the neighboors
            for( int dim = 0; dim < 3; dim++ )
              {
              pIndex1 = *pit;
              PointType _curPoint =  points->GetElement(pIndex1);
              p1[dim] = _curPoint[dim];
              }
            ++pit;
            for( int dim = 0; dim < 3; dim++ )
              {
              pIndex2 = *pit;
              PointType _curPoint =  points->GetElement(pIndex2);
              p2[dim] = _curPoint[dim];
              }
            ++pit;
            for( int dim = 0; dim < 3; dim++ )
              {
              pIndex3 = *pit;
              PointType _curPoint =  points->GetElement(pIndex3);
              p3[dim] = _curPoint[dim];
              }
            // if pIndex1 is the bad point
            if( pIndex1 == i )
              {
              // know if pIndex2 and pIndex3 are already neighboors
              for( int _b = 0; _b < n; _b++ )
                {
                if( pIndex2 == neighboor[_b] )
                  {
                  p2neighboor = 1;
                  }
                if( pIndex3 == neighboor[_b] )
                  {
                  p3neighboor = 1;
                  }
                }
              // nothing change
              if( p2neighboor && p3neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] = point[dim];
                  }
                }
              // if not calculate the new sum of neighboors, and the new num of neighboors
              if( !p3neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] += p3[dim];
                  }
                ++numNeighboor;
                }
              if( !p2neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] += p2[dim];
                  }
                ++numNeighboor;
                }
              neighboor[n] = pIndex2;
              ++n;
              neighboor[n] = pIndex3;
              ++n;
              }
            else if( pIndex2 == i )
              {
              // know if pIndex1 and pIndex3 are already neighboors
              for( int _b = 0; _b < n; _b++ )
                {
                if( pIndex1 == neighboor[_b] )
                  {
                  p2neighboor = 1;
                  }
                if( pIndex3 == neighboor[_b] )
                  {
                  p3neighboor = 1;
                  }
                }
              // nothing change
              if( p2neighboor && p3neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] = point[dim];
                  }
                }
              // if not calculate the new sum of neighboors, and the num of neighboors
              if( !p3neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] += p3[dim];
                  }
                ++numNeighboor;
                }
              if( !p2neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] += p1[dim];
                  }
                ++numNeighboor;
                }
              neighboor[n] = pIndex1;
              ++n;
              neighboor[n] = pIndex3;
              ++n;
              }
            else if( pIndex3 == i )
              {
              // know if pIndex1 and pIndex3 are already neighboors
              for( int _b = 0; _b < n; _b++ )
                {
                if( pIndex1 == neighboor[_b] )
                  {
                  p2neighboor = 1;
                  }
                if( pIndex2 == neighboor[_b] )
                  {
                  p3neighboor = 1;
                  }
                }
              // nothing change
              if( p2neighboor && p3neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] = point[dim];
                  }
                }
              // if not calculate the new sum of neighboors, and the num of neighboors
              if( !p3neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] += p2[dim];
                  }
                ++numNeighboor;
                }
              if( !p2neighboor )
                {
                for( int dim = 0; dim < 3; dim++ )
                  {
                  point[dim] += p1[dim];
                  }
                ++numNeighboor;
                }
              neighboor[n] = pIndex1;
              ++n;
              neighboor[n] = pIndex2;
              ++n;
              }
            ++_cellIterator;
            ++_num;
            }
          // calculate the new point : the average of the neighboors
          for( int dim = 0; dim < 3; dim++ )
            {
            point[dim] /= numNeighboor;
            }

          mesh->SetPoint(i, point);

          // write a new Mesh with correct triangles
          if( correctbadtriangleMeshOn )
            {
            SOMesh->SetMesh(mesh);
            MeshWriterType::Pointer writer = MeshWriterType::New();
            writer->SetInput(SOMesh);
            writer->SetFileName(correctFilename);
            writer->Update();
            }
          // know if the mesh has changed
          if( curPoint[0] != point[0] || curPoint[1] != point[1] || curPoint[2] != point[2] )
            {
            meshChanged = true;
            }
          }
        }
      ++z;
      }

    // write the correct Mesh with all correct triangles
    if( correctbadtriangleMeshOn )
      {
      SOMesh->SetMesh(mesh);
      MeshWriterType::Pointer writer = MeshWriterType::New();
      writer->SetInput(SOMesh);
      writer->SetFileName(correctFilename);
      writer->Update();
      }

    // write the correct values of the average of triangles in a text file
    if( !correctbadtriangleMeshOn )
      {
      // open outputFile
      std::ofstream outfileNor;
      outfileNor.open(outputFilename );
      outfileNor << "NUMBER_OF_POINTS=" << numPoints << std::endl;
      outfileNor << "DIMENSION=" << 1 << std::endl;
      outfileNor << "TYPE=Scalar" << std::endl;
      for( int i = 0; i < numPoints; i++ )
        {
        outfileNor << avgTriangleArea[i] << std::endl;
        }
      // close outputFile
      outfileNor.close();
      }

    delete[] avgTriangleArea;
    delete[] numTriangles;
    delete[] neighboor;

    }
  else if( valueOn )
    {

    char line[1000];

    // open inputFile
    std::ifstream input, input2;
    input2.open( inputFilename );
    input2.seekg(0, std::ios::beg);
    // calculate number of lines of the inputFilename
    int numberOfLines = 0;
    while( (input2.getline(line, 1000) ) != NULL )
      {
      numberOfLines++;
      }

    input2.close();
    // numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
    numberOfLines = numberOfLines - 3;

    // open inputFile
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);

    // create output (KWMeshVisu textFile)
    std::ofstream output;
    output.open( outputFilename );
    output << "NUMBER_OF_POINTS= " << numberOfLines << std::endl;
    output << "DIMENSION=1" << std::endl;
    output << "TYPE=Scalar" << std::endl;

    // extract column 5 and write the values
    while( (input.getline(line, 1000) ) != NULL )
      {
      float dummy, extractValue, point;
      sscanf(line, " %f %f %f %f %f ", &point, &dummy, &dummy, &dummy, &extractValue);
      // ////////uncomment if you want to select with the first column (like only one label)
      // if(point == Index)
//  output <<extractValue << std::endl;
      }

    // close input&output
    input.close();
    output.close();

    }
  else if( labelOn )
    {

    char   line[1000], line1[1000], line2[1000];
    char * valuePtr = NULL;
    int    nPts;

    std::ifstream labelFile;
    labelFile.open( labelFilename );
    labelFile.seekg(0, std::ios::beg);
    bool found = false;
    while( !found && !labelFile.eof() )
      {
      labelFile.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "NUMBER_OF_POINTS" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nPts);
    labelFile.close();

    std::ofstream output;
    output.open( outputFilename );

    std::ifstream inputFile;
    inputFile.open( inputFilename );
    inputFile.seekg(0, std::ios::beg);

    std::ifstream labelFile1;
    labelFile1.open( labelFilename );
    labelFile1.seekg(0, std::ios::beg);
    std::vector<float> label(nPts);
    std::vector<float> thick(nPts);
    int   i = 0;

    found = false;
    while( !found && !labelFile1.eof() )
      {
      labelFile1.getline(line1, 1000);
      if( strstr( line1, "TYPE=Scalar" ) )
        {
        found = true;
        }
      }

    while( (labelFile1.getline(line1, 1000) ) != NULL )
      {
      float Value;
      sscanf(line1, " %f ", &Value);
      label[i] = Value;
      i++;
      }

    found = false;
    while( !found && !inputFile.eof() )
      {
      inputFile.getline(line2, 1000);
      if( strstr( line2, "TYPE=Scalar" ) )
        {
        found = true;
        }
      }

    i = 0;
    while( (inputFile.getline(line2, 1000) ) != NULL )
      {
      float Value1;
      sscanf(line2, " %f ", &Value1);
      thick[i] = Value1;
      i++;
      }

    // /////separeted all labels in 1file
    //     int nu=0;
    //     float label_new[nPts];
    //     int d=0;
    //     float num[nPts];
    //     for (int c = 0; c < nPts; c++)
    //       {
    //  for(int a =0; a<29;a++)
    //    {
    //      if(label[c]==a)
    //        {
    //    num[a]++;
    //        }
    //    }
    //       }

    //     for(int b =0; b<29;b++)
    //       {
    //  if(num[b]!=0)
    //    {
    //      label_new[d]=b;
    //      d++;
    //    }
    //       }

    //       ///////num labels
    //       //     for(int g =0; g<25;g++)
    //       //       {
    //       //   std::cout <<g <<" " << label_new[g] <<std::endl;
    //       //       }

    //     for(int a =0; a<25;a++)
    //       {
    //  nu=0;
    //  output << label_new[a] <<std::endl;
    //  for (int c = 0; c < nPts; c++)
    //    {
    //      if(label[c]==label_new[a])
    //        {
    //    //std::cout << thick[c] <<",";
    //    output << thick[c] <<",";
    //    nu++;
    //        }
    //    }
    //  output << "" <<std::endl;
    //  output << nu <<std::endl;
    //  output << "" <<std::endl;
    //       }

    // //// average all labels :
    //     float thicknew[nPts];
    //     float num[nPts];
    //     for (int c = 0; c < nPts; c++)
    //       {
    //  for(int a =0; a<29;a++)
    //    {
    //      if(label[c]==a)
    //        {
    //    thicknew[a]=thick[c]+thicknew[a];
    //    num[a]++;
    //        }
    //    }
    //       }
    //     for(int b =0; b<29;b++)
    //       {
    //  if(num[b]!=0)

    //    {
    //      thicknew[b]=thicknew[b]/num[b];
    //      output <<thicknew[b] <<std::endl;
    //    }
    //       }

    // ////mean for one label : label == 11
    output << "NUMBER_OF_POINTS=" << 1 << std::endl;
    output << "DIMENSION=1 " << std::endl;
    output << "TYPE=Scalar " << std::endl;

    float thicknew = 0;
    {
    int   _num = 0;
    for( int c = 0; c < nPts; c++ )
      {
      if( label[c] == 11 )
        {
        thicknew = thick[c] + thicknew;
        _num++;
        // every values of one label : label == 11
        // output <<thick[c] << std::endl;
        }
      }
    thicknew = thicknew / _num;
    }
    output << thicknew << std::endl;

    inputFile.close();
    labelFile1.close();
    output.close();

    }
  else if( colorOn )
    {
    char line[1000];
    bool found = false;

    std::ifstream input;
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);

    std::ofstream output;
    output.open( outputFilename );

    while( !found && !input.eof() )
      {
      input.getline(line, 1000);
      output << line << std::endl;
      if( strstr( line, "LOOKUP_TABLE default" ) )
        {
        found = true;
        }
      }

    float *val = new float[30];
    int    max = 0;
    {
    int    i = 1;
    for( int y = 1; y <= num[0]; y++ )
      {
      val[i] = num[y];
      i++;
      }

    val[i] = num[1];
    i++;
    val[i] = num[1];
    i++;
    val[i] = num[1];
    i++;
    for( int z = 1; z <= old_num[0]; z++ )
      {
      val[i] = old_num[z];
      i++;
      max = i;
      }
    }
    while( (input.getline(line, 1000) ) != NULL )
      {
      float a, b, c, d, e, f, g, h, i;
      bool  Afound = false, Bfound = false, Cfound = false, Dfound = false, Efound = false, Ffound = false, Gfound =
        false, Hfound = false, Ifound = false;
      sscanf(line, " %f %f %f %f %f %f %f %f %f", &a, &b, &c, &d, &e, &f, &g, &h, &i );
      for( int q = 1; q < max; q++ )
        {
        if( a == val[q] && Afound == false )
          {
          a = q;
          Afound = true;
          }
        if( b == val[q] && Bfound == false )
          {
          b = q;
          Bfound = true;
          }
        if( c == val[q] && Cfound == false )
          {
          c = q;
          Cfound = true;
          }
        if( d == val[q] && Dfound == false )
          {
          d = q;
          Dfound = true;
          }
        if( e == val[q] && Efound == false )
          {
          e = q;
          Efound = true;
          }
        if( f == val[q] && Ffound == false )
          {
          f = q;
          Ffound = true;
          }
        if( g == val[q] && Gfound == false )
          {
          g = q;
          Gfound = true;
          }
        if( h == val[q] && Hfound == false )
          {
          h = q;
          Hfound = true;
          }
        if( i == val[q] && Ifound == false )
          {
          i = q;
          Ifound = true;
          }
        }

      if( Afound == false )
        {
        a = 0;
        }
      if( Bfound == false )
        {
        b = 0;
        }
      if( Cfound == false )
        {
        c = 0;
        }
      if( Dfound == false )
        {
        d = 0;
        }
      if( Efound == false )
        {
        e = 0;
        }
      if( Ffound == false )
        {
        f = 0;
        }
      if( Gfound == false )
        {
        g = 0;
        }
      if( Hfound == false )
        {
        h = 0;
        }
      if( Ifound == false )
        {
        i = 0;
        }

      output << a << " " << b << " " << c << " " << d << " " << e << " " << f << " " << g << " " << h << " " << i
             << std::endl;
      }

    delete [] val;
    delete [] num;
    delete [] old_num;
    input.close();
    output.close();

    }
  else if( MaxColorOn )
    {

    char line[1000], line2[1000];
    bool found = false;

    std::ifstream in;
    in.open( inputFilename );
    in.seekg(0, std::ios::beg);
    int numberOfLines = 0;
    // find number of lines
    while( (in.getline(line, 1000) ) != NULL )
      {
      numberOfLines++;
      }

    in.close();
    // numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
    numberOfLines = numberOfLines - 3;

    std::ifstream input;
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);
    // pass the three first lines of KWMeshVisu textfile
    while( !found && !input.eof() )
      {
      input.getline(line, 1000);
      if( strstr( line, "TYPE=Scalar" ) )
        {
        found = true;
        }
      }

    // numFiles corresponds to number of files -1 (without the inputfile)
    const int numFiles = MaxFiles.size();

    std::vector<std::ifstream *> attrFile(numFiles);

    double * *value = NULL;
    value = new double *[numberOfLines];
    for( int c = 0; c < numberOfLines; c++ )
      {
      value[c] = new double[numFiles + 1];
      for( int d = 0; d < numFiles + 1; d++ )
        {
        value[c][d] = 0.0;
        }
      }

    int val = 0;
    std::ofstream output;
    std::ofstream output2;
    for( int _num = 0; _num < numFiles; _num++ )
      {
      attrFile[_num]= new std::ifstream();
      attrFile[_num]->open(MaxFiles[_num].c_str());
      attrFile[_num]->seekg(0, std::ios::beg);
      for( val = 0; val < 2; val++ )
        {
        attrFile[_num]->getline(line2, 1000);
        }
      }

    int index = 0;
    while( (input.getline(line, 1000) ) != NULL )
      {
      float point;
      sscanf(line, " %f ", &point);
      value[index][0] = point;
      for( int _num = 0; _num < numFiles; _num++ )
        {
        float point2;
        attrFile[_num]->getline(line2, 1000);
        sscanf(line2, " %f ", &point2);
        value[index][_num + 1] = point2;
        }
      ++index;
      }
    for( int _num = 0; _num < numFiles; _num++ )
      {
      attrFile[_num]->close();
      delete attrFile[_num];
      }
    input.close();
    for( int a = 0; a < numberOfLines; a++ )
      {
      float max = 0.0, maxI = 0.0;
      for( int b = 0; b < numFiles + 1; b++ )
        {
        if( value[a][b] > max )
          {
          max = value[a][b];
          }
        }

      maxI = max - max * 5 / 100;
      for( int e = 0; e < numFiles + 1; e++ )
        {
        if( value[a][e] <= maxI )
          {
          value[a][e] = 0.0;
          }
        }
      }

    // write in the first file -> inputfile
    output.open(inputFilename);
    output << "NUMBER_OF_POINTS=" << numberOfLines << std::endl;
    output << "DIMENSION=1" << std::endl;
    output << "TYPE=Scalar" << std::endl;
    for( int z = 0; z < numberOfLines; z++ )
      {
      output << value[z][0] << std::endl;
      }
    output.close();
    // write in the other files
    for( int x = 0; x < numFiles + 1; x++ )
      {
      output2.open(MaxFiles[x].c_str() );
      output2 << "NUMBER_OF_POINTS=" << numberOfLines << std::endl;
      output2 << "DIMENSION=1" << std::endl;
      output2 << "TYPE=Scalar" << std::endl;
      for( int y = 0; y < numberOfLines; y++ )
        {
        output2 << value[y][x + 1] << std::endl;
        }
      output2.close();
      }
    for( int c = 0; c < numberOfLines; c++ )
      {
      delete [] value[c];
      }
    delete [] value;

    }
  else if( distAbsOn )
    {

    char line[1000], line2[1000];
    bool found = false;

    std::ifstream in;
    in.open( inputFilename );
    in.seekg(0, std::ios::beg);
    int numberOfLines = 0;
    // find number of lines
    while( (in.getline(line, 1000) ) != NULL )
      {
      numberOfLines++;
      }

    in.close();
    // numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
    numberOfLines = numberOfLines - 3;

    std::ifstream input;
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);
    // pass the three first lines of KWMeshVisu textfile

    while( !found && !input.eof() )
      {
      input.getline(line, 1000);
      if( strstr( line, "TYPE=Scalar" ) )
        {
        found = true;
        }
      }

    int       index = 0;
    int       val = 0;
    const int       numFiles = distFiles.size();
    double * *value = NULL;
    value = new double *[numberOfLines];
    for( int c = 0; c < numberOfLines; c++ )
      {
      value[c] = new double[numFiles + 1];
      for( int d = 0; d < numFiles + 1; d++ )
        {
        value[c][d] = 0.0;
        }
      }
    double * *value_new = NULL;
    value_new = new double *[numberOfLines];
    for( int u = 0; u < numberOfLines; u++ )
      {
      value_new[u] = new double[numFiles + 1];
      for( int k = 0; k < numFiles + 1; k++ )
        {
        value_new[u][k] = 0.0;
        }
      }
    std::vector<std::ifstream *> attrFile(numFiles);
    for( int _num = 0; _num < numFiles; _num++ )
      {
      attrFile[_num] = new std::ifstream();
      attrFile[_num]->open(distFiles[_num].c_str() );
      attrFile[_num]->seekg(0, std::ios::beg);
      for( val = 0; val < 2; val++ )
        {
        attrFile[_num]->getline(line2, 1000);
        }
      }

    while( (input.getline(line, 1000) ) != NULL )
      {
      float point;
      sscanf(line, " %f ", &point);
      value[index][0] = point;
      for( int _num = 0; _num < numFiles; _num++ )
        {
        float point2;
        attrFile[_num]->getline(line2, 1000);
        sscanf(line2, " %f ", &point2);
        value[index][_num + 1] = point2;
        }
      ++index;
      }
    for( int _num = 0; _num < numFiles; _num++ )
      {
      attrFile[_num]->close();
      delete attrFile[_num];
      }
    input.close();
    for( int a = 0; a < numberOfLines; a++ )
      {
      for( int b = 1; b < numFiles; b++ )
        {
        value_new[a][b] = value[a][b + 1] - value[a][b - 1];
        }
      }

    // write output file
    std::ofstream output;
    for( int x = 0; x < numFiles; x++ )
      {
      output.open(resultFiles[x].c_str() );
      output << "NUMBER_OF_POINTS=" << numberOfLines << std::endl;
      output << "DIMENSION=1" << std::endl;
      output << "TYPE=Scalar" << std::endl;
      for( int p = 0; p < numberOfLines; p++ )
        {
        output << value_new[p][x + 1] << std::endl;
        }

      output.close();
      }
    for( int c = 0; c < numberOfLines; c++ )
      {
      delete [] value[c];
      delete [] value_new[c];
      }
    delete [] value;
    delete [] value_new;

    }
  else if( distRelOn )
    {

    char line[1000], line2[1000];
    bool found = false;

    std::ifstream in;
    in.open( inputFilename );
    in.seekg(0, std::ios::beg);
    int numberOfLines = 0;
    // find number of lines
    while( (in.getline(line, 1000) ) != NULL )
      {
      numberOfLines++;
      }

    in.close();
    // numberOfLines-3 corresponds to the number of points (without the three first lines of KWMeshVisu textfile)
    numberOfLines = numberOfLines - 3;

    std::ifstream input;
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);
    // pass the three first lines of KWMeshVisu textfile

    while( !found && !input.eof() )
      {
      input.getline(line, 1000);
      if( strstr( line, "TYPE=Scalar" ) )
        {
        found = true;
        }
      }

    int       index = 0;
    int       val = 0;
    int       numFiles = distReFiles.size();
    double * *value = NULL;
    value = new double *[numberOfLines];
    for( int c = 0; c < numberOfLines; c++ )
      {
      value[c] = new double[numFiles + 1];
      for( int d = 0; d < numFiles + 1; d++ )
        {
        value[c][d] = 0.0;
        }
      }
    double * *value_new = NULL;
    value_new = new double *[numberOfLines];
    for( int u = 0; u < numberOfLines; u++ )
      {
      value_new[u] = new double[numFiles + 1];
      for( int k = 0; k < numFiles + 1; k++ )
        {
        value_new[u][k] = 0.0;
        }
      }
    std::vector<std::ifstream *> attrFile(numFiles);
    for( int _num = 0; _num < numFiles; _num++ )
      {
      attrFile[_num] = new std::ifstream();
      attrFile[_num]->open(distReFiles[_num].c_str() );
      attrFile[_num]->seekg(0, std::ios::beg);
      for( val = 0; val < 2; val++ )
        {
        attrFile[_num]->getline(line2, 1000);
        }
      }

    while( (input.getline(line, 1000) ) != NULL )
      {
      float point;
      sscanf(line, " %f ", &point);
      value[index][0] = point;
      for( int _num = 0; _num < numFiles; _num++ )
        {
        float point2;
        attrFile[_num]->getline(line2, 1000);
        sscanf(line2, " %f ", &point2);
        value[index][_num + 1] = point2;
        }
      ++index;
      }
    for( int _num = 0; _num < numFiles; _num++ )
      {
      attrFile[_num]->close();
      delete attrFile[_num];
      }
    input.close();
    for( int a = 0; a < numberOfLines; a++ )
      {
      for( int b = 1; b < numFiles; b++ )
        {
        value_new[a][b] = value[a][b + 1] - value[a][b - 1];
        }

      float max = 0.0;
      float min = 10.0;
      for( int h = 1; h < numFiles; h++ )
        {
        if( value_new[a][h] > max )
          {
          max = value_new[a][h];
          }
        if( value[a][h] < min )
          {
          min = value_new[a][h];
          }
        }
      for( int g = 1; g < numFiles; g++ )
        {
        value_new[a][g] = (value_new[a][g] - min) / (max - min);
        }
      }

    // write output file
    std::ofstream output;
    for( int x = 0; x < numFiles; x++ )
      {
      output.open(restFiles[x].c_str() );
      output << "NUMBER_OF_POINTS=" << numberOfLines << std::endl;
      output << "DIMENSION=1" << std::endl;
      output << "TYPE=Scalar" << std::endl;
      for( int p = 0; p < numberOfLines; p++ )
        {
        output << value_new[p][x + 1] << std::endl;
        }

      output.close();
      }
    for( int c = 0; c < numberOfLines; c++ )
      {
      delete [] value[c];
      delete [] value_new[c];
      }
    delete [] value;
    delete [] value_new;

    }
  else if( firstOn )
    {

    char line[1000], line2[1000];
    int  numfile = firstfiles.size();

    std::ifstream input;
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);

    std::ofstream output;
    output.open( outputFilename );
    output << "NUMBER_OF_POINTS=" << numfile + 2 << std::endl;
    output << "DIMENSION=1" << std::endl;
    output << "TYPE=Scalar" << std::endl;
    const int numEntries = 3;
    int       counter = 0;

    while( counter < numEntries && !input.eof() )
      {
      input.getline( line, 1000 );
      if( (line[0] != '#') )
        {
        counter++;
        }
      }

    while( (input.getline(line, 1000) ) != NULL )
      {
      float point;
      sscanf(line, " %f ", &point);
      output << point << ",";
      }

    input.close();
    for( int index = 0; index < numfile; index++ )
      {

      std::ifstream input2;
      counter = 0;
      input2.open(firstfiles[index].c_str() );
      input2.seekg(0, std::ios::beg);

      while( counter < numEntries && !input2.eof() )
        {
        input2.getline( line2, 1000 );
        if( (line2[0] != '#') )
          {
          counter++;
          }
        }

      while( (input2.getline(line2, 1000) ) != NULL )
        {
        float point2;
        sscanf(line2, " %f ", &point2);
        output << point2 << ",";
        }

      input2.close();
      }
    output << " " << std::endl;
    output.close();

    }
  else if( subKWMOn )
    {

    char   line1[1000], line2[1000];
    float  value1 = 0.0;
    float  value2 = 0.0;
    float  diff = 0.0;
    bool   found;
    char * valuePtr = NULL;
    int    nPts1, nPts2, nDim;
    char   typeString[1001], line[1001];

    // open input
    std::ifstream input1, input2;
    input1.open( inputFilename );
    input1.seekg(0, std::ios::beg);
    found = false;
    // find number of points
    while( !found && !input1.eof() )
      {
      input1.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "NUMBER_OF_POINTS" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nPts1);

    input1.seekg(0, std::ios::beg);
    found = false;
    while( !found && !input1.eof() )
      {
      input1.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "DIMENSION" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nDim);

    input1.seekg(0, std::ios::beg);
    found = false;
    while( !found && !input1.eof() )
      {
      input1.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "TYPE" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %s ", typeString);
    assert( nDim == 1 );
    assert( strcmp( typeString, "Scalar" ) == 0 );

    // open subKWMFile
    input2.open( subKWMFilename );
    input2.seekg(0, std::ios::beg);
    found = false;
    // find the number of points
    while( !found && !input2.eof() )
      {
      input2.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "NUMBER_OF_POINTS" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nPts2);

    input2.seekg(0, std::ios::beg);
    found = false;
    while( !found && !input2.eof() )
      {
      input2.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "DIMENSION" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nDim);

    input2.seekg(0, std::ios::beg);
    found = false;
    while( !found && !input2.eof() )
      {
      input2.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "TYPE" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %s ", typeString);
    // number of points of the first file = number of points of the second file
    assert( nPts2 == nPts1 );
    assert( nDim == 1 );
    assert( strcmp( typeString, "Scalar" ) == 0 );

    // create ouput file
    std::ofstream output;
    output.open( outputFilename );
    output << "NUMBER_OF_POINTS=" << nPts1 << std::endl;
    output << "DIMENSION=" << 1 << std::endl;
    output << "TYPE=Scalar" << std::endl;
    // find the values of the two files, make the difference and write them into the ouputFile
    while( ( (input1.getline(line1, 1000) ) != NULL) && (input2.getline(line2, 1000) != NULL) )
      {
      value1 = atof(line1);
      value2 = atof(line2);
      diff = value1 - value2;
      output << diff << std::endl;
      }

    // close inputs and ouput
    input1.close();
    input2.close();
    output.close();

    }
  else if( avgGaussKWMOn )
    {

    float  sum = 0;
    float  gaussKWM = 0;
    int    nPts;
    int    nPtsFile, nDim;
    bool   found;
    char * valuePtr = NULL;
    char   typeString[1000], line[1000];
    double Gaussvalue;

    // open file(input)
    std::ifstream input;
    input.open( inputFilename );
    input.seekg(0, std::ios::beg);
    found = false;
    // find the number of points
    while( !found && !input.eof() )
      {
      input.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "NUMBER_OF_POINTS" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nPtsFile);
    nPts = atoi(valuePtr);
    input.seekg(0, std::ios::beg);
    found = false;
    while( !found && !input.eof() )
      {
      input.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "DIMENSION" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %d ", &nDim);
    input.seekg(0, std::ios::beg);
    found = false;
    while( !found && !input.eof() )
      {
      input.getline( line, 1000 );
      if( line[0] != '#' && strstr( line, "TYPE" ) )
        {
        found = true;
        }
      }

    valuePtr = strchr(line, '=');
    if( !valuePtr )
      {
      return -1;
      }
    valuePtr++;
    sscanf(valuePtr, " %s ", typeString);
    assert( nPtsFile == nPts );
    assert( nDim == 1 );
    assert( strcmp( typeString, "Scalar" ) == 0 );

    float *avgGaussKWM = new float[nPts];
    float  val;

    input.seekg(0, std::ios::beg);
    const int numEntries = 3;
    int       counter = 0;
    while( counter < numEntries && !input.eof() )
      {
      input.getline( line, 1000 );
      if( (line[0] != '#') )
        {
        counter++;
        }
      }

    // calculate gauss coefficient
    gaussKWM =
      (1
       / (GaussKWMParams[1]
          * sqrt(2
                 * M_PI) )
       * exp(-(AgeKWM[0]
               - GaussKWMParams[0])
             * (AgeKWM[0] - GaussKWMParams[0]) / (GaussKWMParams[1] * GaussKWMParams[1] * 2) ) ) / 0.0797885;
    sum = gaussKWM;
    // find values and calculate the gauss values
    for( int i = 0; i < nPts; i++ )
      {
      input.getline(line, 1000);
      val = atof(line);
      avgGaussKWM[i] = gaussKWM * val;
      }
    // close file
    input.close();

    int   numFiles = avgGaussKWMFiles.size();
    float newval;
    // for all files
    for( int index = 0; index < numFiles; index++ )
      {

      // open file
      input.open(avgGaussKWMFiles[index].c_str() );
      input.seekg(0, std::ios::beg);
      found = false;
      // find the number of points
      while( !found && !input.eof() )
        {
        input.getline( line, 1000 );
        if( line[0] != '#' && strstr( line, "NUMBER_OF_POINTS" ) )
          {
          found = true;
          }
        }

      valuePtr = strchr(line, '=');
      if( !valuePtr )
        {
        return -1;
        }
      valuePtr++;
      sscanf(valuePtr, " %d ", &nPtsFile);
      nPts = atoi(valuePtr);
      input.seekg(0, std::ios::beg);
      found = false;
      while( !found && !input.eof() )
        {
        input.getline( line, 1000 );
        if( line[0] != '#' && strstr( line, "DIMENSION" ) )
          {
          found = true;
          }
        }

      valuePtr = strchr(line, '=');
      if( !valuePtr )
        {
        return -1;
        }
      valuePtr++;
      sscanf(valuePtr, " %d ", &nDim);
      input.seekg(0, std::ios::beg);
      found = false;
      while( !found && !input.eof() )
        {
        input.getline( line, 1000 );
        if( line[0] != '#' && strstr( line, "TYPE" ) )
          {
          found = true;
          }
        }

      valuePtr = strchr(line, '=');
      if( !valuePtr )
        {
        return -1;
        }
      valuePtr++;
      sscanf(valuePtr, " %s ", typeString);
      assert( nPtsFile == nPts );
      assert( nDim == 1 );
      assert( strcmp( typeString, "Scalar" ) == 0 );
      input.seekg(0, std::ios::beg);
      const int _numEntries = 3;
      int       _counter = 0;
      while( _counter < _numEntries && !input.eof() )
        {
        input.getline( line, 1000 );
        if( (line[0] != '#') )
          {
          _counter++;
          }
        }

      // calculate gauss coefficient
      gaussKWM =
        (1
         / (GaussKWMParams[1]
            * sqrt(2
                   * M_PI) )
         * exp(-(AgeKWM[index
                        + 1]
                 - GaussKWMParams[0])
               * (AgeKWM[index + 1] - GaussKWMParams[0]) / (GaussKWMParams[1] * GaussKWMParams[1] * 2) ) ) / 0.0797885;
      sum = sum + gaussKWM;
      // find values and calculate the gauss values
      for( int i = 0; i < nPts; i++ )
        {
        input.getline(line, 1000);
        newval = atof(line);
        Gaussvalue = gaussKWM * newval;
        avgGaussKWM[i] = avgGaussKWM[i] + Gaussvalue;
        }

      // close file
      input.close();

      }
    // calculate average (divide the values by sum)
    for( int i = 0; i < nPts; i++ )
      {
      avgGaussKWM[i] = avgGaussKWM[i] / sum;
      }
    // write KWMeshVisu textFile with new values
    std::ofstream output;
    output.open( outputFilename );
    output << "NUMBER_OF_POINTS=" << nPts << std::endl;
    output << "DIMENSION=" << 1 << std::endl;
    output << "TYPE=Scalar" << std::endl;
    for( int i = 0; i < nPts; i++ )
      {
      output << avgGaussKWM[i] << std::endl;
      }
    // close output
    output.close();

    delete [] GaussKWMParams;
    delete [] AgeKWM;
    delete [] avgGaussKWM;

    }
  else if( avgGaussMeshOn )
    {

    float sum = 0;
    float gauss = 0;

    MeshConverterType::Pointer converter = MeshConverterType::New();
    // read input
    MeshSOType::Pointer        SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer          surfaceMesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer avgPoints = surfaceMesh->GetPoints();
    int                        numberOfPoints = avgPoints->Size();
      {
      MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
      // calculate the gauss coefficient with mean,std,age...
      gauss =
        (1
         / (GaussParams[1]
            * sqrt(2
                   * M_PI) )
         * exp(-(Age[0]
                 - GaussParams[0]) * (Age[0] - GaussParams[0]) / (GaussParams[1] * GaussParams[1] * 2) ) ) / 0.0797885;
      // calculate the average values
      for( unsigned int pointID = 0; pointID < avgPoints->Size(); pointID++ )
        {
        PointType vert;
        PointType curPointAvg =  avgPoints->GetElement(pointID);
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          vert[dim] = gauss * curPointAvg[dim];
          }
        pointsTmp->InsertElement(pointID, vert);
        }
      sum = gauss;
      avgPoints = pointsTmp;
      }

    // for all meshes
    int numMeshes = AvgGaussMeshFiles.size();
    for( int index = 0; index < numMeshes; index++ )
      {
      try
        {
        MeshSOType::Pointer              inputMeshSO =
          dynamic_cast<MeshSOType *>(converter->ReadMeta(AvgGaussMeshFiles[index].c_str() ).GetPointer() );
        MeshType::Pointer                inputMesh = inputMeshSO->GetMesh();
        MeshType::PointsContainerPointer points = inputMesh->GetPoints();
        MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
        // calculate the gauss coefficient with mean,std,age...
        gauss =
          (1
           / (GaussParams[1]
              * sqrt(2
                     * M_PI) )
           * exp(-(Age[index
                       + 1]
                   - GaussParams[0])
                 * (Age[index + 1] - GaussParams[0]) / (GaussParams[1] * GaussParams[1] * 2) ) ) / 0.0797885;
        sum = sum + gauss;

        if( points->Size() != (unsigned int) numberOfPoints )
          {
          std::cerr << "Number of Points do not agree, expected : " << numberOfPoints << " , got " << points->Size()
                    << " File: " << AvgGaussMeshFiles[index] <<  std::endl;
          exit(-1);
          }
        // calculate the average values
        for( unsigned int pointID = 0; pointID < points->Size(); pointID++ )
          {
          PointType curPoint =  points->GetElement(pointID);
          PointType vert;
          PointType curPointAvg =  avgPoints->GetElement(pointID);
          for( unsigned int dim = 0; dim < 3; dim++ )
            {
            vert[dim] = gauss * curPoint[dim] + curPointAvg[dim];
            }
          pointsTmp->InsertElement(pointID, vert);
          }
        avgPoints = pointsTmp;
        }
      catch( itk::ExceptionObject ex )
        {
        std::cerr << "Error reading meshfile:  " << AvgGaussMeshFiles[index] << std::endl << "ITK error: "
                  << ex.GetDescription() << std::endl;
        exit(-3);
        }
      }
    // calculate the average (values divide by the sum)
    MeshType::PointsContainerPointer pointsTmp = MeshType::PointsContainer::New();
    for( unsigned int pointID = 0; pointID < avgPoints->Size(); pointID++ )
      {
      PointType curPoint =  avgPoints->GetElement(pointID);
      PointType vert;
      for( unsigned int dim = 0; dim < 3; dim++ )
        {
        vert[dim] = curPoint[dim] / sum;
        }
      pointsTmp->InsertElement(pointID, vert);
      }
    avgPoints = pointsTmp;

    // create new Mesh with the new values
    surfaceMesh->SetPoints(avgPoints);
    SOMesh->SetMesh(surfaceMesh);
    MeshWriterType::Pointer writer = MeshWriterType::New();
    writer->SetInput(SOMesh);
    writer->SetFileName(outputFilename);
    writer->Update();

    delete [] GaussParams;
    delete [] Age;

    }
  else if( averageOn )
    {

    int NbAveFile = AveFiles.size();

    // read input
    MeshConverterType::Pointer converter = MeshConverterType::New();
    if( debug )
      {
      std::cout << "Reading input mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer inputMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer   inputMesh = inputMeshSO->GetMesh();
    int                 NbOfPoint = inputMesh->GetNumberOfPoints();

    // Convert the meta data into the VTK Poly Data format
    if( debug )
      {
      std::cout << "Computing normals" << std::endl;
      }
    TriangleMeshConverterType::Pointer Triangleconverter = TriangleMeshConverterType::New();
    TriangleMeshSOType::Pointer inputTriangleMeshSO =
     dynamic_cast<TriangleMeshSOType *>(Triangleconverter->ReadMeta(inputFilename).GetPointer());
    TriangleMeshType::Pointer   inputTriangleMesh = inputTriangleMeshSO->GetMesh();

    itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
    convertMeshToVTK->SetInput(inputTriangleMesh);

    vtkPolyData *       vtkMesh = convertMeshToVTK->GetOutput();
    vtkPolyDataNormals *MeshNormals = vtkPolyDataNormals::New();

    MeshNormals->SetComputePointNormals(1);
    MeshNormals->SetComputeCellNormals(0);
    MeshNormals->SetSplitting(0);
    #if VTK_MAJOR_VERSION > 5
    MeshNormals->SetInputData(vtkMesh);
    #else
    MeshNormals->SetInput(vtkMesh);
    MeshNormals->Update();//no more Update() function for vtkPolyData: http://www.vtk.org/Wiki/VTK/VTK_6_Migration/Removal_of_Update
    #endif
    vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    vtkMeshNormals->Update();
    #endif
    vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
    vtkDataArray * ArrayNormal = NormalPoints->GetNormals();

    char   output[512];
    float *AverageX = new float[NbOfPoint];
    float *AverageY = new float[NbOfPoint];
    float *AverageZ = new float[NbOfPoint];
    if( debug )
      {
      std::cout << "Vector initialization" << std::endl;
      }
    for( int i = 0; i < NbOfPoint; i++ )
      {
      AverageX[i] = AverageY[i] = AverageZ[i] = 0.0;
      }
    for( int FileId = 0; FileId < (int)(AveFiles.size() ); FileId++ )
      {
      if( debug )
        {
        std::cout << "Reading file number " << FileId << std::endl;
        }
      // Open the file
      std::ifstream avefile;
      avefile.open(AveFiles[FileId].c_str() );
      // Get the number of point in the file
      avefile.getline(output, 32);
      int Val = GetIntFromString(output);
      if( Val != NbOfPoint )
        {
        std::cerr << "The file must have the same number of points as the meta file" << std::endl;
        exit(-1);
        }
      else
        {

        // Read the two next lines
        avefile.getline(output, 64);
        avefile.getline(output, 64);
        // Then get the numbers
        for( int Line = 0; Line < NbOfPoint; Line++ )
          {
          float X, Y, Z;
          avefile.getline(output, 64);
          sscanf(output, "%f %f %f", &X, &Y, &Z);
          AverageX[Line] += X;
          AverageY[Line] += Y;
          AverageZ[Line] += Z;
          }
        }
      }

    if( debug )
      {
      std::cout << "Writing outputs" << std::endl;
      }
    // Computing the average of each vector and writing it down in the file along with the norm values
    std::ofstream outfileVec;
    outfileVec.open( outputFilename );
    // print the header
    outfileVec << "NUMBER_OF_POINTS = " << NbOfPoint << std::endl;
    outfileVec << "DIMENSION = " << 3 << std::endl;
    outfileVec << "TYPE = Vector" << std::endl;

    outputFilename2 = std::string(outputFilename);
    DotLoc = outputFilename2.find(".txt", 0);
    if( DotLoc != -1 )
      {
      outputFilename2.insert(DotLoc, ".mag");
      }
    else
      {
      outputFilename2 = std::string(outputFilename) + std::string(".mag");
      }

    std::ofstream outfileNor;
    outfileNor.open( outputFilename2.c_str() );
    // print the header
    outfileNor << "NUMBER_OF_POINTS=" << NbOfPoint << std::endl;
    outfileNor << "DIMENSION=" << 1 << std::endl;
    outfileNor << "TYPE=Scalar" << std::endl;

    double *Normal = new double[3];
    for( int Line = 0; Line < NbOfPoint; Line++ )
      {
      float MeanX = AverageX[Line] / NbAveFile;
      float MeanY = AverageY[Line] / NbAveFile;
      float MeanZ = AverageZ[Line] / NbAveFile;
      // if the mean vector is projected on the normal at each point
      if( normaveOn )
        {
        Normal = ArrayNormal->GetTuple3(Line);
        float DotProd = (float)Normal[0] * MeanX + (float)Normal[1] * MeanY + (float)Normal[2] * MeanZ;
        float MeanXProj = (float)Normal[0] * DotProd;
        float MeanYProj = (float)Normal[1] * DotProd;
        float MeanZProj = (float)Normal[2] * DotProd;
        outfileVec << MeanXProj << " " << MeanYProj << " " << MeanZProj << std::endl;
        float Norme = sqrt(MeanXProj * MeanXProj + MeanYProj * MeanYProj + MeanZProj * MeanZProj);
        outfileNor << Norme << std::endl;
        }
      else
        {
        outfileVec << MeanX << " " << MeanY << " " << MeanZ << std::endl;
        float Norme = sqrt(MeanX * MeanX + MeanY * MeanY + MeanZ * MeanZ);
        outfileNor << Norme << std::endl;
        }
      }

    delete [] AverageX;
    delete [] AverageY;
    delete [] AverageZ;
    delete [] Normal;
    outfileVec.close();
    outfileNor.close();

    }
  else if( invvectOn )
    {

    std::ifstream vectorFile;
    vectorFile.open(VectFile);
    std::ofstream outfileVec;
    outfileVec.open( outputFilename );
    // Get the number of point in the file
    char output[64];
    vectorFile.getline(output, 32);
    outfileVec << output << std::endl;
    int Val = GetIntFromString(output);
    vectorFile.getline(output, 32);
    outfileVec << output << std::endl;
    vectorFile.getline(output, 32);
    outfileVec << output << std::endl;

    float X, Y, Z;
    for( int Line = 0; Line < Val; Line++ )
      {
      vectorFile.getline(output, 64);
      sscanf(output, "%f %f %f", &X, &Y, &Z);
      X = -1 * X;
      Y = -1 * Y;
      Z = -1 * Z;
      outfileVec << X << " " << Y << " " << Z << std::endl;
      }
    outfileVec.close();

    }
  else if( magdirOn )
    {

    // Read the meta file
    MeshConverterType::Pointer converter = MeshConverterType::New();
    if( debug )
      {
      std::cout << "Reading input mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer inputMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer() );
    MeshType::Pointer   inputMesh = inputMeshSO->GetMesh();
    int                 NbOfPoint = inputMesh->GetNumberOfPoints();

    // Convert the meta data into the VTK Poly Data format
    if( debug )
      {
      std::cout << "Computing normals" << std::endl;
      }
    TriangleMeshConverterType::Pointer Triangleconverter = TriangleMeshConverterType::New();
    TriangleMeshSOType::Pointer inputTriangleMeshSO =
     dynamic_cast<TriangleMeshSOType *>(Triangleconverter->ReadMeta(inputFilename).GetPointer());
    TriangleMeshType::Pointer   inputTriangleMesh = inputTriangleMeshSO->GetMesh();

    itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
    convertMeshToVTK->SetInput(inputTriangleMesh);

    // compute the normal
    vtkPolyData *       vtkMesh = convertMeshToVTK->GetOutput();
    vtkPolyDataNormals *MeshNormals = vtkPolyDataNormals::New();

    MeshNormals->SetComputePointNormals(1);
    MeshNormals->SetComputeCellNormals(0);
    MeshNormals->SetSplitting(0);
    #if VTK_MAJOR_VERSION > 5
    MeshNormals->SetInputData(vtkMesh);
    #else
    MeshNormals->SetInput(vtkMesh);
    #endif
    MeshNormals->Update();
    vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    vtkMeshNormals->Update();
    #endif
    vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
    vtkDataArray * ArrayNormal = NormalPoints->GetNormals();
    double *       MagList = new double[NbOfPoint];

    // Read the vector field file
    char out[64];
    if( debug )
      {
      std::cout << "Reading vector field " << VectFile << std::endl;
      }
    // Open the file
    std::ifstream avefile;
    avefile.open(VectFile);
    // Get the number of point in the file
    avefile.getline(out, 32);
    int Val = GetIntFromString(out);
    if( Val != NbOfPoint )
      {
      std::cerr << "The file must have the same number of point as the meta file" << std::endl;
      exit(-1);
      }
    else
      {
      double *Normal = new double[3];

      // Read the two next lines
      avefile.getline(out, 64);
      avefile.getline(out, 64);
      // Then get the numbers
      for( int Line = 0; Line < NbOfPoint; Line++ )
        {
        float X, Y, Z;
        avefile.getline(out, 64);
        sscanf(out, "%f %f %f", &X, &Y, &Z);
        double magnitude = sqrt( (double)X * (double)X + (double)Y * (double)Y + (double)Z * (double)Z);
        MagList[Line] = magnitude;
        Normal = ArrayNormal->GetTuple3(Line);
        double DotProduct = X * Normal[0] + Y * Normal[1] + Z * Normal[2];
        if( DotProduct < 0 )
          {
          MagList[Line] = -1 * MagList[Line];
          }
        }
      delete [] Normal;
      }
    if( debug )
      {
      std::cout << "Write the magnitude file" << std::endl;
      }
    std::ofstream outfileNor;
    outfileNor.open( outputFilename );
    // print the header
    outfileNor << "NUMBER_OF_POINTS=" << NbOfPoint << std::endl;
    outfileNor << "DIMENSION=" << 1 << std::endl;
    outfileNor << "TYPE=Scalar" << std::endl;
    for( int OutLine = 0; OutLine < NbOfPoint; OutLine++ )
      {
      outfileNor << MagList[OutLine] << std::endl;
      }
    outfileNor.close();

    delete [] MagList;

    }
  else if( magNormDirOn )
    {

    // Read the meta file
    MeshConverterType::Pointer converter = MeshConverterType::New();
    if( debug )
      {
      std::cout << "Reading input mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer inputMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer   inputMesh = inputMeshSO->GetMesh();
    int                 NbOfPoint = inputMesh->GetNumberOfPoints();

    // Convert the meta data into the VTK Poly Data format
    if( debug )
      {
      std::cout << "Computing normals" << std::endl;
      }
    TriangleMeshConverterType::Pointer Triangleconverter = TriangleMeshConverterType::New();
    TriangleMeshSOType::Pointer inputTriangleMeshSO =
      dynamic_cast<TriangleMeshSOType *>(Triangleconverter->ReadMeta(inputFilename).GetPointer());
    TriangleMeshType::Pointer   inputTriangleMesh = inputTriangleMeshSO->GetMesh();

    itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
    convertMeshToVTK->SetInput(inputTriangleMesh);

    // compute the normal
    vtkPolyData *       vtkMesh = convertMeshToVTK->GetOutput();
    vtkPolyDataNormals *MeshNormals = vtkPolyDataNormals::New();

    MeshNormals->SetComputePointNormals(1);
    MeshNormals->SetComputeCellNormals(0);
    MeshNormals->SetSplitting(0);
    #if VTK_MAJOR_VERSION > 5
    MeshNormals->SetInputData(vtkMesh);
    #else
    MeshNormals->SetInput(vtkMesh);
    #endif
    MeshNormals->Update();
    vtkPolyData * vtkMeshNormals = MeshNormals->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    vtkMeshNormals->Update();
    #endif

    vtkPointData * NormalPoints = vtkMeshNormals->GetPointData();
    vtkDataArray * ArrayNormal = NormalPoints->GetNormals();
    double *       MagList = new double[NbOfPoint];

    // Read the vector field file
    char out[64];
    if( debug )
      {
      std::cout << "Reading vector field " << VectFile << std::endl;
      }
    // Open the file
    std::ifstream avefile;
    avefile.open(VectFile);
    // Get the number of point in the file
    avefile.getline(out, 32);
    int         Val = GetIntFromString(out);
    double      Mean = 0.0;
    if( Val != NbOfPoint )
      {
      std::cerr << "The file must have the same number of point as the meta file" << std::endl;
      exit(-1);
      }
    else
      {
      double *Normal = new double[3];

      // Read the next two lines
      avefile.getline(out, 64);
      avefile.getline(out, 64);
      // Then get the numbers
      for( int Line = 0; Line < NbOfPoint; Line++ )
        {
        float X, Y, Z;
        avefile.getline(out, 64);
        sscanf(out, "%f %f %f", &X, &Y, &Z);
        Normal = ArrayNormal->GetTuple3(Line);
        double DotProduct = (double)Normal[0] * (double)X + (double)Normal[1] * (double)Y + (double)Normal[2]
          * (double)Z;
        double XProj = (double)Normal[0] * DotProduct;
        double YProj = (double)Normal[1] * DotProduct;
        double ZProj = (double)Normal[2] * DotProduct;

        double magnitude = sqrt(XProj * XProj + YProj * YProj + ZProj * ZProj);
        MagList[Line] = magnitude;
        if( DotProduct < 0 )
          {
          MagList[Line] = -1 * MagList[Line];
          }
        Mean += MagList[Line];
        }
      Mean /= (double)NbOfPoint;
      delete [] Normal;
      }

    std::ofstream outfileNor, outfileNorCenter;
    if( debug )
      {
      std::cout << "Write the signed magnitude file" << std::endl;
      }
    outfileNor.open( outputFilename );
    // print the header
    outfileNor << "NUMBER_OF_POINTS=" << NbOfPoint << std::endl;
    outfileNor << "DIMENSION=" << 1 << std::endl;
    outfileNor << "TYPE=Scalar" << std::endl;
    for( int OutLine = 0; OutLine < NbOfPoint; OutLine++ )
      {
      outfileNor << MagList[OutLine] << std::endl;
      }
    outfileNor.close();

    outputFilename2 = std::string(outputFilename);
    DotLoc = outputFilename2.find(".txt", 0);
    if( DotLoc != -1 )
      {
      outputFilename2.insert(DotLoc, "_centered");
      }
    else
      {
      outputFilename2 = std::string(outputFilename) + std::string(".centered");
      }

    if( debug )
      {
      std::cout << "Write the centered signed magnitude file" << std::endl;
      }
    outfileNorCenter.open( outputFilename2.c_str() );
    // print the header
    outfileNorCenter << "NUMBER_OF_POINTS=" << NbOfPoint << std::endl;
    outfileNorCenter << "DIMENSION=" << 1 << std::endl;
    outfileNorCenter << "TYPE=Scalar" << std::endl;
    for( int OutLine = 0; OutLine < NbOfPoint; OutLine++ )
      {
      outfileNorCenter << MagList[OutLine] - Mean << std::endl;
      }
    outfileNorCenter.close();

    delete [] MagList;

    }
  else if( applyvecOn )
    {

    MeshConverterType::Pointer converter = MeshConverterType::New();
    if( debug )
      {
      std::cout << "Reading input mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer         inputMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer() );
    MeshType::Pointer           inputMesh = inputMeshSO->GetMesh();
    TriangleMeshConverterType::Pointer Triangleconverter = TriangleMeshConverterType::New();
    TriangleMeshSOType::Pointer outputTriangleMeshSO =
      dynamic_cast<TriangleMeshSOType *>(Triangleconverter->ReadMeta(inputFilename).GetPointer());
    TriangleMeshType::Pointer   outputTriangleMesh = outputTriangleMeshSO->GetMesh();

    char          output[64];
    std::ifstream vectorFile;
    vectorFile.open(VectFile);
    if( debug )
      {
      std::cout << "Reading the vector file " << VectFile << std::endl;
      }
    vectorFile.getline(output, 64);
    unsigned int Val = GetIntFromString(output);

    // make sure the two meshes have the same number of verts
    if( inputMesh->GetNumberOfPoints() != Val )
      {
      std::cout << "Mesh and vector field must have the same number of points" << std::endl;
      std::cout << "The mesh has " << inputMesh->GetNumberOfPoints() << " points and the vector field has " << Val
                << " points." << std::endl;
      exit(-1);
      }
    // Skipping the second and third lines of the file
    vectorFile.getline(output, 32);
    vectorFile.getline(output, 32);

    float              X, Y, Z;
    std::vector<float> Xval;
    std::vector<float> Yval;
    std::vector<float> Zval;
    for( unsigned int Line = 0; Line < Val; Line++ )
      {
      vectorFile.getline(output, 64);
      sscanf(output, "%f %f %f", &X, &Y, &Z);
      Xval.push_back(X);
      Yval.push_back(Y);
      Zval.push_back(Z);
      }

    PointType *        point1 = new PointType;
    PointTriangleType *Tripoint = new PointTriangleType;
    for( unsigned int i = 0; i < inputMesh->GetNumberOfPoints(); i++ )
      {
      if( inputMesh->GetPoint( i, point1 ) && &Xval[i] != NULL )
        {
        if( debug )
          {
          std::cout << "Original Point[" << i << "]= " << point1->GetElement(0) << " | " << point1->GetElement(1)
                    << " | " << point1->GetElement(2) << std::endl;
          }
        Tripoint->SetElement(0, point1->GetElement(0) + Xval[i]);
        Tripoint->SetElement(1, point1->GetElement(1) + Yval[i]);
        Tripoint->SetElement(2, point1->GetElement(2) + Zval[i]);
        if( debug )
          {
          std::cout << "New Point[" << i << "]= " << Tripoint->GetElement(0) << " | " << Tripoint->GetElement(1)
                    << " | " << Tripoint->GetElement(2) << std::endl;
          }
        outputTriangleMesh->SetPoint(i, *Tripoint);
        }
      else
        {
        return -1;
        }
      }

    // Convert the outputTriangleMesh to a vtk poly data.
    itkMeshTovtkPolyData *convertMeshToVTK = new itkMeshTovtkPolyData();
    convertMeshToVTK->SetInput(outputTriangleMesh);
    vtkPolyData * vtkPolyData = convertMeshToVTK->GetOutput();

    // Convert the vtk poly data to itk mesh data structure
    vtkPolyDataToitkMesh *vtkItkConverter = new vtkPolyDataToitkMesh();
    vtkItkConverter->SetInput(vtkPolyData);

    // Convert the itk mesh data in Spatial Object mesh
    // write out the itk spatial object meta mesh file
    TriangleMeshSOType::Pointer meshSO = TriangleMeshSOType::New();
    meshSO->SetMesh(vtkItkConverter->GetOutput() );
    Triangleconverter->WriteMeta(meshSO, outputFilename);

    }
  else if( MC2OriginOn ) // cchou
    {
    MeshConverterType::Pointer converter = MeshConverterType::New();
    MeshSOType::Pointer        SOMesh =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer                surfaceMesh = SOMesh->GetMesh();
    MeshType::PointsContainerPointer points = surfaceMesh->GetPoints();

    // Sum up Original Points
    double sum[3];
    sum[0] = 0;
    sum[1] = 0;
    sum[2] = 0;

    for( unsigned int pointID = 0; pointID < points->Size(); pointID++ )
      {
      PointType curPoint = points->GetElement(pointID);
      for( unsigned int dim = 0; dim < 3; dim++ )
        {
        sum[dim] += curPoint[dim];
        }
      }
    // Calculate the Center of Mass
    double MC[3];
    for( unsigned int dim = 0; dim < 3; dim++ )
      {
      MC[dim] = (double) sum[dim] / (points->Size() + 1);
      }
   // if (debug) 
    //  {
        std::cout << "Mesh " << inputFilename << " MC " << MC[0] << " " << MC[1] << " " << MC[2] << std::endl;
    //  }
    // Create a New Point Set "newpts" with the Shifted Values
    MeshType::PointsContainerPointer newpts = MeshType::PointsContainer::New();
    for( unsigned int pointID = 0; pointID < points->Size(); pointID++ )
      {
      PointType curPoint = points->GetElement(pointID);
      PointType sftPoint;
      for( unsigned int dim = 0; dim < 3; dim++ )
        {
        sftPoint[dim] = curPoint[dim] - MC[dim];
        }
      newpts->InsertElement(pointID, sftPoint);
      }
    // Set the Shifted Points back to "points"
    points = newpts;
    surfaceMesh->SetPoints(points);
    SOMesh->SetMesh(surfaceMesh);
    MeshWriterType::Pointer writer = MeshWriterType::New();
    writer->SetInput(SOMesh);
    writer->SetFileName(outputFilename);
    writer->Update();

    }
  else if( testFillOn )  // bp2009
    {
    int NBorders;

    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(inputFilename);
    meshin->Update();

    if( debug )
      {
      std::cout << "Number of points before cleaning...  " << meshin->GetOutput()->GetNumberOfPoints() << std::endl;
      }

    vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    meshinC->SetInputData(meshin->GetOutput() );
    #else
    meshinC->SetInput(meshin->GetOutput() );
    #endif
    meshinC->Update();

    if( debug )
      {
      std::cout << "Number of points before patching...  " << meshinC->GetOutput()->GetNumberOfPoints() << std::endl;
      }

    // detect borders
    vtkFeatureEdges *boundaryEdges = vtkFeatureEdges::New();
    #if VTK_MAJOR_VERSION > 5
    boundaryEdges->SetInputData(meshinC->GetOutput() );
    #else
    boundaryEdges->SetInput(meshinC->GetOutput() );
    #endif
    boundaryEdges->BoundaryEdgesOn();       // only border points
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
    boundaryEdges->Update();

    // "NBorders" gives us the total number of border edges
    NBorders = (boundaryEdges->GetOutput()->GetNumberOfLines() );

    if( debug )
      {
      std::cout << "Number of Border Edges...  " << NBorders << std::endl;
      }

    vtkStripper *patchSkel = vtkStripper::New();
    #if VTK_MAJOR_VERSION > 5
    patchSkel->SetInputData(boundaryEdges->GetOutput() );
    #else
    patchSkel->SetInput(boundaryEdges->GetOutput() );
    #endif

    vtkPolyData *      patchPoly = vtkPolyData::New();
    vtkTriangleFilter *patchTri = vtkTriangleFilter::New();

    // ORIGINAL
    patchPoly->Initialize();
    #if VTK_MAJOR_VERSION > 5
    patchTri->SetInputData(patchPoly);
    #else
    patchTri->SetInput(patchPoly);
    #endif

    vtkAppendPolyData *meshAppend = vtkAppendPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    meshAppend->AddInputConnection(patchTri->GetOutputPort() );
    meshAppend->AddInputConnection(meshinC->GetOutputPort() );
    #else
    meshAppend->AddInput(patchTri->GetOutput() );
    meshAppend->AddInput(meshinC->GetOutput() );
    #endif
    vtkCleanPolyData *poly = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    poly->SetInputData(meshAppend->GetOutput() );
    #else
    poly->SetInput(meshAppend->GetOutput() );
    #endif
    // Creating the patch!
    patchSkel->Update();
    patchPoly->Initialize();
    patchPoly->SetPoints(patchSkel->GetOutput()->GetPoints() );
    patchPoly->SetPolys(patchSkel->GetOutput()->GetLines() );
    #if VTK_MAJOR_VERSION <= 5
    patchPoly->Update();
    #endif

    // Triangulate the patches
    patchTri->Update();

    if( debug )
      {
      std::cout << "Number of points in the patch...  " << patchTri->GetOutput()->GetNumberOfPoints() << std::endl;
      }

    // Add the patch in the new mesh
    meshAppend->Update();

    // Update the new mesh, including the patches
    poly->Update();

    vtkCleanPolyData *polyC = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    polyC->SetInputData(poly->GetOutput() );
    #else
    polyC->SetInput(poly->GetOutput() );
    #endif
    polyC->Update();

    // Verify if there are still open borders in the new patched mesh
    #if VTK_MAJOR_VERSION > 5
    boundaryEdges->SetInputData(polyC->GetOutput() );
    #else
    boundaryEdges->SetInput(polyC->GetOutput() );
    #endif
    // boundaryEdges->SetInput(normals->GetOutput());
    boundaryEdges->Update();

    NBorders = ( (boundaryEdges->GetOutput() )->GetNumberOfLines() );
    if( debug )
      {
      if( NBorders == 0 )
        {
        std::cout << "MESH COMPLETELY CLOSED..." << std::endl;
        }
      else
        {
        std::cout << NBorders << " open borders still remain... " << std::endl;
        }
      std::cout << "Number of points after patching...  " << polyC->GetOutput()->GetNumberOfPoints() << std::endl;
      }

    // project week update Arnaud Gelas
    vtkPolyDataNormals* normals = vtkPolyDataNormals::New();
    #if VTK_MAJOR_VERSION > 5
    normals->SetInputData( polyC->GetOutput() );
    #else
    normals->SetInput( polyC->GetOutput() );
    #endif
    normals->SetConsistency( true );
    normals->SetFlipNormals( true );
    normals->GetFlipNormals();
    normals->Update();

    // Writing the new mesh, with the patched hole
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    // SurfaceWriter->SetInput(polyC->GetOutput());
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(normals->GetOutput() );
    #else
    SurfaceWriter->SetInput(normals->GetOutput() );
    #endif
    // SurfaceWriter->SetInput(patchTri->GetOutput());
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();
    if( debug )
      {
      std::cout << "Writing new mesh " << outputFilename << std::endl;
      }

    }
  else if( BordersOutOn ) // bp2009
    {
    int NBorders;
    int TOTALBorders;

    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(inputFilename);
    meshin->Update();

    vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    meshinC->SetInputData(meshin->GetOutput() );
    #else
    meshinC->SetInput(meshin->GetOutput() );
    #endif
    meshinC->Update();

    // detect borders
    vtkFeatureEdges *boundaryEdges = vtkFeatureEdges::New();
    #if VTK_MAJOR_VERSION > 5
    boundaryEdges->SetInputData(meshinC->GetOutput() );
    #else
    boundaryEdges->SetInput(meshinC->GetOutput() );
    #endif
    boundaryEdges->BoundaryEdgesOn();       // only border points
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
    boundaryEdges->Update();

    // "NBorders" gives us the total number of border edges
    NBorders = (boundaryEdges->GetOutput()->GetNumberOfLines() );
    TOTALBorders = NBorders;

    if( debug )
      {
      std::cout << "Number of Border Edges...  " << TOTALBorders << std::endl;
      }

    // Writing the new mesh, only detected borders
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(boundaryEdges->GetOutput() );
    #else
    SurfaceWriter->SetInput(boundaryEdges->GetOutput() );
    #endif
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();
    if( debug )
      {
      std::cout << "Writing new mesh " << outputFilename << std::endl;
      }

    }
  else if( IsOpenOn ) // bp2009
    {
    int NBorders;

    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(inputFilename);
    meshin->Update();

    vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    meshinC->SetInputData(meshin->GetOutput() );
    #else
    meshinC->SetInput(meshin->GetOutput() );
    #endif
    meshinC->Update();

    // detect borders
    vtkFeatureEdges *boundaryEdges = vtkFeatureEdges::New();
    #if VTK_MAJOR_VERSION > 5
    boundaryEdges->SetInputData(meshinC->GetOutput() );
    #else
    boundaryEdges->SetInput(meshinC->GetOutput() );
    #endif
    boundaryEdges->BoundaryEdgesOn();       // only border points
    boundaryEdges->FeatureEdgesOff();
    boundaryEdges->NonManifoldEdgesOff();
    boundaryEdges->ManifoldEdgesOff();
    boundaryEdges->Update();

    // "NBorders" gives us the total number of border edges
    NBorders = (boundaryEdges->GetOutput()->GetNumberOfLines() );

    if( NBorders == 0 )
      {
      std::cout << "0" << std::endl;
      }
    else
      {
      std::cout << "1" << std::endl;
      }

    }
  else if( avgOneKWMOn ) // bp2009
    {
    char     line[70];
    ifstream input;
    ofstream output;
    int      NPoints;
    float    value;
    float    sum = 0;
    char *   aux;
    float    avgOne = 0;

    // Start reading the Input
    input.open(inputFilename, ios::in);
    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    input.getline(line, 70, '\n');
    input.getline(line, 70, '\n');
    while( !input.getline(line, 70, '\n').eof() )
      {
      value = atof(line);
      sum += value;
      }

    input.close();
    // End reading the Input
    // Calculate the avg
    avgOne = sum / NPoints;

    // Start writing the Output
    output.open(outputFilename, ios::out);
    output << avgOne << endl;
    output.close();
    // End writing the Output
    }
  else if( minOneKWMOn ) // bp2009
    {
    char     line[70];
    ifstream input;
    ofstream output;
    int      NPoints;
    float    value;
    char *   aux;
    float    minOne = 9999999;

    // Start reading the Input
    input.open(inputFilename, ios::in);
    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    input.getline(line, 70, '\n');
    input.getline(line, 70, '\n');

    while( !input.getline(line, 70, '\n').eof() )
      {
      value = atof(line);
      if( value < minOne )
        {
        minOne = value;
        }
      }

    input.close();
    // End reading the Input

    // Start writing the Output
    output.open(outputFilename, ios::out);
    output << minOne << endl;
    output.close();
    // End writing the Output
    }
  else if( maxOneKWMOn ) // bp2009
    {
    char     line[70];
    ifstream input;
    ofstream output;
    int      NPoints;
    float    value;
    char *   aux;
    float    maxOne = 0;

    // Start reading the Input
    input.open(inputFilename, ios::in);
    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    input.getline(line, 70, '\n');
    input.getline(line, 70, '\n');

    while( !input.getline(line, 70, '\n').eof() )
      {
      value = atof(line);
      if( value > maxOne )
        {
        maxOne = value;
        }
      }

    input.close();
    // End reading the Input

    // Start writing the Output
    output.open(outputFilename, ios::out);
    output << maxOne << endl;
    output.close();
    // End writing the Output
    }
  else if( medianOneKWMOn ) // bp2009
    {
    char     line[70];
    ifstream input;
    ofstream output;
    int      NPoints;
    float    value;
    char *   aux;
    float    meanOne = 0;

    // Start reading the Input
    input.open(inputFilename, ios::in);
    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    input.getline(line, 70, '\n');
    input.getline(line, 70, '\n');

    std::vector<float> vector_qsort(NPoints);
    int   cont = 0;
    while( !input.getline(line, 70, '\n').eof() )
      {
      value = atof(line);
      vector_qsort[cont] = value;
      cont++;
      }

    input.close();
    // End reading the Input

    /* sort array using qsort functions */
    qsort(&(vector_qsort[0]), NPoints, sizeof(float), oldStyleCompare);

    int meanIndex = NPoints / 2;

    meanOne = vector_qsort[meanIndex];

    // Start writing the Output
    output.open(outputFilename, ios::out);
    output << meanOne << endl;
    output.close();
    // End writing the Output
    }
  else if( per1OneKWMOn ) // bp2009
    {
    char     line[70];
    ifstream input;
    ofstream output;
    int      NPoints;
    float    value;
    char *   aux;
    float    per1One = 0;

    // Start reading the Input
    input.open(inputFilename, ios::in);
    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    input.getline(line, 70, '\n');
    input.getline(line, 70, '\n');

    std::vector<float> vector_qsort(NPoints);
    int   cont = 0;
    while( !input.getline(line, 70, '\n').eof() )
      {
      value = atof(line);
      vector_qsort[cont] = value;
      cont++;
      }

    input.close();
    // End reading the Input

    /* sort array using qsort functions */
    qsort(&(vector_qsort[0]), NPoints, sizeof(float), oldStyleCompare);

    int per1Index = (1 * NPoints) / 100;

// std::cout << per1I/ndex << std::endl;
    per1One = vector_qsort[per1Index];

    // Start writing the Output
    output.open(outputFilename, ios::out);
    output << per1One << endl;
    output.close();
    // End writing the Output
    }
  else if( per99OneKWMOn ) // bp2009
    {
    char     line[70];
    ifstream input;
    ofstream output;
    int      NPoints;
    float    value;
    char *   aux;
    float    per99One = 0;

    // Start reading the Input
    input.open(inputFilename, ios::in);
    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    input.getline(line, 70, '\n');
    input.getline(line, 70, '\n');

    std::vector<float> vector_qsort(NPoints);
    int   cont = 0;
    while( !input.getline(line, 70, '\n').eof() )
      {
      value = atof(line);
      vector_qsort[cont] = value;
      cont++;
      }

    input.close();
    // End reading the Input

    /* sort array using qsort functions */
    qsort(&(vector_qsort[0]), NPoints, sizeof(float), oldStyleCompare);

    int per99Index = (99 * NPoints) / 100;

// std::cout << per1I/ndex << std::endl;
    per99One = vector_qsort[per99Index];

    // Start writing the Output
    output.open(outputFilename, ios::out);
    output << per99One << endl;
    output.close();
    // End writing the Output
    }
  else if( cleanMeshOn ) // bp2009
    {
    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(inputFilename);
    meshin->Update();

    vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    meshinC->SetInputData(meshin->GetOutput() );
    #else
    meshinC->SetInput(meshin->GetOutput() );
    #endif
    meshinC->Update();

    // Writing the new mesh, with not degenerated triangles
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(meshinC->GetOutput() );
    #else
    SurfaceWriter->SetInput(meshinC->GetOutput() );
    #endif
    // SurfaceWriter->SetInput(boundaryEdges->GetOutput());
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();
    if( debug )
      {
      std::cout << "Writing new mesh " << outputFilename << std::endl;
      }

    }
  else if( smoothMeshOn ) // bp2009
    {
    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(inputFilename);
    meshin->Update();

    vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
    #if VTK_MAJOR_VERSION > 5
    meshinC->SetInputData(meshin->GetOutput() );
    #else
    meshinC->SetInput(meshin->GetOutput() );
    #endif
    meshinC->Update();

    vtkSmoothPolyDataFilter *smoother = vtkSmoothPolyDataFilter::New();
    smoother->SetInputConnection(meshinC->GetOutputPort() );
    smoother->SetNumberOfIterations(smoothIterationNb);

    // Writing the new mesh, with not degenerated triangles
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(smoother->GetOutput() );
    #else
    SurfaceWriter->SetInput(smoother->GetOutput() );
    #endif
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();

    if( debug )
      {
      std::cout << "Writing new mesh " << outputFilename << std::endl;
      }
    }
  else if( statsROIOn ) // bp2009
    {
    char     line[70];
    char     line2[70];
    ifstream inputROI, inputFile;
    ofstream output;
    int      NPoints, Dim;
    float    valueROI, valueFile;
    char *   aux;

    // Start reading the InputMASK
    if( debug )
      {
      std::cout << "Reading the mask " << ROIFileName << std::endl;
      }
    inputROI.open(ROIFileName, ios::in);
    inputFile.open(inputFilename, ios::in);
    output.open(outputFilename, ios::out);
    inputROI.getline(line, 70, '\n');
    inputFile.getline(line2, 70, '\n');
    output << line2 << std::endl;
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atoi(aux);
    inputROI.getline(line, 70, '\n');
    inputFile.getline(line2, 70, '\n');
    output << line2 << std::endl;
    aux = strtok(line2, " = ");
    aux = strtok(NULL, " = ");
    Dim = atoi(aux);
    // std::cout << "Dimensions ..." << Dim << std::endl;
    inputROI.getline(line, 70, '\n');
    inputFile.getline(line2, 70, '\n');
    output << line2 << std::endl;

    std::vector<float>  vectorOut(NPoints);
    unsigned int pointNumber = 0;
    int          zeroes = 0;
    while( !inputROI.getline(line, 70, '\n').eof() )
      {
      valueROI = atof(line);
      inputFile.getline(line2, 70, '\n');
      valueFile = atof(line2);
      if( valueROI == 1 )
        {
        if( Dim == 1 )
          {
          vectorOut[pointNumber] = valueFile;
          output <<  valueFile << std::endl;
          }
        else
          {
          output <<  line2 << std::endl;
          }
        }
      else
        {
        if( Dim == 1 )
          {
          vectorOut[pointNumber] = valueROI;
          output <<  valueROI << std::endl;
          }
        else
          {
          output << "0.00 0.00 0.00" << std::endl;
          }
        zeroes++;
        }
      pointNumber++;
      }

    if( debug )
      {
      std::cout << "Mask ROI with a total of "
                << (pointNumber - zeroes) << " correspondent points selected !!" << std::endl;
      }

    inputROI.close();
    inputFile.close();
    output.close();
    // End reading the InputMASK
    }
  else if( filterNormalsOn ) // bp2009
    {

    std::cout << "Arguments " << " " << argv[0] << " " << inputFilename << " " << argv[2] << " " << argv[3] << std::endl;

    std::cout << "Reading mesh " << " " << argv[2] << std::endl;
    vtkPolyDataReader *meshin = vtkPolyDataReader::New();
    meshin->SetFileName(argv[2]);
    meshin->Update();
    vtkPolyData *polydata = meshin->GetOutput();
    // unsigned int nPoints = polydata->GetNumberOfPoints();
    // std::cout << nPoints << std::endl ;

    // generate normals
    vtkPolyDataNormals* normals = vtkPolyDataNormals::New();
    #if VTK_MAJOR_VERSION > 5
    normals->SetInputData( polydata );
    #else
    normals->SetInput( polydata );
    #endif
    // normals->FlipNormalsOn();
    normals->AutoOrientNormalsOn();
    normals->SetConsistency(1);
    normals->SplittingOff();
    normals->Update();
    polydata = normals->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    polydata->Update();
    #endif

/*
  //optional settings
  normalGenerator->SetFeatureAngle(0.1);
  normalGenerator->SetSplitting(1);
  normalGenerator->SetConsistency(0);
  normalGenerator->SetAutoOrientNormals(0);
  normalGenerator->SetComputePointNormals(1);
  normalGenerator->SetComputeCellNormals(0);
  normalGenerator->SetFlipNormals(0);
  normalGenerator->SetNonManifoldTraversal(1);
*/

     // Writing the new mesh, with the patched hole
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    // SurfaceWriter->SetInput(polyC->GetOutput());
    // SurfaceWriter->SetInput(normals->GetOutput());
    // SurfaceWriter->SetInput(patchTri->GetOutput());
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(polydata);
    #else
    SurfaceWriter->SetInput(polydata);
    #endif
    SurfaceWriter->SetFileName(argv[3]);
    SurfaceWriter->Update();
    std::cout << "Writing new mesh " << argv[3] << std::endl;

    }
  else if( KWMtoPolyDataOn ) // bp2009
    {
    vtkPolyDataReader *polyIn = vtkPolyDataReader::New();
    polyIn->SetFileName(inputFilename);
    polyIn->Update();
    vtkPolyData* polydataAtt = polyIn->GetOutput();

    if( debug )
      {
      std::cout << "Input Mesh read ...  " << inputFilename << std::endl;
      }

    // *** START PARSING KWMeshVisu file
    char     line[70];
    ifstream input;
    int      NPoints, NDimension;
    char *   aux;
    input.open(files[0], ios::in);

    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NPoints = atof(aux);

    input.getline(line, 70, '\n');
    aux = strtok(line, " = ");
    aux = strtok(NULL, " = ");
    NDimension = atof(aux);

    input.getline(line, 70, '\n'); // read type line

    vtkFloatArray *scalars = vtkFloatArray::New();
    scalars->SetNumberOfComponents(NDimension);
    scalars->SetName(files[1]);

    int cont = 0;
    for( int i = 0; i < NPoints; i++ )
      {
          input.getline(line, 70, '\n');
	  float value=0;
	  if (NDimension == 1)
          { 
            value = atof(line);
	    scalars->InsertNextValue(value);
          }
          else
	  {
            std::string proc_string = line;
	    std::istringstream iss(proc_string);

	    do{
              std::string sub;
              iss >> sub;
 	      if (sub!= "")
              {
	        value = atof(sub.c_str());
                scalars->InsertNextValue(value);
	      }	
	    }while (iss);
          } 
	  cont++;
      }

    input.close();
    // End reading the Input
    // *** END PARSING KWMeshVisu file

    polydataAtt->GetPointData()->AddArray(scalars);
        if( debug )
      {
      std::cout << "Scalar map added to the mesh" << std::endl;
      }

    // Writing the new mesh
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(polydataAtt);
    #else
    SurfaceWriter->SetInput(polydataAtt);
    #endif
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();
    if( debug )
      {
      std::cout << "Writing new mesh " << outputFilename << std::endl;
      }

    scalars->Delete();
    }
  else if( surfaceAreaOn )
    {
    // Reading the input mesh
    MeshConverterType::Pointer converter = MeshConverterType::New();
    if( debug )
      {
      std::cout << "Reading input mesh " << inputFilename << std::endl;
      }
    MeshSOType::Pointer inputMeshSO =
      dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer());
    MeshType::Pointer                inputMesh = inputMeshSO->GetMesh();
    MeshType::PointsContainerPointer points = inputMesh->GetPoints();

    // Reading the input attribute file
    if( debug )
      {
      std::cout << "Reading attribute file " << inputFilename << std::endl;
      }
    std::ifstream    Infile;
    char             Line[40];
    int              CurrentAttribute;
    std::vector<int> v_VerticesAttributes; // vector with attribute values
    std::vector<int> v_Attributes;         // Attributes vector
    Infile.open(AttributeFileName);
    while( strncmp(Line, "NUMBER_OF_POINTS =", 18) && strncmp(Line, "NUMBER_OF_POINTS=", 17) )
      {
      Infile.getline(Line, 40);
      }

    unsigned int NbVertices = atoi(strrchr(Line, '=') + 1);
    if( inputMesh->GetNumberOfPoints() != NbVertices )
      {
      std::cerr << "Input mesh and attribute file must have the same number of vertices!" << std::endl;
      exit(-1);
      }
    Infile.getline( Line, 40);
    Infile.getline( Line, 40);
    for( unsigned int i = 0; i < NbVertices; i++ )
      {
      Infile >> CurrentAttribute;
      v_VerticesAttributes.push_back(CurrentAttribute);
      bool NewAttribute = true;
      for( unsigned int j = 0; j < v_Attributes.size(); j++ )
        {
        if( CurrentAttribute == v_Attributes[j] )
          {
          NewAttribute = false;
          }
        }
      if( NewAttribute )
        {
        v_Attributes.push_back(CurrentAttribute);
        }
      }
    Infile.close();

    // Sort attributes vector into ascending order
    sort(v_Attributes.begin(), v_Attributes.end() );

    // Computing surface area
    typedef CellType::PointIdIterator PointIdIterator;
    CellIterator        cellIterator = inputMesh->GetCells()->Begin();
    CellIterator        cellEnd = inputMesh->GetCells()->End();
    PointType           curPoint;
    std::vector<double> v_Area(*(std::max_element(v_Attributes.begin(), v_Attributes.end() ) ) + 1);

    while( cellIterator != cellEnd )
      {
      CellType *                cell = cellIterator.Value();
      TriangleType *            line = dynamic_cast<TriangleType *>(cell);
      LineType::PointIdIterator pit = line->PointIdsBegin();
      int                       pIndex1, pIndex2, pIndex3;
      PointType                 p1, p2, p3;
      itk::Vector<float, 3>     v1, v2;

      pIndex1 = *pit;
      p1 = points->GetElement(pIndex1);
      ++pit;
      pIndex2 = *pit;
      p2 = points->GetElement(pIndex2);
      ++pit;
      pIndex3 = *pit;
      p3 = points->GetElement(pIndex3);

      // Computing vectors of the current triangle
      v1 = p2 - p1;
      v2 = p3 - p1;

      // Computing triangle area
      double area =
        sqrt(pow(v1[1] * v2[2] - v2[1] * v1[2],
                 2) + pow(v1[2] * v2[0] - v2[2] * v1[0], 2) + pow(v1[0] * v2[1] - v2[0] * v1[1], 2) ) / 2.0;
      if( (v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex2]) &&
          (v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex3]) )
        {
        v_Area[v_VerticesAttributes[pIndex1]] += area;
        }
      else if( (v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex2]) &&
               (v_VerticesAttributes[pIndex1] != v_VerticesAttributes[pIndex3]) )
        {
        v_Area[v_VerticesAttributes[pIndex1]] += 2.0 * area / 3.0;
        v_Area[v_VerticesAttributes[pIndex3]] += area / 3.0;
        }
      else if( (v_VerticesAttributes[pIndex1] != v_VerticesAttributes[pIndex2]) &&
               (v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex3]) )
        {
        v_Area[v_VerticesAttributes[pIndex1]] += 2.0 * area / 3.0;
        v_Area[v_VerticesAttributes[pIndex2]] += area / 3.0;
        }
      else if( (v_VerticesAttributes[pIndex1] != v_VerticesAttributes[pIndex2]) &&
               (v_VerticesAttributes[pIndex2] == v_VerticesAttributes[pIndex3]) )
        {
        v_Area[v_VerticesAttributes[pIndex1]] += area / 3.0;
        v_Area[v_VerticesAttributes[pIndex2]] += 2.0 * area / 3.0;
        }
      else
        {
        v_Area[v_VerticesAttributes[pIndex1]] += area / 3.0;
        v_Area[v_VerticesAttributes[pIndex2]] += area / 3.0;
        v_Area[v_VerticesAttributes[pIndex3]] += area / 3.0;
        }
      ++cellIterator;
      }

    if( debug )
      {
      std::cout << "Writing output file..." << std::endl;
      }
    std::ofstream OutputFile(outputFilename, ios::out);
    if( !OutputFile )
      {
      std::cerr << "Error: cannot open " << outputFilename << std::endl;
      exit(-1);
      }
    OutputFile << "##########################################################################" << std::endl;
    OutputFile << "#File format :" << std::endl;
    OutputFile << "#LABEL \t AREA" << std::endl;
    OutputFile << "#Fields :" << std::endl;
    OutputFile << "#\t LABEL         Label number" << std::endl;
    OutputFile << "#\t AREA          Surface area in square mm" << std::endl;
    OutputFile << "##########################################################################" << std::endl
               << std::endl;
    for( unsigned int i = 0; i < v_Attributes.size(); i++ )
      {
      OutputFile << v_Attributes[i] << "\t" << v_Area[v_Attributes[i]] << std::endl;
      }
    OutputFile.close();

    }
  else if( varianceOn )
    {
    float MeanVariance = 0.0;
    if( debug )
      {
      std::cout << "NbFiles: " << nbfile << std::endl;
      }
    for( int file = 0; file < nbfile; file++ )
      {
      // Reading attribute file
      std::ifstream      Infile;
      char               Line[40];
      float              CurrentAttribute, Mean, Variance;
      std::vector<float> v_Attribute;
      if( debug )
        {
        std::cout << "Reading attribute file " << VarFiles[file] << std::endl;
        }
      Infile.open(VarFiles[file].c_str() );
      while( strncmp(Line, "NUMBER_OF_POINTS =", 18) && strncmp(Line, "NUMBER_OF_POINTS=", 17) )
        {
        Infile.getline(Line, 40);
        }

      unsigned int NbVertices = atoi(strrchr(Line, '=') + 1);
      Infile.getline( Line, 40);
      Infile.getline( Line, 40);

      Mean = 0.0;
      for( unsigned int i = 0; i < NbVertices; i++ )
        {
        Infile >> CurrentAttribute;
        v_Attribute.push_back(CurrentAttribute);
        Mean += CurrentAttribute;
        }
      Infile.close();
      Mean /= NbVertices;
      if( debug )
        {
        std::cout << "\tMean: " << Mean << std::endl;
        }

      Variance = 0.0;
      for( unsigned int i = 0; i < NbVertices; i++ )
        {
        Variance += pow(v_Attribute[i] - Mean, 2);
        }
      Variance /= NbVertices;
      if( debug )
        {
        std::cout << "\tVariance: " << Variance << std::endl;
        }
      MeanVariance += Variance;
      }
    MeanVariance /= nbfile;
    if( debug )
      {
      std::cout << "Mean variance: " << MeanVariance << std::endl;
      }

    if( debug )
      {
      std::cout << "Writing output file..." << std::endl;
      }
    std::ofstream OutputFile(outputFilename, ios::out);
    if( !OutputFile )
      {
      std::cerr << "Error: cannot open " << outputFilename << std::endl;
      exit(-1);
      }
    std::cout << MeanVariance << std::endl;
    OutputFile << MeanVariance << std::endl;
    OutputFile.close();
    }
  else if( GetCurvaturesOn )   // bp2010
    {
    std::cout << "Arguments " << " " << argv[0] << " " << inputFilename << " " << argv[2] << " " << argv[3] << " "
              << argv[4] << std::endl;

    // Read PolyData info
    vtkPolyDataReader *polyIn = vtkPolyDataReader::New();
    polyIn->SetFileName(inputFilename); polyIn->Update();
    vtkPolyData* polydata = polyIn->GetOutput();

    // Calculate Curvedness and Shape Index
    vtkCurvatures* curveMax = vtkCurvatures::New();
    #if VTK_MAJOR_VERSION > 5
    curveMax->SetInputData(polydata);
    #else
    curveMax->SetInput(polydata);
    #endif
    curveMax->SetCurvatureTypeToMaximum();
    curveMax->Update();
    vtkPolyData* polydataCurvMax = curveMax->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    polydataCurvMax->Update();
    #endif
    vtkCurvatures* curveMin = vtkCurvatures::New();
    #if VTK_MAJOR_VERSION > 5
    curveMin->SetInputData(polydata);
    #else
    curveMin->SetInput(polydata);
    #endif
    curveMin->SetCurvatureTypeToMinimum();
    curveMin->Update();
    vtkPolyData* polydataCurvMin = curveMin->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    polydataCurvMin->Update();
    #endif
    // Writing curvatures out
    std::ofstream curvedness( argv[3] );
    std::ofstream shapeIndex( argv[4] );
    std::ofstream gauss( argv[5] );
    std::ofstream mean( argv[6] );

    unsigned int    nPoints = polydataCurvMax->GetNumberOfPoints();
    vtkDoubleArray *ArrayCurvMax = vtkDoubleArray::SafeDownCast(polydataCurvMax->GetPointData()->GetScalars() );
    vtkDoubleArray *ArrayCurvMin = vtkDoubleArray::SafeDownCast(polydataCurvMin->GetPointData()->GetScalars() );

    curvedness << "NUMBER_OF_POINTS=" << nPoints << std::endl;
    curvedness << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl;

    shapeIndex << "NUMBER_OF_POINTS=" << nPoints << std::endl;
    shapeIndex << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl;

    gauss << "NUMBER_OF_POINTS=" << nPoints << std::endl;
    gauss << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl;

    mean << "NUMBER_OF_POINTS=" << nPoints << std::endl;
    mean << "DIMENSION=1" << std::endl << "TYPE=Scalar" << std::endl;

    if( ArrayCurvMax && ArrayCurvMin )
      {
      // std::cout << "Got arrays" << std::endl;
      for( unsigned int i = 0; i < nPoints; i++ )
        {
        //
        float curvmax, curvmin;
        curvmax = ArrayCurvMax->GetValue(i);
        curvmin = ArrayCurvMin->GetValue(i);

        float aux = sqrt( (pow(curvmax, 2) + pow(curvmin, 2) ) / 2);
        curvedness << aux << std::endl;

        aux = (2 / M_PI) * (atan( (curvmax + curvmin) / (curvmax - curvmin) ) );
        if( aux != aux )
          {
          aux = 0;
          }
        shapeIndex << aux << std::endl;

        aux = curvmax * curvmin;
        gauss << aux << std::endl;

        aux = (curvmax + curvmin) / 2;
        mean << aux << std::endl;
        // std::cout << "Curvature: " << curv << std::endl;
        }
      }
    else
      {
      // std::cout << "Got no array?" << std::endl;

      curvedness.close();
      }
    shapeIndex.close();
    gauss.close();
    mean.close();

    std::cout << "Curvatures computed " << std::endl;

    }
  else if( GetDirectionFeaturesOn ) // bp2010
    {
    // Reading the first mesh (template)
    // std::cout << "Reading template mesh " << inputFilename << std::endl;
    vtkPolyDataReader *inputTemplate = vtkPolyDataReader::New();
    inputTemplate->SetFileName(inputFilename);
    inputTemplate->Update();
    vtkPolyData* polydata = inputTemplate->GetOutput();

    // Reading the point to create the directionality vectors
    vtkIdType numPoints = polydata->GetNumberOfPoints();
    // cout << "There are " << numPoints << " points." << endl;

    // Feature computing
    double point[3], vectorx[3], vectory[3], vectorz[3];
    vectorx[0] = 1; vectorx[1] = 0; vectorx[2] = 0;
    vectory[0] = 0; vectory[1] = 1; vectory[2] = 0;
    vectorz[0] = 0; vectorz[1] = 0; vectorz[2] = 1;

    // Computing the normals in the original mesh
    vtkPolyDataNormals *meshNormals = vtkPolyDataNormals::New();
    meshNormals->SetComputePointNormals(1);
    meshNormals->SetComputeCellNormals(0);
    meshNormals->SetSplitting(0);
    #if VTK_MAJOR_VERSION > 5
    meshNormals->SetInputData(inputTemplate->GetOutput() );
    #else
    meshNormals->SetInput(inputTemplate->GetOutput() );
    #endif
    meshNormals->Update();
    vtkPolyData * vtkMeshNormals = meshNormals->GetOutput();
    #if VTK_MAJOR_VERSION <= 5
    vtkMeshNormals->Update();
    #endif

    // Write normals out
    vtkDataArray *Array = vtkDataArray::SafeDownCast(vtkMeshNormals->GetPointData()->GetNormals() );

    // Creating feature vectors
    std::ofstream outfileVecX;
    std::ofstream outfileVecY;
    std::ofstream outfileVecZ;

    outfileVecX.open( "feature_unitDirx.txt" );
    outfileVecY.open( "feature_unitDiry.txt" );
    outfileVecZ.open( "feature_unitDirz.txt" );

    outfileVecX << "NUMBER_OF_POINTS = " << numPoints << std::endl;
    outfileVecX << "DIMENSION = " << 1 << std::endl;
    outfileVecX << "TYPE = Scalar" << std::endl;

    outfileVecY << "NUMBER_OF_POINTS = " << numPoints << std::endl;
    outfileVecY << "DIMENSION = " << 1 << std::endl;
    outfileVecY << "TYPE = Scalar" << std::endl;

    outfileVecZ << "NUMBER_OF_POINTS = " << numPoints << std::endl;
    outfileVecZ << "DIMENSION = " << 1 << std::endl;
    outfileVecZ << "TYPE = Scalar" << std::endl;
    for( vtkIdType i = 0; i < numPoints; i++ )
      {
      // Retrieve the particle normal from template and mesh
      Array->GetTuple(i, point);

      double dotprodX = point[0] * vectorx[0] + point[1] * vectorx[1] + point[2] * vectorx[2];
      double dotprodY = point[0] * vectory[0] + point[1] * vectory[1] + point[2] * vectory[2];
      double dotprodZ = point[0] * vectorz[0] + point[1] * vectorz[1] + point[2] * vectorz[2];

      outfileVecX << dotprodX << std::endl;
      outfileVecY << dotprodY << std::endl;
      outfileVecZ << dotprodZ << std::endl;
      }

    outfileVecX.close();
    outfileVecY.close();
    outfileVecZ.close();

    }
  else if( particleOn )
    {
    int numFiles = TestMeshFiles.size();

    // Reading the first mesh (template)
    std::cout << "Reading template mesh " << TestMeshFiles[0].c_str() << std::endl;
    vtkPolyDataReader *inputTemplate = vtkPolyDataReader::New();
    inputTemplate->SetFileName(TestMeshFiles[0].c_str() );
    inputTemplate->Update();

    // Computing the normals in the original mesh
    vtkDataArray *ArrayTemplate;
    {
    vtkPolyData *vtkMeshNormals;
    vtkPolyDataNormals *meshNormals = vtkPolyDataNormals::New();
    meshNormals->SetComputePointNormals(1);
    meshNormals->SetComputeCellNormals(0);
    meshNormals->SetSplitting(0);
    #if VTK_MAJOR_VERSION > 5
    meshNormals->SetInputData(inputTemplate->GetOutput() );
    #else
    meshNormals->SetInput(inputTemplate->GetOutput() );
    #endif
    meshNormals->Update();
    vtkMeshNormals = meshNormals->GetOutput();
   	#if VTK_MAJOR_VERSION <= 5
    vtkMeshNormals->Update();
    #endif
    // Write normals out
    unsigned int nPoints = vtkMeshNormals->GetNumberOfPoints();
    std::cout << nPoints << std::endl;
    ArrayTemplate = vtkDataArray::SafeDownCast(vtkMeshNormals->GetPointData()->GetNormals() );
    }
    // Creating the point locator to compute the closest point to
    vtkPointLocator *pointLocatorT = vtkPointLocator::New();
    pointLocatorT->SetDataSet(inputTemplate->GetOutput() );
    pointLocatorT->BuildLocator();

    int             counter = 0; unsigned int ptID;
    vtkDoubleArray *particlesT = vtkDoubleArray::New();
    particlesT->SetNumberOfComponents( 3 );

    std::ifstream in( TestMeshFiles[1].c_str() );
    while( in )
      {
      double pt[3], normal[3];
      in >> pt[0] >> pt[1] >> pt[2];
      // ptID=pointLocatorT->FindClosestPoint(pt[0],pt[1],pt[2]);
      ptID = pointLocatorT->FindClosestPoint(pt);
      ArrayTemplate->GetTuple(ptID, normal);
      particlesT->InsertTupleValue( counter, normal );
      counter++;
      }

    std::cout << particlesT->GetNumberOfTuples() << " particles read from " << TestMeshFiles[1] << std::endl;
    in.close();

    vtkIntArray *flags = vtkIntArray::New();
    flags->SetNumberOfComponents(1);
    flags->SetNumberOfValues(counter);
    for( int i = 0; i < counter; i++ )
      {
      flags->SetValue(i, 0);
      }
    // for the rest of the files
    for( int index = 2; index < numFiles; index = index + 2 )
      {

      std::cout << "Reading  mesh " << TestMeshFiles[index] << std::endl;
      vtkPolyDataReader *input = vtkPolyDataReader::New();
      input->SetFileName(TestMeshFiles[index].c_str() );
      input->Update();

      // **** START STEP 1 -> Computing the normals in the new mesh
      vtkPolyData * vtkMeshNormals;
      {
      vtkPolyDataNormals *meshNormals = vtkPolyDataNormals::New();
      meshNormals->SetComputePointNormals(1);
      meshNormals->SetComputeCellNormals(0);
      meshNormals->SetSplitting(0);
    	#if VTK_MAJOR_VERSION > 5
      meshNormals->SetInputData(input->GetOutput() );
      #else
      meshNormals->SetInput(input->GetOutput() );
      #endif
      meshNormals->Update();
      vtkMeshNormals = meshNormals->GetOutput();
      }
    	#if VTK_MAJOR_VERSION <= 5
      vtkMeshNormals->Update();
      #endif

      // Write normals out
      unsigned int nPoints = vtkMeshNormals->GetNumberOfPoints();
      std::cout << nPoints << std::endl;
      vtkDataArray *Array = vtkDataArray::SafeDownCast(vtkMeshNormals->GetPointData()->GetNormals() );
      // **** END STEP 1

      // *** START STEP 2 -> PARSING particle file
      vtkPointLocator *pointLocator = vtkPointLocator::New();
      pointLocator->SetDataSet(input->GetOutput() );
      pointLocator->BuildLocator();

      int             _counter = 0;
      vtkDoubleArray *particles = vtkDoubleArray::New();
      particles->SetNumberOfComponents( 3 );
      double ptT[3];

      std::ifstream _in( TestMeshFiles[index + 1].c_str() );
      while( _in )
        {
        double pt[3], normal[3];
        _in >> pt[0] >> pt[1] >> pt[2];
        // _ptID=pointLocator->FindClosestPoint(pt[0],pt[1],pt[2]);
        unsigned int _ptID = pointLocator->FindClosestPoint(pt);
        Array->GetTuple(_ptID, normal);
        particles->InsertTupleValue( _counter, normal );
        _counter++;
        }

      std::cout << particles->GetNumberOfTuples() << " particles read from " << TestMeshFiles[index + 1] << std::endl;
      in.close();
      // std::cout << "Number of particles...  " << _counter << std::endl;
      for( int i = 0; i < _counter; i++ )
        {
        double pt[3];
        // Retrieve the particle normal from template and mesh
        particlesT->GetTupleValue(i, ptT);
        particles->GetTupleValue(i, pt);

        double dotprod = pt[0] * ptT[0] + pt[1] * ptT[1] + pt[2] * ptT[2];;

        if( dotprod < 0.0 )
          {
          flags->SetValue(i, 1);
          // std::cout << dotprod << std::endl;
          }
        // else
        //	flags->InsertNextValue(0);
        }

      // *** END STEP 3
      }

    for(int index = 1; index < numFiles; index = index + 2 )
      {

      std::cout << "Processing " << TestMeshFiles[index].c_str() << std::endl;

      // string str ("Please split this phrase into tokens");
      std::string outFile = TestMeshFiles[index].c_str();
      outFile += "new.lpts";

      std::ofstream out( outFile.c_str() );
      std::ifstream inputFile( TestMeshFiles[index].c_str() );
      int _counter = 0;
      while( inputFile )
        {
        double pt[3];
        inputFile >> pt[0] >> pt[1] >> pt[2];
        // cout << pt[0] << " " << pt[1] << " " << pt[2] << endl;
        if( flags->GetValue(_counter) != 1 )
          {
          out << pt[0] << " " << pt[1] << " " << pt[2] << endl;
          }
        else
          {
          std::cout << "Out particle " << _counter << std::endl;
          }
        _counter++;
        }

      inputFile.close();
      out.close();
      _counter = 0;
      }

    // std::cout << "Se fini!" << std::endl;
    }
  else if( closestPointOn )
    {
    // Reading input mesh 1
    if( debug )
      {
      std::cout << "Reading input mesh 1..." << std::endl;
      }
    vtkPolyDataReader *meshReader1 = vtkPolyDataReader::New();
    meshReader1->SetFileName(inputFilename);
    meshReader1->Update();
    vtkPolyData *mesh1 = meshReader1->GetOutput();

    // Reading input mesh 2
    if( debug )
      {
      std::cout << "Reading input mesh 2..." << std::endl;
      }
    vtkPolyDataReader *meshReader2 = vtkPolyDataReader::New();
    meshReader2->SetFileName(closestPointFiles[1].c_str() );
    meshReader2->Update();
    vtkPolyData *mesh2 = meshReader2->GetOutput();

    // Reading attribute file
    if( debug )
      {
      std::cout << "Reading attribute file..." << std::endl;
      }
    std::ifstream      infile;
    char               Line[40];
    float              CurrentAttribute;
    std::vector<float> v_AttributeIn, v_AttributeOut;
    infile.open(closestPointFiles[0].c_str() );
    if( !infile )
      {
      std::cerr << "Unable to open file: " << argv[2] << std::endl;
      exit(1);
      }
    while( std::strncmp(Line, "NUMBER_OF_POINTS =", 18) && strncmp(Line, "NUMBER_OF_POINTS=", 17) )
      {
      infile.getline(Line, 40);
      }

    unsigned int NbVertices = atoi(strrchr(Line, '=') + 1);
    infile.getline(Line, 40);
    infile.getline(Line, 40);
    for( unsigned int i = 0; i < NbVertices; i++ )
      {
      infile >> CurrentAttribute;
      v_AttributeIn.push_back(CurrentAttribute);
      }
    infile.close();

    // Computing closest points
    if( debug )
      {
      std::cout << "Computing closest points..." << std::endl;
      }
    double    x[3];
    vtkIdType closestVertId;
    for( int PointId = 0; PointId < mesh2->GetNumberOfPoints(); PointId++ )
      {
      mesh2->GetPoint(PointId, x);
      closestVertId = mesh1->FindPoint(x);
      if( closestVertId >= 0 )
        {
        v_AttributeOut.push_back(v_AttributeIn[(int)closestVertId]);
        }
      else
        {
        std::cerr << "Error: no closest point for vertex" << PointId << std::endl;
        return EXIT_FAILURE;
        }
      }

    // Writing attribute file
    if( debug )
      {
      std::cout << "Writing attribute file..." << std::endl;
      }
    std::ofstream outfile;
    outfile.open(outputFilename);
    outfile << "NUMBER_OF_POINTS=" << v_AttributeOut.size() << std::endl;
    outfile << "DIMENSION=" << 1 << std::endl;
    outfile << "TYPE=Scalar" << std::endl;
    for( unsigned int i = 0; i < v_AttributeOut.size(); i++ )
      {
      outfile << v_AttributeOut[i] << std::endl;
      }
    outfile.close();

    }
  else if( extractVerticesOn )
    {
    // Reading input mesh
    if( debug )
      {
      std::cout << "Reading input mesh 1..." << std::endl;
      }
    vtkPolyDataReader *meshReader = vtkPolyDataReader::New();
    meshReader->SetFileName(inputFilename);
    meshReader->Update();
    vtkPolyData *mesh = meshReader->GetOutput();

    std::ofstream outfile0, outfile1, outfile2;

    outfile0.open(outputFilename);
    outfile1.open(extractVerticesFiles[0].c_str() );
    outfile2.open(extractVerticesFiles[1].c_str() );

    int NbPoints = mesh->GetNumberOfPoints();
    outfile0 << "NUMBER_OF_POINTS=" << NbPoints << std::endl << "DIMENSION=1" << std::endl << "TYPE=Scalar"
             << std::endl;
    outfile1 << "NUMBER_OF_POINTS=" << NbPoints << std::endl << "DIMENSION=1" << std::endl << "TYPE=Scalar"
             << std::endl;
    outfile2 << "NUMBER_OF_POINTS=" << NbPoints << std::endl << "DIMENSION=1" << std::endl << "TYPE=Scalar"
             << std::endl;

    double x[3];
    for( int PointId = 0; PointId < NbPoints; PointId++ )
      {
      mesh->GetPoint(PointId, x);
      outfile0 << x[0] << std::endl;
      outfile1 << x[1] << std::endl;
      outfile2 << x[2] << std::endl;
      }
    outfile0.close();
    outfile1.close();
    outfile2.close();

    } else if (lobarSurfaceAreaOn)
    {
      // Reading the input mesh
      MeshConverterType::Pointer converter = MeshConverterType::New();
      if (debug) std::cout << "Reading input mesh " << inputFilename << std::endl;
      MeshSOType::Pointer inputMeshSO = dynamic_cast<MeshSOType *>(converter->ReadMeta(inputFilename).GetPointer() );
      MeshType::Pointer   inputMesh = inputMeshSO->GetMesh();
      MeshType::PointsContainerPointer points = inputMesh->GetPoints();

      // Reading the input parcellation attribute file
      if (debug) std::cout << "Reading attribute file " << inputFilename << std::endl;
      std::ifstream Infile;
      char Line[40];
      int CurrentAttribute;
      std::vector<int> v_VerticesAttributes; // vector with attribute values
      std::vector<int> v_Attributes; //Attributes vector
      Infile.open(AttributeFileName);
      while ( strncmp (Line, "NUMBER_OF_POINTS =", 18) && strncmp (Line, "NUMBER_OF_POINTS=", 17))
	Infile.getline (Line, 40);
      unsigned int NbVertices = atoi(strrchr(Line,'=')+1);
      if (inputMesh->GetNumberOfPoints() != NbVertices)
	{
	  std::cerr<<"Input mesh and attribute file must have the same number of vertices!"<<std::endl;
	  exit(-1);
	}
      Infile.getline ( Line, 40);
      Infile.getline ( Line, 40);

      for (unsigned int i = 0; i < NbVertices; i++ )
   	{
	  Infile >> CurrentAttribute;
	  v_VerticesAttributes.push_back(CurrentAttribute);
	  bool NewAttribute = true;
	  for (unsigned int j = 0; j < v_Attributes.size(); j++)
	    {
	      if (CurrentAttribute == v_Attributes[j])
		NewAttribute = false;
	    }
	  if (NewAttribute)
	    v_Attributes.push_back(CurrentAttribute);
	}
      Infile.close();      

      // Sort attributes vector into ascending order
      sort(v_Attributes.begin(), v_Attributes.end());
     
      // Computing surface area
      typedef CellType::PointIdIterator PointIdIterator;
      CellIterator cellIterator = inputMesh->GetCells()->Begin();
      CellIterator cellEnd = inputMesh->GetCells()->End();
      PointType curPoint;
      std::vector<double> v_Area(*(std::max_element(v_Attributes.begin(),v_Attributes.end()))+1);

      while( cellIterator != cellEnd )
	{
	  CellType * cell = cellIterator.Value();
	  TriangleType * line = dynamic_cast<TriangleType *> (cell);
	  LineType::PointIdIterator pit = line->PointIdsBegin();
	  int pIndex1,pIndex2,pIndex3;
	  PointType p1, p2, p3;
	  itk::Vector<float, 3> v1, v2;

	  pIndex1= *pit;
	  p1 = points->GetElement(pIndex1);
	  ++pit;
	  pIndex2= *pit;
	  p2 = points->GetElement(pIndex2);
	  ++pit;
	  pIndex3= *pit;
	  p3 = points->GetElement(pIndex3);

	  //Computing vectors of the current triangle
	  v1 = p2-p1;
	  v2 = p3-p1;
	  
	  // Computing triangle area
	  double area = sqrt(pow(v1[1]*v2[2]-v2[1]*v1[2],2)+pow(v1[2]*v2[0]-v2[2]*v1[0],2)+pow(v1[0]*v2[1]-v2[0]*v1[1],2))/2.0;
	  if ((v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex2]) && (v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex3]))
	    v_Area[v_VerticesAttributes[pIndex1]] += area;
	  else if((v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex2]) && (v_VerticesAttributes[pIndex1] != v_VerticesAttributes[pIndex3]))
	    {
	      v_Area[v_VerticesAttributes[pIndex1]] += 2.0*area/3.0;
	      v_Area[v_VerticesAttributes[pIndex3]] += area/3.0;
	    }
	  else if((v_VerticesAttributes[pIndex1] != v_VerticesAttributes[pIndex2]) && (v_VerticesAttributes[pIndex1] == v_VerticesAttributes[pIndex3]))
	    {
	      v_Area[v_VerticesAttributes[pIndex1]] += 2.0*area/3.0;
	      v_Area[v_VerticesAttributes[pIndex2]] += area/3.0;
	    }
	  else if((v_VerticesAttributes[pIndex1] != v_VerticesAttributes[pIndex2]) && (v_VerticesAttributes[pIndex2] == v_VerticesAttributes[pIndex3]))
	    {
	      v_Area[v_VerticesAttributes[pIndex1]] += area/3.0;
	      v_Area[v_VerticesAttributes[pIndex2]] += 2.0*area/3.0;
	    }
	  else
	    {
	      v_Area[v_VerticesAttributes[pIndex1]] += area/3.0;
	      v_Area[v_VerticesAttributes[pIndex2]] += area/3.0;
	      v_Area[v_VerticesAttributes[pIndex3]] += area/3.0;
	    }
	  ++cellIterator;
	}

      if (debug) std::cout<<"Writing output file..."<<std::endl;
      std::ofstream OutputFile(outputFilename, ios::out);  
      if (!OutputFile) 
	{
	  std::cerr<<"Error: cannot open "<<outputFilename<<std::endl;
	  exit(-1);
	}
      OutputFile<<"##########################################################################"<<std::endl;
      OutputFile<<"# File format :"<<std::endl;
      OutputFile<<"# LABEL \t AREA"<<std::endl;
      OutputFile<<"# Fields :"<<std::endl;
      OutputFile<<"#\t LABEL         Label number"<<std::endl;
      OutputFile<<"#\t AREA          Surface area in square mm"<<std::endl;
      OutputFile<<"##########################################################################"<<std::endl<<std::endl; 
      for (unsigned int i = 0; i < v_Attributes.size(); i++)
	OutputFile<<v_Attributes[i]<<"\t"<<v_Area[v_Attributes[i]]<<std::endl;
      OutputFile.close();

    } else if (meanOn)
    {
      std::ifstream inFile;
      std::ofstream outFile;
      float Attribute;
      std::vector<float> v_Attribute;
      char Line[40];
      unsigned int NbVertices = 0;
      
      for (int FileNb = 0; FileNb < nbfile; FileNb++)
	{
	  // Reading file
	  if (debug) std::cout<<"Reading file: "<<meanFiles[FileNb]<<std::endl;	  
	  inFile.open(meanFiles[FileNb].c_str());
	  if (!inFile)
	    {
	      std::cerr << "Unable to open file"<<std::endl;
	      exit(1);
	    }	  
	  while ( std::strncmp (Line, "NUMBER_OF_POINTS =", 18) && strncmp (Line, "NUMBER_OF_POINTS=", 17))
	    inFile.getline (Line, 40);	  
	  if (FileNb == 0 )
	    {
	      NbVertices = atoi(strrchr(Line,'=')+1);
	      inFile.getline(Line, 40);
	      inFile.getline(Line, 40);      
	      
	      while(!inFile.eof())
		{
		  inFile >> Attribute;
		  v_Attribute.push_back(Attribute);
		}
	    }
	  else
	    {
	      unsigned int NbPoints = atoi(strrchr(Line,'=')+1);
	      if (NbPoints != NbVertices)
		{
		  std::cerr << "Attribute Files don't have the same number of vertices!"<<std::endl;
		  exit(1);
		}
	      inFile.getline(Line, 40);
	      inFile.getline(Line, 40);      
	  
	      unsigned int Point = 0;
	      while(!inFile.eof())
		{
		  inFile >> Attribute;
		  v_Attribute[Point] += Attribute;
		  Point++;
		}
	    }
	  inFile.close();
	}
      
      //Computing mean attribute
      if (debug) std::cout<<"Computing mean attribute..."<<std::endl;
      for (unsigned int i = 0; i < NbVertices; i++)
	v_Attribute[i] /= nbfile;
      
      // Writing output file
      if (debug) std::cout<<"Writing output file..."<<std::endl;
      outFile.open(outputFilename);
      if (!outFile)
	{
	  std::cerr << "Unable to open file";
	  exit(1);
	}
      outFile<<"NUMBER_OF_POINTS="<<NbVertices<<std::endl;
      outFile<<"DIMENSION="<<1<<std::endl;
      outFile<<"TYPE=Scalar"<<std::endl;
      for (unsigned int i = 0; i < NbVertices; i++)
	outFile<<v_Attribute[i]<<std::endl;
      outFile.close();
      
    } else if(cart2baryOn)
    {
      // std::cout<<inputFilename << "  " << outputFilename << "   " << vtkPointFile.c_str() << std::endl;

      // Reading input mesh 1
      vtkPolyDataReader *meshReader1 = vtkPolyDataReader::New();
      meshReader1->SetFileName(inputFilename);
      meshReader1->Update();
      vtkPolyData *mesh1 = meshReader1->GetOutput();

      // Reading input mesh 2 that contains the point list
      vtkPolyDataReader *meshReader2 = vtkPolyDataReader::New();
      meshReader2->SetFileName(vtkPointFile.c_str());
      meshReader2->Update();
      vtkPolyData *mesh2 = meshReader2->GetOutput();

      vtkPointLocator* pointLocator = vtkPointLocator::New();
      pointLocator->SetDataSet(mesh1);
      pointLocator->BuildLocator();

      vtkPoints* pointList = mesh2->GetPoints();
      mesh1->BuildLinks(); // BuildLinks needs to be called before GetPointCells

      double currLandmark[3], closestPoint[3];

      std::ofstream outfile (outputFilename);
      outfile << pointList->GetNumberOfPoints() << std::endl;

      for (int i = 0; i < pointList->GetNumberOfPoints(); i++) {

	pointList->GetPoint(i, currLandmark);
	int id = pointLocator->FindClosestPoint(currLandmark);
	mesh1->GetPoint(id, closestPoint);
	//std::cout << "Closest vertex: " << closestPoint[0] << "  "  << closestPoint[1] << "  "  << closestPoint[2] << ", ";

	vtkIdList* cellIds = vtkIdList::New();
	mesh1->GetPointCells(id, cellIds);
    
	/**
	 * projection a point onto cell surface
	 * 1) find cell
	 * 2) compute normal 
	 * 3) construct implicit plane function
	 * 4) projection onto the plane
	 * 5) determine the projected point is inside the cell
	 * 6) find the point with minimum distance
	 */
	double trueClosestPoint[3], projectedPoint[3], dist2 = vtkMath::Distance2BetweenPoints(currLandmark, closestPoint); //std::cout<<dist2<<std::endl;
	memcpy(trueClosestPoint, closestPoint  , sizeof(closestPoint  ));

	int minJ = cellIds->GetId(0);
	for (int j = 0; j < cellIds->GetNumberOfIds(); j++) {

	  vtkCell* cell = mesh1->GetCell(cellIds->GetId(j));

	  double t1[3], t2[3], t3[3], normal[3];
	  cell->GetPoints()->GetPoint(0, t1); //std::cout<< t1[0] << " " << t1[1] << " " << t1[2] << ", ";
	  cell->GetPoints()->GetPoint(1, t2); //std::cout<< t2[0] << " " << t2[1] << " " << t2[2] << ", ";
	  cell->GetPoints()->GetPoint(2, t3); //std::cout<< t3[0] << " " << t3[1] << " " << t3[2] << ", ";

	  vtkTriangle::ComputeNormal(t1, t2, t3, normal);
	  vtkPlane::GeneralizedProjectPoint(currLandmark, t1, normal, projectedPoint);

	  int inside = vtkTriangle::PointInTriangle(projectedPoint, t1, t2, t3, 0); //std::cout<< inside << ", ";
	  double ndist2 = vtkMath::Distance2BetweenPoints(currLandmark, projectedPoint); //std::cout<<ndist2<<std::endl;
	  if (inside == 1) {
	    if (dist2 > ndist2) {
	      memcpy(trueClosestPoint,projectedPoint , sizeof(projectedPoint));
	      cout << i << " vertex found closer projection onto a polygon (" << ndist2 << " < " << dist2 << ")" << std::endl;
	      dist2 = ndist2;
	      minJ = cellIds->GetId(j);
	    }

	  } else {
	    if (dist2 > ndist2) {
	    }
	  }
	}

	mesh2->GetPoints()->SetPoint(i, trueClosestPoint);

	std::cout <<"Closest point: " << trueClosestPoint[0] << "  "  << trueClosestPoint[1] << "  "  << trueClosestPoint[2] << "; Cell id: " << minJ << std::endl;
	vtkCell* closestCell = mesh1->GetCell(minJ);
	double t1[3], t2[3], t3[3], t1_2d[2], t2_2d[2], t3_2d[2], trueClosestPt_2d[2], bcoords[3];
	closestCell->GetPoints()->GetPoint(0, t1);
	closestCell->GetPoints()->GetPoint(1, t2);
	closestCell->GetPoints()->GetPoint(2, t3);

	vtkTriangle::ProjectTo2D(t1,t2,t3,t1_2d,t2_2d,t3_2d);
	//std::cout<< "2d coordinates: " << t1_2d[0] << " " << t1_2d[1] << ",  " << t2_2d[0] << " " << t2_2d[1] << "; " << t3_2d[0] << " " << t3_2d[1] << std::endl;
	vtkTriangle::ProjectTo2D(t1,t2,trueClosestPoint,t1_2d,t2_2d,trueClosestPt_2d);
	//std::cout<< "2d coordinates: " << t1_2d[0] << " " << t1_2d[1] << ",  " << t2_2d[0] << " " << t2_2d[1] << "; " << trueClosestPt_2d[0] << " " << trueClosestPt_2d[1] << std::endl;

	vtkTriangle::BarycentricCoords(trueClosestPt_2d, t1_2d, t2_2d, t3_2d, bcoords);
	std::cout<< "Barycentric coordinates: " << bcoords[0] << " " << bcoords[1] << " " << bcoords[2] << std::endl;

	vtkIdList * closestCellIds = closestCell->GetPointIds();
	outfile << minJ << " " << closestCellIds->GetId(0) << " " << closestCellIds->GetId(1) << " " <<  closestCellIds->GetId(2) << " " << bcoords[0] << " "  << bcoords[1] << " " << bcoords[2] << std::endl;
      }

      outfile.close();

      vtkPolyDataWriter *meshWriter = vtkPolyDataWriter::New();
    	#if VTK_MAJOR_VERSION > 5
      meshWriter->SetInputData(mesh2);
      #else
      meshWriter->SetInput(mesh2);
      #endif
      std::string projectedFilename = outputFilename;
      size_t found = projectedFilename.rfind(".");
      if(found == std::string::npos)
	projectedFilename = projectedFilename + ".projected.vtk";
      else
	projectedFilename = projectedFilename.substr(0,found) + ".projected.vtk";
     
      meshWriter->SetFileName(projectedFilename.c_str());
      meshWriter->Update();
      if(debug) std::cout<< "Writing new projected landmark points: " << std::endl;


    } else if (attSTDOn)
   {
	//std::cout << inputFilename << std::endl;
	//for (int i = 0 ; i < nbfile ; i ++)
	//{
		//std::cout << VectorStdFiles[i] << std::endl;	
	//}



	std::ifstream input;
	char line[70];
	char *aux;
	int NPoints, Dimension;
	input.open(inputFilename, ios::in);
    	input.getline(line,70,'\n');
    	aux=strtok(line, " = ");
    	aux=strtok(NULL, " = ");
    	NPoints=atoi(aux);
    	input.getline(line,70,'\n');
    	aux=strtok(line, " = ");
    	aux=strtok(NULL, " = ");
    	Dimension=atoi(aux);
    	input.getline(line,70,'\n');

	if (Dimension == 3)
	{
		std::vector < std::vector <itk::Point<double,3> > > vectorData (nbfile+1);
		std::vector <itk::Point<double,3> > vectorAvg (NPoints);
		itk::Point<double,3> pt;
		pt[0] = 0; pt[1] = 0; pt[2] = 0;
		std::vector <itk::Point<double,3> > vectorStd (NPoints,pt);
		
	
		std::cout << inputFilename << std::endl;
		vectorData[0].resize(NPoints);
		for (int i = 0 ; i < NPoints ; i ++)
		{ 
			input >> pt[0] >> pt[1] >> pt[2];
			vectorData[0][i]=pt;
			vectorAvg[i]=pt;
			//std::cout << vectorData[0][i][0] << " " << vectorData[0][i][1] << " " << vectorData[0][i][2] << std::endl;
		}
		input.close();
	
		for (int i = 0 ; i < nbfile ; i ++)
		{
	
			vectorData[i+1].resize(NPoints);
			input.open(VectorStdFiles[i].c_str(), ios::in);
			std::cout << VectorStdFiles[i] << std::endl;
			input.getline(line,70,'\n');
			input.getline(line,70,'\n');
			input.getline(line,70,'\n');
		
			for (int j = 0 ; j < NPoints ; j ++)
			{ 
				input >> pt[0] >> pt[1] >> pt[2];
				vectorData[i+1][j]=pt;
				
				for (int k=0 ; k <3 ; k ++)
				{
					double value = vectorAvg[j][k] + pt[k];
					vectorAvg[j][k] = value;
					//std::cout << value << std::endl;	 
				}
				
				//std::cout << vectorData[i+1][j][0] << " " << vectorData[i+1][j][1] << " " << vectorData[i+1][j][2] << std::endl;
			}
			input.close();	
		}
	
		std::ofstream output;
		output.open ( "vectorSTD_average.txt" ) ;
		output << "NUMBER_OF_POINTS = " << NPoints << std::endl;
		output << "DIMENSION = 3" << std::endl;
		output << "TYPE = Vector" << std::endl;
		for (int j = 0 ; j < NPoints ; j ++)
		{ 
			for (int k=0 ; k <3 ; k ++)
			{
				double value = vectorAvg[j][k] / (nbfile+1);
				vectorAvg[j][k] = value;
			}
			output << vectorAvg[j][0] << " " << vectorAvg[j][1] << " " << vectorAvg[j][2] << std::endl;
		}
		output.close();
	
		output.open ( "vectorSTD_std.txt" ) ;
		output << "NUMBER_OF_POINTS = " << NPoints << std::endl;
		output << "DIMENSION = 3" << std::endl;
		output << "TYPE = Vector" << std::endl;
		for (int j = 0 ; j < NPoints ; j ++)
		{ 
			for (int i = 0 ; i < nbfile+1 ; i ++)
			{	
				for (int k=0 ; k <3 ; k ++)
				{
					double value = (vectorData[i][j][k] - vectorAvg[j][k]) * (vectorData[i][j][k] - vectorAvg[j][k]);
					double value2 = vectorStd[j][k] + value;
					vectorStd[j][k] = value2;	 
				}
				
			}	
			
			for (int k=0 ; k <3 ; k ++)
			{
				double value = vectorStd[j][k] / (nbfile+1);
				double value2 = sqrt (value);
				vectorStd[j][k] = value2;
			}
	
			//std::cout << vectorStd[j][0] << " " << vectorStd[j][1] << " " << vectorStd[j][2] << std::endl;
			output << vectorStd[j][0] << " " << vectorStd[j][1] << " " << vectorStd[j][2] << std::endl;
		}
		output.close();
	}

	if (Dimension == 1)
	{
		std::vector < std::vector <double > > vectorData (nbfile+1);
		std::vector <double > vectorAvg (NPoints);
		std::vector <double > vectorStd (NPoints,0);
		double point;
		
	
		std::cout << inputFilename << std::endl;
		vectorData[0].resize(NPoints);
		for (int i = 0 ; i < NPoints ; i ++)
		{ 
			input >> point;
			vectorData[0][i]=point;
			vectorAvg[i]=point;
			std::cout << vectorData[0][i] << std::endl;
		}
		input.close();
	
		for (int i = 0 ; i < nbfile ; i ++)
		{
	
			vectorData[i+1].resize(NPoints);
			input.open(VectorStdFiles[i].c_str(), ios::in);
			std::cout << VectorStdFiles[i] << std::endl;
			input.getline(line,70,'\n');
			input.getline(line,70,'\n');
			input.getline(line,70,'\n');
		
			for (int j = 0 ; j < NPoints ; j ++)
			{ 
				input >> point;
				vectorData[i+1][j]=point;
				
				double value = vectorAvg[j] + point;
				vectorAvg[j] = value;
				
				std::cout << vectorData[i+1][j] << std::endl;
			}
			input.close();	
		}
	
		std::ofstream output;
		output.open ( "scalarSTD_average.txt" ) ;
		output << "NUMBER_OF_POINTS = " << NPoints << std::endl;
		output << "DIMENSION = 1" << std::endl;
		output << "TYPE = Scalar" << std::endl;
		for (int j = 0 ; j < NPoints ; j ++)
		{ 
			double value = vectorAvg[j] / (nbfile+1);
			vectorAvg[j] = value;
			output << vectorAvg[j] << std::endl;
		}
		output.close();
	
		output.open ( "scalarSTD_std.txt" ) ;
		output << "NUMBER_OF_POINTS = " << NPoints << std::endl;
		output << "DIMENSION = 1" << std::endl;
		output << "TYPE = Scalar" << std::endl;
		for (int j = 0 ; j < NPoints ; j ++)
		{ 
			for (int i = 0 ; i < nbfile+1 ; i ++)
			{	
				double value = (vectorData[i][j] - vectorAvg[j]) * (vectorData[i][j] - vectorAvg[j]);
				double value2 = vectorStd[j] + value;
				vectorStd[j] = value2;	 
			}	
			
			double value = vectorStd[j] / (nbfile+1);
			double value2 = sqrt (value);
			vectorStd[j] = value2;
	
			output << vectorStd[j] << std::endl;
		}
		output.close();
	}
	
   } else if (lptsOn)
   {
	//std::cout << inputFilename << " " << argv[2] << " " << outputFilename << std::endl;
        int Pointwritten;
	Pointwritten=0;
	int nbPoints=0;

	if(nbPoints ==0)
	{
		vtkPolyDataReader *meshin = vtkPolyDataReader::New();
		meshin->SetFileName(inputFilename);
		meshin->Update();
		nbPoints=meshin->GetOutput()->GetNumberOfPoints();
	}

	std::ofstream lptsfile(outputFilename, std::ios::out | std::ios::trunc);

	vtkPolyDataReader *vtkreader = vtkPolyDataReader::New();
	vtkreader->SetFileName(inputFilename);
	vtkreader->Update();
	double x[3];
	vtkPoints * PointVTK = vtkPoints::New();
	PointVTK=vtkreader->GetOutput()->GetPoints();
	for (int PointId = 0; PointId < (vtkreader->GetOutput()->GetNumberOfPoints()); PointId++)
	{
		PointVTK->GetPoint(PointId,x);
		lptsfile<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
	}
	lptsfile.close();
	
   } else if(relaxPolygonsOn) //bp2012
  {
	vtkPolyDataReader *meshin = vtkPolyDataReader::New();
	meshin->SetFileName(inputFilename);
	meshin->Update();

	vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
	#if VTK_MAJOR_VERSION > 5
	meshinC->SetInputData(meshin->GetOutput());
  #else
	meshinC->SetInput(meshin->GetOutput());
  #endif
	meshinC->Update();

	vtkWindowedSincPolyDataFilter *smoother = vtkWindowedSincPolyDataFilter::New();
	smoother->SetInputConnection(meshinC->GetOutputPort());
	smoother->SetNumberOfIterations(relaxPolygonsNb);
	smoother->BoundarySmoothingOff();
	smoother->FeatureEdgeSmoothingOff();
	smoother->SetFeatureAngle(120.0);
	smoother->SetPassBand(0.001);
	smoother->NonManifoldSmoothingOn();
	smoother->NormalizeCoordinatesOn();
	smoother->Update();


 	//Writing the new mesh, with not degenerated triangles
	vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
	#if VTK_MAJOR_VERSION > 5
	SurfaceWriter->SetInputData(smoother->GetOutput());
  #else
	SurfaceWriter->SetInput(smoother->GetOutput());
  #endif
	SurfaceWriter->SetFileName(outputFilename);
	SurfaceWriter->Update();

	if (debug) std::cout << "Writing new mesh " << outputFilename << std::endl;
  } else if(decimateMeshOn) //bp2012
  {
	vtkPolyDataReader *meshin = vtkPolyDataReader::New();
	meshin->SetFileName(inputFilename);
	meshin->Update();

	vtkCleanPolyData *meshinC = vtkCleanPolyData::New();
	#if VTK_MAJOR_VERSION > 5
	meshinC->SetInputData(meshin->GetOutput());
  #else
	meshinC->SetInput(meshin->GetOutput());
  #endif
	meshinC->Update();

  	vtkSmartPointer<vtkPolyData> input = vtkSmartPointer<vtkPolyData>::New();
  	input->ShallowCopy(meshinC->GetOutput());

	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	#if VTK_MAJOR_VERSION <= 5
    decimate->SetInputConnection(input->GetProducerPort());
	#else
    decimate->SetInputData(input);
	#endif
  	decimate->SetTargetReduction(decimateMeshTR); //10% reduction (if there was 100 triangles, now there will be 90)
  	decimate->Update();
 
  	vtkSmartPointer<vtkPolyData> decimated = vtkSmartPointer<vtkPolyData>::New();
  	decimated->ShallowCopy(decimate->GetOutput());

 	//Writing the new mesh, with not degenerated triangles
	vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
  #if VTK_MAJOR_VERSION > 5
	SurfaceWriter->SetInputData(decimate->GetOutput());  
  #else
	SurfaceWriter->SetInput(decimate->GetOutput());  
  #endif
	SurfaceWriter->SetFileName(outputFilename);
	SurfaceWriter->Update();

	if (debug) std::cout << "Writing new mesh " << outputFilename << std::endl;
  } else if (listPointDataOn) //ms2013
  {
    // Read PolyData info
    vtkPolyDataReader *polyIn = vtkPolyDataReader::New();
    if (debug) std::cout << "Reading vtk mesh " << inputFilename << std::endl;
    polyIn->SetFileName(inputFilename); 
    polyIn->Update();
    vtkPolyData* polydata = polyIn->GetOutput();
    unsigned int    nPoints = polydata->GetNumberOfPoints();;

    unsigned int numberOfArrays = polydata->GetPointData()->GetNumberOfArrays();
    std::cout << "Number of Scalars " << numberOfArrays << ", number of Points " << nPoints << std::endl;
     
    for (unsigned int j = 0 ; j < numberOfArrays ; j ++)
      { 
	std::cout << "  name " << j << ": " << polydata->GetPointData()->GetArrayName(j) <<std::endl;

	polydata->GetPointData()->SetActiveScalars(polydata->GetPointData()->GetArrayName(j));
	vtkDoubleArray *scalarData = (vtkDoubleArray *) polydata->GetPointData()->GetArray(j);

	if (scalarData) 
	  {
	    double min = scalarData->GetComponent(0,0);
	    double max = min;
	    for (unsigned int i = 0 ; i < nPoints; i++)
	      {
		double val = scalarData->GetComponent(i,0);
		if (min > val)  min = val;
		if (max < val)  max = val;
	      }
	
	    std::cout << "     min/max :" << min << " / " << max <<std::endl;
	  }
      }
  } else if (PointDataOpOn)
  {
    
    // Read PolyData info
    vtkPolyDataReader *polyIn = vtkPolyDataReader::New();
    if (debug) std::cout << "Reading vtk mesh " << inputFilename << std::endl;
    polyIn->SetFileName(inputFilename); 
    polyIn->Update();
    vtkPolyData* polydata = polyIn->GetOutput();
    unsigned int    nPoints = polydata->GetNumberOfPoints();;

    unsigned int numberOfArrays = polydata->GetPointData()->GetNumberOfArrays();
    if (debug) std::cout << "Number of Scalars " << numberOfArrays << ", number of Points " << nPoints << std::endl;

    const char * name = pointDataOpCmd[0].c_str();
    polydata->GetPointData()->SetActiveScalars(name);
    vtkDoubleArray *scalarData = (vtkDoubleArray *) polydata->GetPointData()->GetArray(name);
    if (scalarData) 
      {
	double opValue = atof(pointDataOpCmd[2].c_str());
	if (debug)  std::cout << "operation " << pointDataOpCmd[1] << " : " << opValue << std::endl;
     
	if (! pointDataOpCmd[1].compare("threshBelow"))
	{
	  for (unsigned int i = 0 ; i < nPoints; i++)
	    {
	      if (scalarData->GetComponent(i,0) < opValue)  scalarData->SetComponent(i,0, 0);
	    }
	
	} else if (! pointDataOpCmd[1].compare("sub")) 
	{
	  for (unsigned int i = 0 ; i < nPoints; i++)
	    {
	      double curValue = scalarData->GetComponent(i,0);
	      scalarData->SetComponent(i,0, curValue - opValue);
	    }
	} else 
	{
	  std::cout << "No (known) point operation detected -> exiting" << std::endl;
	}
      }
    if (debug)  std::cout << "saving modified vtk in " << outputFilename  << std::endl;
     
    vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
    #if VTK_MAJOR_VERSION > 5
    SurfaceWriter->SetInputData(polydata);
    #else
    SurfaceWriter->SetInput(polydata);
    #endif
    SurfaceWriter->SetFileName(outputFilename);
    SurfaceWriter->Update();

  } else  if( translateOn )
    {

     std::cout << "Displacing mesh by " << TranslateParams[0] << " - "  << TranslateParams[1] << " - " << TranslateParams[2] << std::endl;

     vtkPolyDataReader *meshin = vtkPolyDataReader::New();
     meshin->SetFileName(inputFilename);
     meshin->Update();
     double x[3];
     vtkPoints * PointsVTK = vtkPoints::New();
     PointsVTK = meshin->GetOutput()->GetPoints();

     for( int PointId = 0; PointId < (meshin->GetOutput()->GetNumberOfPoints() ); PointId++ )
        {
        PointsVTK->GetPoint(PointId, x);
        double vert[3];
        for( unsigned int dim = 0; dim < 3; dim++ )
          {
          vert[dim] = x[dim] + TranslateParams[dim];
          }
        PointsVTK->SetPoint(PointId, vert);
        }

     vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
     #if VTK_MAJOR_VERSION > 5
     SurfaceWriter->SetInputData(meshin->GetOutput() );
     #else
     SurfaceWriter->SetInput(meshin->GetOutput() );
     #endif
     SurfaceWriter->SetFileName(outputFilename);
     SurfaceWriter->Update();
     
    }
    else if( lookupPointOn ) // bp2009
  {
  vtkPolyDataReader *polyRead = vtkPolyDataReader::New();
  polyRead->SetFileName(inputFilename);
  polyRead->Update();
  vtkPolyData* polyOut = polyRead->GetOutput();

  if( debug )
    {
    std::cout << "Input Mesh read ...  " << inputFilename << std::endl;
    }

  // *** START PARSING the csv file
  std::ifstream csvfile(files[0]);
  if( !csvfile.good() )
  {
    std::cerr << "Could not open CSV file: " << files[0] << std::endl ;
    std::cerr << "Abort" << std::endl ;
    return 1 ;
  }
  vtkSmartPointer<vtkDelimitedTextReader> CSVreader = vtkSmartPointer<vtkDelimitedTextReader>::New();
  CSVreader->SetFieldDelimiterCharacters(",");
  CSVreader->SetFileName(files[0]);
  CSVreader->SetHaveHeaders(true);
  CSVreader->Update();
  vtkSmartPointer<vtkTable> table = CSVreader->GetOutput();

  vtkDoubleArray* outputScalars = (vtkDoubleArray *) polyOut->GetPointData()->GetArray(files[1]);

  if (outputScalars)
  {
      if (debug)
      {
          std::cout << "Got array!" << std::endl;
      }
  }
  else
  {
    std::cerr << "Could not find field \'"<< files[1] << "\' in VTK file \'" << inputFilename << "\'"<< std::endl ;
    std::cerr << "Abort" << std::endl ;
    return 1 ;
  }

  std::cout << "Processing point data ... " << files[1] << std::endl;

  for (unsigned int idx=0; idx < table->GetNumberOfRows() ; idx++)
  {
    for (unsigned int idxs = 0 ; idxs < polyOut->GetNumberOfPoints() ; idxs++)
    {
        if (outputScalars->GetVariantValue(idxs) == table->GetValue(idx,0))
        {
            outputScalars->SetVariantValue(idxs,table->GetValue(idx,1));
        }
    }
  }
  if( debug )
    {
    std::cout << "Scalar map added to the mesh" << std::endl;
    }


  // Writing the new mesh
  vtkPolyDataWriter *SurfaceWriter = vtkPolyDataWriter::New();
  #if VTK_MAJOR_VERSION > 5
  SurfaceWriter->SetInputData(polyOut);
  #else
  SurfaceWriter->SetInput(polyOut);
  #endif
  SurfaceWriter->SetFileName(outputFilename);
  SurfaceWriter->Update();
  if( debug )
    {
    std::cout << "Writing new mesh " << outputFilename << std::endl;
    }

  }
   else {

    std::cout << "No operation to do -> exiting" << std::endl;
    exit(-1);
    }
  return EXIT_SUCCESS;
}
