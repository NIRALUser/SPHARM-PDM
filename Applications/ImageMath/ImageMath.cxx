/* 
 * compute image math and combinations
 *
 * author:  Martin Styner 
 * 
 * changes:
 *
 */

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#pragma warning ( disable : 4503 )
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <set>


using namespace std;

#include <string.h>
#include <sys/types.h>
#include <stdlib.h>    // for exit, system
#include <math.h>
#include <errno.h>

#include <itkImage.h>
#include <itkImageFileReader.h> 
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>
#include <itkImageToImageFilter.h>

#include <itkThresholdImageFilter.h> 
#include <itkBinaryThresholdImageFilter.h> 
#include <itkBinaryBallStructuringElement.h> 
#include <itkBinaryDilateImageFilter.h> 
#include <itkBinaryErodeImageFilter.h> 

#include <itkImageRegionIterator.h>
#include <itkImageRegionConstIterator.h>
#include <itkNeighborhoodIterator.h>
#include <itkCastImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkAddImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <itkDivideImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>

#include <itkCurvatureFlowImageFilter.h> 
#include <itkDiscreteGaussianImageFilter.h> 
#include <itkGrayscaleDilateImageFilter.h> 
#include <itkGrayscaleErodeImageFilter.h> 
#include <itkGrayscaleFillholeImageFilter.h>
#include <itkMeanImageFilter.h> 
#include <itkGradientAnisotropicDiffusionImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkRelabelComponentImageFilter.h>
#include <itkExtractImageFilter.h>
#include <itkMaximumImageFilter.h>
#include <itkMinimumImageFilter.h>
#include <itkMinimumMaximumImageCalculator.h>

#include <itkDiffusionTensor3D.h>
#include <itkResampleImageFilter.h>
#include <itkAffineTransform.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkImageRegionIterator.h>
#include <itkMetaDataObject.h>

#include "argio.h"
#include "ImageMath.h" 

#define DEFAULT_SAMP 2
// number of samples by default

using namespace std;
using namespace itk;

typedef float PixelType;
typedef short ShortPixelType;
typedef unsigned char BinaryPixelType;
enum { ImageDimension = 3 };
typedef Image<PixelType,ImageDimension>       ImageType;
typedef Image<ShortPixelType,ImageDimension>  ShortImageType;
typedef Image<BinaryPixelType,ImageDimension> BinaryImageType;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageRegionIterator< ImageType >      IteratorType;
typedef ImageRegionIterator< ShortImageType> ShortIteratorType;
typedef ImageRegionConstIterator<ImageType>   ConstIteratorType;
typedef NeighborhoodIterator<ImageType>       NeighborhoodIteratorType;
typedef ImageType::Pointer                    ImagePointer;
typedef ImageToImageFilter<ImageType,BinaryImageType> Conv_float2binType;
typedef ImageToImageFilter<BinaryImageType,ImageType> Conv_bin2floatType;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;
typedef ImageFileWriter< BinaryImageType >    BinaryVolumeWriterType;
typedef ImageFileWriter< ShortImageType >     ShortVolumeWriterType;

typedef CastImageFilter< ImageType, BinaryImageType > castBinaryFilterType;
typedef CastImageFilter< ImageType,  ShortImageType > castShortFilterType; 

typedef ThresholdImageFilter< ImageType > maskThreshFilterType;
typedef BinaryThresholdImageFilter< ImageType , ImageType > threshFilterType;
typedef BinaryThresholdImageFilter< ShortImageType , ImageType > ShortthreshFilterType;
typedef MaskImageFilter< ImageType, ImageType, ImageType >  maskFilterType;
typedef AddImageFilter< ImageType, ImageType,  ImageType > addFilterType;
typedef SubtractImageFilter< ImageType, ImageType,  ImageType > subFilterType;
typedef MultiplyImageFilter< ImageType, ImageType,  ImageType > mulFilterType;
typedef DivideImageFilter< ImageType, ImageType,  ImageType > divFilterType;

typedef BinaryBallStructuringElement<PixelType,ImageDimension> StructuringElementType;

typedef BinaryDilateImageFilter<ImageType, ImageType, StructuringElementType> dilateFilterType;
typedef BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType> erodeFilterType;
typedef HistogramMatchingImageFilter< ImageType , ImageType > matchHistogramFilterType;

typedef MeanImageFilter<ImageType, ImageType> meanFilterType;
typedef CurvatureFlowImageFilter<ImageType, ImageType> curvFilterType;
typedef DiscreteGaussianImageFilter<ImageType, ImageType> gaussFilterType;
typedef GrayscaleErodeImageFilter<ImageType, ImageType, StructuringElementType> erodegrayFilterType;
typedef GrayscaleDilateImageFilter<ImageType, ImageType, StructuringElementType> dilategrayFilterType;
typedef GrayscaleFillholeImageFilter<ImageType, ImageType> fillholegrayFilterType;
typedef GradientAnisotropicDiffusionImageFilter< ImageType, ImageType > anisoDiffFilterType;

//typedef ConnectedComponentImageFilter<ImageType,ImageType> ConnectiveFilterType;
typedef ConnectedComponentImageFilter<ImageType,ShortImageType> ConnectiveFilterType;
typedef RelabelComponentImageFilter<ShortImageType,ShortImageType> RelabelFilterType;

typedef ExtractImageFilter<ImageType, ImageType> CropFilterType;

typedef MaximumImageFilter<ImageType, ImageType> MaximumImageFilterType;
typedef MinimumImageFilter<ImageType, ImageType> MinimumImageFilterType;

typedef itk::MinimumMaximumImageCalculator<ImageType> MaxFilterType;

static int debug = 0;


//What pixeltype is the image 
void GetImageType( char* fileName ,
                   itk::ImageIOBase::IOPixelType &pixelType ,
                   itk::ImageIOBase::IOComponentType &componentType )
{
  typedef itk::Image< unsigned char , 3 > ImageType ;
  itk::ImageFileReader< ImageType >::Pointer imageReader =
    itk::ImageFileReader< ImageType >::New();
  imageReader->SetFileName( fileName ) ;
  imageReader->UpdateOutputInformation() ;
  pixelType = imageReader->GetImageIO()->GetPixelType() ;
  componentType = imageReader->GetImageIO()->GetComponentType() ;
}


int main(const int argc, const char **argv)
{
  if (argc <=1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
    cout << "ImageMath 1.1 version (Dec 2007)" << endl;
    cout << " computes base statistics of an image" << endl;
    cout << "usage: ImageMath infile <options>" << endl;
    cout << endl;
    cout << "infile                 input dataset" << endl;;
    cout << "-outbase outbase       base-outputfilename, if omitted then the same as base of input" << endl; 
    cout << "-outfile outfile       outfilename, will only be applied to main output (if multiple output are given)" << endl;
    cout << "-nocomp                don't automatically compress the output image" << endl;
    cout << "-combine infile2       combine the inputfile by interpreting them as labelfiles. " << endl 
       << "   -relabel              labels in infile2 will be relabeled to succeed labels in infile (no label overlap)" << endl;
    cout << "       labels in infile2 overwrite only the background label in infile1" << endl; 
    cout << "-extractLabel label    extract the mentioned label from the file" << endl;
    cout << "-threshold min,max     threshold: everything I < min , I > max is set to 0, otherwise 1" << endl;
    cout << "-threshMask min,max    mask threshold: everything I < min , I > max is set to 0, otherwise is left as is" << endl;
    cout << "-mask infile2          use infile2 as a binary mask (intensities > 0 ) and combine with inputfile  " << endl ;
    cout << "-constOper opID,val    apply the following operation to the image: I op val, op = +/0, -/1, */2, //3" << endl;
    cout << "-add infile2           apply the following operation to the image: I1 + I2" << endl;
    cout << "-sub infile2           apply the following operation to the image: I1 - I2" << endl;
    cout << "-mul infile2           apply the following operation to the image: I1 * I2" << endl;
    cout << "-div infile2           apply the following operation to the image: I1 / I2" << endl;
    cout << "-pwr arg               compute each voxels to power arg" << endl;
    cout << "-normalizeEMS csf,wm,gmfiles  (infile should be grayscale template) normalizes the EMS prob maps" << endl;
    cout << "-editPixdims px,py,pz   simply change the pixdims (without reslicing) of the image" << endl;
    cout << "-dilate radius,val      apply isotropic dilation with ball element of given radius and value" << endl;
    cout << "-erode radius,val       apply isotropic erosion with ball element of given radius and value" << endl;
    cout << "-matchHistogram infile2  match the image histogram to the one in infile2" << endl;
    cout << "-matchHistoPara bins,points,thresh  optional parameters for matchHistogram (bins= number of histogram bins" << endl;
    cout << "                        points = number of control points, thresh = boolean for mean intensity threshol [0/1])" << endl;
    cout << "-smooth -gauss -curveEvol -grayOpen -grayClose -grayDilate -grayErode -grayFillHole -meanFilter -anisoDiff [-size val] [-iter num]" << endl;
    cout << "                       smoothing of image using any of the mentioned smoothing filters" << endl;
    cout << "                       size is stuctural element size, variance, or timestep depending on the filter choice" << endl;
    cout << "                       iter is number of iterations" << endl;
    cout << "-v                     verbose mode " << endl;
    cout << "-type byte|short|float Type of processing image (computations are always done with float images), default is short" << endl; 
    cout << "-extension ext         Extension of the output (determines output format)" << endl; 
    cout << "-conComp Lbl           For a binary image, rank all the connected components according to their sizes and create a binary image with the 'Lbl' biggest ones" << endl;
    cout << "                       if Lbl=0 outputs the labeled image with all the components labeled according to their size" << endl;
    cout << "-changeOrig px,py,pz   Change the orgine of the image" << endl;
    cout << "-createIm X,Y,Z,px,py,pz Create an empty image with the specified parameters: image dimensions X,Y,Z, image resolution px,py,pz" << endl;
    cout << " -crop px,py,pz,w,h,d   cropimage: origin px,py,pz (startindex is 0) dimensions width(w),height(h), depth(d)" << endl;
    cout << "-changeSp spx,spy,spz  Change the spacing of the image" << endl;
    cout << "-max infile2           Compute the maximum between the corresponding pixels of the two input images" << endl;
    cout << "-min infile2           Compute the minimum between the corresponding pixels of the two input images" << endl;
    cout << "-avg infile2 infile3...       Compute the average image" << endl;
    cout << "-majorityVoting infile2 infile3...       Compute an accurate parcellation map considering a majority voting process" << endl;
    cout << "-center                Center image" << endl;
    cout << "-flip [x,][y,][z]      Flip image" << endl;
    cout << endl << endl;
    exit(0);
  }

  const int BGVAL = 0;
  const int FGVAL = 1;

  char *inputFileName = strdup(argv[1]);
  char *outputFileName = ipGetStringArgument(argv, "-outfile", NULL);  
  char *outbase    = ipGetStringArgument(argv, "-outbase", NULL);  
  char * base_string;
  if (!outbase) {
    if (!outputFileName) {
      base_string = strdup(ipGetBaseName(inputFileName));
    } else {
      base_string = strdup(ipGetBaseName(outputFileName));
    }
  } else {
    base_string = outbase;
  }
  string outFileName ("dummy");
  char *typeChat       = ipGetStringArgument(argv, "-type", NULL);
  bool writeFloat= false;
  bool writeByte = false;
  if (typeChat && !strcmp(typeChat,"byte")) writeByte = true;
  if (typeChat && !strcmp(typeChat,"float"))writeFloat = true;

  char * formatChar = ipGetStringArgument(argv, "-extension", ".gipl");
  string format;
  if (! strchr(formatChar, '.')) {
    format = string(".") + string(formatChar);
  } else {
    format = string(formatChar);
  }

  debug      = ipExistsArgument(argv, "-v");

  bool nocompOn = false;
  nocompOn = ipExistsArgument(argv, "-nocomp");
  bool center = false;
  center = ipExistsArgument(argv, "-center");

  char *combineFile    = ipGetStringArgument(argv, "-combine", NULL);  
  bool relabelOn = ipExistsArgument(argv, "-relabel");

  bool extractLabelOn   = ipExistsArgument(argv, "-extractLabel"); 
  int extractLabel   = ipGetIntArgument(argv, "-extractLabel", 0); 

  bool maskOn   = ipExistsArgument(argv, "-mask"); 
  char *maskFile    = ipGetStringArgument(argv, "-mask", NULL);  

  bool addOn   = ipExistsArgument(argv, "-add"); 
  char *addFile    = ipGetStringArgument(argv, "-add", NULL); 
  bool subOn   = ipExistsArgument(argv, "-sub"); 
  char *subFile    = ipGetStringArgument(argv, "-sub", NULL);  
  bool mulOn   = ipExistsArgument(argv, "-mul"); 
  char *mulFile    = ipGetStringArgument(argv, "-mul", NULL); 
  bool divOn   = ipExistsArgument(argv, "-div"); 
  char *divFile    = ipGetStringArgument(argv, "-div", NULL); 

  bool thresholdOn    = ipExistsArgument(argv, "-threshold"); 
  char * tmp_str      = ipGetStringArgument(argv, "-threshold", NULL);
  PixelType tmin = 0;
  PixelType tmax = 0;
  float textend[2];
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "threshold needs 2 comma separated entries: min,max" << endl;
      exit(1);
    } else {
      tmin = (PixelType) textend[0];
      tmax = (PixelType) textend[1];
    }
  }

  bool threshMaskOn    = ipExistsArgument(argv, "-threshMask"); 
  tmp_str      = ipGetStringArgument(argv, "-threshMask", NULL);
  PixelType tmaskmin = 0;
  PixelType tmaskmax = 0;
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "mask threshold needs 2 comma separated entries: min,max" << endl;
      exit(1);
    } else {
      tmaskmin = (PixelType) textend[0];
      tmaskmax = (PixelType) textend[1];
    }
  }

  int erodeRadius, dilateRadius;
  erodeRadius = dilateRadius = 1;
  PixelType erodeVal, dilateVal;
  bool dilateOn   = ipExistsArgument(argv, "-dilate");  
  tmp_str      = ipGetStringArgument(argv, "-dilate", NULL); 
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "dilate needs 2 comma separated entries: radius,value" << endl;
      exit(1);
    } else {
      dilateRadius = (int) textend[0];
      dilateVal = (PixelType) textend[1];
    }
  }
  bool erodeOn   = ipExistsArgument(argv, "-erode");  
  tmp_str      = ipGetStringArgument(argv, "-erode", NULL); 
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "erode needs 2 comma separated entries: radius,value" << endl;
      exit(1);
    } else {
      erodeRadius = (int) textend[0];
      erodeVal = (PixelType) textend[1];
    }
  }

  bool editPixdimsOn    = ipExistsArgument(argv, "-editPixdims"); 
  tmp_str      = ipGetStringArgument(argv, "-editPixdims", NULL);
  float pixdims[3];
  if (tmp_str) { 
    int num = ipExtractFloatTokens(pixdims, tmp_str, 3);
    if (3 != num) {
      cerr << "editPixdims needs 3 comma separated entries: px,py,pz " << endl;
      exit(1); 
    } 
  }  

  bool constOperOn    = ipExistsArgument(argv, "-constOper"); 
  tmp_str      = ipGetStringArgument(argv, "-constOper", NULL);
  int operID = 0;
  PixelType operVal = 0;
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 2);
    if (2 != num) {
      cerr << "oper needs 2 comma separated entries: opID,val" << endl;
      exit(1);
    } else {
      operID = (int) textend[0];
      operVal = (PixelType) textend[1];
    }
  }

  bool connectiveCompOn  =  ipExistsArgument(argv, "-conComp");
  tmp_str =     ipGetStringArgument(argv, "-conComp", NULL);
  int Lbl = 1;
  if(tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 1);
    if(1 != num){
      cerr << "conCompt only needs 1 value" << endl;
      exit(1);
    } else {
      Lbl = (int) textend[0];
    }
  }

  bool normalizeEMSOn    = ipExistsArgument(argv, "-normalizeEMS"); 
  tmp_str      = ipGetStringArgument(argv, "-normalizeEMS", NULL);
  char * probFiles[3];
  string  csfFile, wmFile, gmFile;
  if (tmp_str) {
    int num = ipExtractStringTokens(probFiles, tmp_str, 3);
    if (3 != num) {
      cerr << "normalizeEMS needs 3 comma separated entries: csf,wm,gm" << endl;
      exit(1);
    } else {
      csfFile = string(probFiles[0]);
      wmFile = string(probFiles[1]);
      gmFile = string(probFiles[2]);
    }
  }
  
  bool changeOrigOn    = ipExistsArgument(argv, "-changeOrig"); 
  tmp_str      = ipGetStringArgument(argv, "-changeOrig", NULL); 
  float origCoor[3];
  if (tmp_str) { 
    int num = ipExtractFloatTokens(textend, tmp_str, 3);
    if (3 != num) {
      cerr << "changeOrig needs 3 comma separated entries: px,py,pz " << endl;
      exit(1); 
    } else {
      origCoor[0] = (int)textend[0];
      origCoor[1] = (int)textend[1];
      origCoor[2] = (int)textend[2];
    } 
  }
  bool changeSpOn    = ipExistsArgument(argv, "-changeSp");
  tmp_str      = ipGetStringArgument(argv, "-changeSp", NULL);
  float spacingval[3];
  if (tmp_str) {
    int num = ipExtractFloatTokens(textend, tmp_str, 3);
    if (3 != num) {
      cerr << "changeSp needs 3 comma separated entries: spx,spy,spz " << endl;
      exit(1);
    } else {
      spacingval[0] = static_cast<float>(textend[0]);
      spacingval[1] = static_cast<float>(textend[1]);
      spacingval[2] = static_cast<float>(textend[2]);
    }
  }

  bool matchHistogramOn = ipExistsArgument(argv, "-matchHistogram");
  char *matchHistogramFile    = ipGetStringArgument(argv, "-matchHistogram", NULL); 
  tmp_str      = ipGetStringArgument(argv, "-matchHistoPara", NULL);
  int matchHistoPara[3];
  int matchHistoNumBins, matchHistoNumPoints, matchHistoThresh;
  if (tmp_str) {
    int num = ipExtractIntTokens(matchHistoPara, tmp_str, 3);
    if (3 != num) {
      cerr << "matchHistoPara needx 3 comma separated entries: bins,points,threshBool" << endl;
      exit(1);
    } else {
      matchHistoNumBins = matchHistoPara[0];
      matchHistoNumPoints = matchHistoPara[1];
      matchHistoThresh = matchHistoPara[2];
    }
  } else {
      matchHistoNumBins = 1024;
      matchHistoNumPoints = 50;
      matchHistoThresh = 1;
  }

  bool imageCreationOn = ipExistsArgument(argv, "-createIm");
  tmp_str = ipGetStringArgument(argv, "-createIm", NULL); 
  float Dims[6]; 
  if(tmp_str) {
    int num = ipExtractFloatTokens(Dims, tmp_str, 6);
    if(6 != num){
      cerr << "createIm needs 6 parameters separated by commas" << endl;
    } else { 
    
      std::cout << "Val: " << Dims[0] << " | " << Dims[1] << " | " << Dims[2] << " | " << Dims[3] << " | " << Dims[4] << " | " << Dims[5] << std::endl;
    }
  }

  const int numCropParam = 6;
  int cropParam[numCropParam];
  bool cropOn = ipExistsArgument(argv,"-crop");
  if (cropOn) {  
    char *tmp_str    = ipGetStringArgument(argv, "-crop", NULL);
    int numDim       = ipExtractIntTokens(cropParam, tmp_str, numCropParam);
    if (numDim != numCropParam) {              
      cerr << argv[0] << ": crop needs "<< numCropParam << " parameters.\n";
      exit(-1);
    }
    free(tmp_str);
  } 

  bool smoothOn =  ipExistsArgument(argv,"-smooth"); 
  bool gaussianOn =  ipExistsArgument(argv,"-gauss"); 
  bool curveEvolOn = ipExistsArgument(argv,"-curveEvol"); 
  bool grayOpenOn = ipExistsArgument(argv,"-grayOpen"); 
  bool grayCloseOn = ipExistsArgument(argv,"-grayClose"); 
  bool grayDilateOn = ipExistsArgument(argv,"-grayDilate"); 
  bool grayErodeOn = ipExistsArgument(argv,"-grayErode"); 
  bool grayFillHoleOn = ipExistsArgument(argv, "-grayFillHole");
  bool meanFilterOn = ipExistsArgument(argv,"-meanFilter"); 
  bool anisoDiffOn = ipExistsArgument(argv,"-anisoDiff"); 
  if (smoothOn && !gaussianOn && !curveEvolOn 
		&& !grayOpenOn && !grayCloseOn && !grayDilateOn && !grayErodeOn && !grayFillHoleOn
		&& !meanFilterOn && !anisoDiffOn) { 
    curveEvolOn = true; 
  };
  double smoothSize = ipGetDoubleArgument(argv,"-size",-1);
  unsigned int numIter = ipGetIntArgument(argv,"-iter",1);

  bool MaxOn   = ipExistsArgument(argv, "-max"); 
  char *MaxFile    = ipGetStringArgument(argv, "-max", NULL);
  
  bool MinOn   = ipExistsArgument(argv, "-min"); 
  char *MinFile    = ipGetStringArgument(argv, "-min", NULL);

  bool pwrOn   = ipExistsArgument(argv, "-pwr");
  tmp_str = ipGetStringArgument(argv, "-pwr", NULL);
  float pwrval = 1;
  if (tmp_str) {
    int num = ipExtractFloatTokens(&pwrval, tmp_str, 1);
    if (1 != num) {
      cerr << "pwr option requires one entry: the value of the power" << endl;
      exit(1);
    }
  }

  const int MaxNumFiles = 1000;
  int NbFiles = 0;
  char *files[MaxNumFiles];
  vector<string> InputFiles;

  bool AvgOn = ipExistsArgument(argv, "-avg");
  if (AvgOn)
    {
      NbFiles = ipGetStringMultipArgument(argv, "-avg", files, MaxNumFiles);
      for(int i = 0 ; i < NbFiles ; i++)
	InputFiles.push_back(files[i]);
    }
  
  bool MajorityVotingOn = ipExistsArgument(argv, "-majorityVoting");
  if (MajorityVotingOn)
    {
      NbFiles = ipGetStringMultipArgument(argv, "-majorityVoting", files, MaxNumFiles);
      for(int i = 0 ; i < NbFiles ; i++)
	InputFiles.push_back(files[i]);
    }
  
  bool flip   = ipExistsArgument(argv, "-flip"); 
  tmp_str      = ipGetStringArgument(argv, "-flip", NULL);
  itk::Matrix< double , 3 , 3 > flipMatrix ;
  flipMatrix.SetIdentity() ;
  if( tmp_str )
  {
    bool x = false ;
    bool y = false ;
    bool z = false ;
    char * flipAxes[3];
    int num = ipExtractStringTokens(flipAxes, tmp_str, 3);
    for( int i = 0 ; i < num ; i++ )
    {
      if( !strcmp( flipAxes[ i ] , "x" ) )
      {
        if( x == true )//'x' appears multiple times
        {
          std::cout<<"x appears multiple times"<<std::endl;
          return EXIT_FAILURE ;
        }
        flipMatrix[ 0 ][ 0 ] = -1 ;
        x = true ;
      }
      else if( !strcmp( flipAxes[ i ] , "y" ) )
      {
        if( y == true )//'y' appears multiple times
        {
          std::cout<<"y appears multiple times"<<std::endl;
          return EXIT_FAILURE ;
        }
        flipMatrix[ 1 ][ 1 ] = -1 ;
        y = true ;
      }
      else if( !strcmp( flipAxes[ i ] , "z" ) )
      {
        if( z == true )//'z' appears multiple times
        {
          std::cout<<"z appears multiple times"<<std::endl;
          return EXIT_FAILURE ;
        }
        flipMatrix[ 2 ][ 2 ] = -1 ;
        z = true ;
      }
      else
      {
        std::cout<<"Error: can only flip along x, y, z or any combination of them"<<std::endl;
        return EXIT_FAILURE ;        
      }
    }
  }

  // **********************************************
  // **********************************************
  // **********************************************
  // Cmd line parsing done
  // **********************************************
  // **********************************************
  // **********************************************
  ImagePointer inputImage ;
  typedef itk::ImageBase< 3 > ImageBaseType ;
  ImageBaseType::Pointer inputBaseImage ;
  // Check the type of image that is loaded
  itk::ImageIOBase::IOPixelType pixelType ;
  itk::ImageIOBase::IOComponentType componentType ;
  GetImageType( inputFileName , pixelType , componentType ) ;
  bool diffusionImage = false ;
  //If tensor image
  typedef itk::Image< itk::DiffusionTensor3D< PixelType > , 3 > DiffusionImageType ;
  if (debug) cout << "Loading file " << inputFileName << endl;
  if( pixelType == itk::ImageIOBase::SYMMETRICSECONDRANKTENSOR
   || pixelType == itk::ImageIOBase::DIFFUSIONTENSOR3D
   || pixelType == itk::ImageIOBase::VECTOR
    )
  {
    if( !changeSpOn && !changeOrigOn && !editPixdimsOn && !cropOn && !center  && !maskOn )
    {
      std::cerr << "The only operations supported on Diffusion Tensor Images are: -editPixdims, -changeSp, -changeOrig, -crop , -mask and -center"<< std::endl ;
      return EXIT_FAILURE;
    }
    typedef itk::ImageFileReader< DiffusionImageType > DiffusionReaderType ;
    DiffusionReaderType::Pointer diffusionReader = DiffusionReaderType::New() ;
    diffusionReader->SetFileName( inputFileName ) ;
    try
    {
      diffusionReader->Update();
    }
    catch (ExceptionObject & err)
    {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }
    inputBaseImage = diffusionReader->GetOutput() ;
    diffusionImage = true ;
  }
  else
  {
 // load image
//  if (debug) cout << "Loading file " << inputFileName << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(inputFileName) ;
    try
    {
      imageReader->Update();
    }
    catch (ExceptionObject & err)
    {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }    
    inputImage = imageReader->GetOutput();
    inputBaseImage = inputImage ;
  }

 
  
  // do something to InputImage
  if (combineFile) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_comb");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << combineFile << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(combineFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }    
    inputImage2 = imageReader->GetOutput();
    
    IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
    IteratorType iterImage2 (inputImage2, inputImage->GetBufferedRegion());
    
    if (relabelOn) {
      if (debug) cout << "Relabeling image2 " << endl;
      // relabel image2 to contain consecutive labels starting at the max label in image1
      // get max in image1
      PixelType max = iterImage1.Get();
      while ( !iterImage1.IsAtEnd() )  {
      PixelType value =  iterImage1.Get();
      if (max < value) max = value;
      ++iterImage1;
      }
      iterImage1.GoToBegin();
      // analyze labels in image2
      typedef set<PixelType> labelSetType;
      labelSetType labelSet; // sets are automatically sorted and contain only unique entries
      labelSetType::iterator pos;
      while ( !iterImage2.IsAtEnd() )  {
        PixelType value =  iterImage2.Get();
        if (value) {
          labelSet.insert(value);
        }
        ++iterImage2;
      }
      iterImage2.GoToBegin();

      // relabel labels in image2
      if (debug) {
	std::cout << " Relabeling map : " << std::endl;
        pos = labelSet.begin();
        int intpos = 1;
        while (pos != labelSet.end()){
	  std::cout << *pos << " ---> " << max + intpos << std::endl;
          intpos++;
          pos++;
        }
        
      }
      while ( !iterImage2.IsAtEnd() )  {
      PixelType value =  iterImage2.Get();
      if (value) {
        pos = labelSet.begin();
        int intpos = 1;
        while (pos != labelSet.end() && *pos != value){
          intpos++;
          pos++;
        }
        iterImage2.Set(max + intpos);
        
      }
      ++iterImage2;
      }
      iterImage2.GoToBegin();
    }

    // combine them now
    if (debug) cout << "Combining images " << endl;
    while ( !iterImage1.IsAtEnd() )  {
      PixelType value1 =  iterImage1.Get();
      PixelType value2 =  iterImage2.Get();
      
      if (!value1 && value2) {
       iterImage1.Set(value2);
      }
      ++iterImage1;
      ++iterImage2;
    }
  } else if (extractLabelOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_label");
    if (debug) cout << "extracting object " << extractLabel << endl; 
    
    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetLowerThreshold(extractLabel);
    threshFilter->SetUpperThreshold(extractLabel);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->SetInsideValue (FGVAL);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }        
    inputImage = threshFilter->GetOutput();
    
  }  else if (thresholdOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_thresh");
    if (debug) cout << "thresholding image  " << tmin << "," << tmax << endl; 
    
    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetLowerThreshold(tmin);
    threshFilter->SetUpperThreshold(tmax);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->SetInsideValue (FGVAL);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }        
    inputImage = threshFilter->GetOutput();
  }  else if (threshMaskOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_threshMask");
    if (debug) cout << "threshold/mask image  " << tmaskmin << "," << tmaskmax << endl; 
    
    maskThreshFilterType::Pointer threshFilter = maskThreshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->ThresholdOutside(tmaskmin,tmaskmax);
    try {
      threshFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }        
    inputImage = threshFilter->GetOutput();
  } else if (maskOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_mask");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << maskFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(maskFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;
    } 
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "masking images  " << endl;

    if( diffusionImage )
    {
      DiffusionImageType::Pointer inputDiffusionImage =
          dynamic_cast< DiffusionImageType* >( inputBaseImage.GetPointer() ) ;
      ImageType::SizeType maskSize ;
      maskSize = inputImage2->GetLargestPossibleRegion().GetSize() ; 
      ImageType::SizeType size ;
      size = inputDiffusionImage->GetLargestPossibleRegion().GetSize() ;
      //Check that diffusion input volume and mask volume have the same size
      for( int i = 0 ; i < 3 ; i++ )
      {
        if( size[ i ] != maskSize[ i ] )
        {
          std::cout<<"Mask and input diffusion volume must have the same size"<<std::endl ;
          return EXIT_FAILURE ;
        }
      }
      //Create output volume and fill it with null tensors
      DiffusionImageType::Pointer outputVolume = DiffusionImageType::New() ;
      outputVolume->CopyInformation( inputDiffusionImage ) ;
      outputVolume->SetRegions( size ) ;
      outputVolume->Allocate() ;
      itk::DiffusionTensor3D< PixelType > pixel ;
      pixel.Fill( (PixelType)0.0 ) ;
      outputVolume->FillBuffer( pixel ) ;
      //Create iterators
      typedef itk::ImageRegionIterator< DiffusionImageType > DiffusionIterator ;
      typedef itk::ImageRegionIterator< ImageType > MaskIterator ;
      DiffusionIterator it( inputDiffusionImage , inputDiffusionImage->GetLargestPossibleRegion() ) ;
      MaskIterator maskIt(  inputImage2 ,  inputImage2->GetLargestPossibleRegion() ) ;
      DiffusionIterator out( outputVolume , outputVolume->GetLargestPossibleRegion() ) ;
      //Copy tensors where mask is not null
      for( it.GoToBegin() , maskIt.GoToBegin() , out.GoToBegin() ;
           !it.IsAtEnd() ; ++it , ++maskIt , ++out )
      {
        if( maskIt.Get() )
        {
          out.Set( it.Get() ) ;
        }
      }
      outputVolume->SetMetaDataDictionary( inputBaseImage->GetMetaDataDictionary() ) ;
      inputBaseImage = outputVolume ;
    }
    else
    {
      maskFilterType::Pointer maskFilter = maskFilterType::New() ;
      maskFilter->SetInput1( inputImage ) ;
      maskFilter->SetInput2( inputImage2 ) ;
      try
      {
        maskFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl ;
        cerr << err << endl ;
        return EXIT_FAILURE ;
      }
      inputImage = maskFilter->GetOutput();
    }    
  } else if (addOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_add");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << addFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(addFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }    
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "adding images  " << endl;

    addFilterType::Pointer addFilter = addFilterType::New();
    addFilter->SetInput1(inputImage);
    addFilter->SetInput2(inputImage2);
    try {
      addFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage = addFilter->GetOutput();
    
  } else if (subOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_sub");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << subFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(subFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "subing images  " << endl;

    subFilterType::Pointer subFilter = subFilterType::New();
    subFilter->SetInput1(inputImage);
    subFilter->SetInput2(inputImage2);
    try {
      subFilter->Update();    
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage = subFilter->GetOutput();

  } else if (mulOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_mul");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << mulFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(mulFile) ;
    try {
      imageReader->Update();   
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "muling images  " << endl;

    mulFilterType::Pointer mulFilter = mulFilterType::New();
    mulFilter->SetInput1(inputImage);
    mulFilter->SetInput2(inputImage2);
    try {
       mulFilter->Update();  
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage = mulFilter->GetOutput();

  } else if (divOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_div");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << divFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(divFile) ;
    try {
      imageReader->Update();   
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "diving images  " << endl;

    divFilterType::Pointer divFilter = divFilterType::New();
    divFilter->SetInput1(inputImage);
    divFilter->SetInput2(inputImage2);
    try {
      divFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage = divFilter->GetOutput();

  } else if (matchHistogramOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_matchHisto");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << matchHistogramFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(matchHistogramFile) ;
    try {
      imageReader->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "matching images  " << endl;

    matchHistogramFilterType::Pointer matchHistogramFilter = matchHistogramFilterType::New();
    matchHistogramFilter->SetSourceImage(inputImage);
    matchHistogramFilter->SetReferenceImage(inputImage2);
    matchHistogramFilter->SetThresholdAtMeanIntensity(matchHistoThresh);
    matchHistogramFilter->SetNumberOfHistogramLevels(matchHistoNumBins); // number of bins
    matchHistogramFilter->SetNumberOfMatchPoints(matchHistoNumPoints); // number of equally distributed match points
    try {
      matchHistogramFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    inputImage = matchHistogramFilter->GetOutput();

  } else if (constOperOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_oper");
    if (debug) cout << "Operation  " << operID << "," << operVal << endl; 
    
    IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
    while ( !iterImage1.IsAtEnd() )  {
      PixelType value1 =  iterImage1.Get();
      PixelType value2;
      if (operID == 0) { value2 = value1 + operVal;
      } else if (operID == 1) { value2 = value1 - operVal;
      } else if (operID == 2) { value2 = value1 * operVal;
      } else if (operID == 3) { value2 = value1 / operVal;
      }
      
      iterImage1.Set(value2);
      ++iterImage1;
    }
  } else if (connectiveCompOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_conComp"); 

    if (debug) cout << "Get all the  " <<Lbl  << " bigger elements " << endl; 

    ConnectiveFilterType::Pointer Connective = ConnectiveFilterType::New();
    RelabelFilterType::Pointer relabelFilter = RelabelFilterType::New();
    ShortthreshFilterType::Pointer ThresFilter = ShortthreshFilterType::New();
    //Get the connectivity map of the image
    Connective->SetInput(inputImage);
    try {
      Connective->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    } 
    
    ShortImageType::Pointer binImage = Connective->GetOutput();
    ShortIteratorType it (binImage, binImage->GetLargestPossibleRegion());
    IteratorType itout (inputImage, inputImage->GetLargestPossibleRegion());      
    it.GoToBegin();
    itout.GoToBegin();
    while(! it.IsAtEnd())
      {
	ShortPixelType pix = it.Get();
	PixelType outpix = static_cast<PixelType>(pix);	   
	if(Lbl == 0)
	  itout.Set(outpix);
	else {	      
	  if(pix <= Lbl && pix > 0)
	    {
	      itout.Set(outpix);
	      if(pix != 0)
		std::cout << "Pix: " << pix << " Lbl: " << Lbl << " Outputval: " << outpix << std::endl;
	    }
	  else
	    itout.Set(0);
	}
	
	++it;
	++itout;
      }       
  } else if (dilateOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_dilate");

    if (debug) cout << "dilate ball radius  " <<dilateRadius  << ", val " << dilateVal << endl; 
    StructuringElementType structuringElement;
    structuringElement.SetRadius( dilateRadius );  // 3x3x3 structuring element
    structuringElement.CreateStructuringElement( );
      
    dilateFilterType::Pointer dilateFilter = dilateFilterType::New(); 
    dilateFilter->SetInput(inputImage);
    dilateFilter->SetDilateValue (dilateVal);
    dilateFilter->SetKernel( structuringElement );
    try {
      dilateFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }   
    
    inputImage = dilateFilter->GetOutput();
  } else if (changeOrigOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_newOrig");
//    inputImage->SetOrigin(origCoor);
    inputBaseImage->SetOrigin(origCoor);
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if (changeSpOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_newSpacing");
    std::cout << "SPACING: " << spacingval[0] << " | " << spacingval[1] << " | " << spacingval[2]  << std::endl;
//    inputImage->SetSpacing(spacingval);
    inputBaseImage->SetSpacing(spacingval);
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if (center) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_centered");
    ImageType::SizeType size ;
    size = inputBaseImage->GetLargestPossibleRegion().GetSize() ;
    ImageType::PointType origin ;
    origin = inputBaseImage->GetOrigin() ;
    itk::Index< 3 > index ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      index[ i ] = size[ i ] - 1 ;
    }
    itk::Point< double , 3 > corner ;
    inputBaseImage->TransformIndexToPhysicalPoint( index , corner ) ;
    itk::Point< double, 3 > newOrigin ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      newOrigin[ i ] = ( origin[ i ] - corner[ i ] ) / 2.0 ;
    }
    inputBaseImage->SetOrigin( newOrigin ) ;
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if (erodeOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_erode");

    if (debug) cout << "dilate ball radius  " << erodeRadius << ", val " << erodeVal << endl; 
    StructuringElementType structuringElement;
    structuringElement.SetRadius( erodeRadius );  // 3x3x3 structuring element
    structuringElement.CreateStructuringElement( );
      
    erodeFilterType::Pointer erodeFilter = erodeFilterType::New(); 
    erodeFilter->SetInput(inputImage);
    erodeFilter->SetErodeValue (erodeVal);
    erodeFilter->SetKernel( structuringElement );
    try {
      erodeFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }    
    
    inputImage = erodeFilter->GetOutput();
    
  } else if(cropOn)
    {
    ImageRegionType extractRegion;
    extractRegion.SetIndex(0,cropParam[0]);
    extractRegion.SetIndex(1,cropParam[1]);
    extractRegion.SetIndex(2,cropParam[2]);
    extractRegion.SetSize(0,cropParam[3]);
    extractRegion.SetSize(1,cropParam[4]);
    extractRegion.SetSize(2,cropParam[5]);
    
    int dim[3];
    ImageRegionType imageRegion = inputBaseImage->GetLargestPossibleRegion();
    dim[0] = imageRegion.GetSize(0);
    dim[1] = imageRegion.GetSize(1);
    dim[2] = imageRegion.GetSize(2);
    if (debug) cout << "size of the original image " << dim[0] << "," << dim[1] << "," << dim[2]<< endl; 
    if (debug) cout << "cropping (x,y,z,w,h,d) " << extractRegion.GetIndex(0) << "," 
		    << extractRegion.GetIndex(1) << "," << extractRegion.GetIndex(2) << "," 
		    << extractRegion.GetSize(0) << "," << extractRegion.GetSize(1) << "," 
		    << extractRegion.GetSize(2) <<  endl;
    if( diffusionImage )
    {
      typedef ExtractImageFilter< DiffusionImageType , DiffusionImageType > DiffusionCropFilterType ;
      DiffusionCropFilterType::Pointer cropFilter = DiffusionCropFilterType::New() ;
      cropFilter->SetInput( dynamic_cast< DiffusionImageType* >( inputBaseImage.GetPointer() ) ) ;
      cropFilter->SetExtractionRegion( extractRegion ) ;
      try
      {
        cropFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;	
      }
      cropFilter->GetOutput()->SetMetaDataDictionary( inputBaseImage->GetMetaDataDictionary() ) ;
      inputBaseImage = cropFilter->GetOutput() ;

    }
    else
    {
      CropFilterType::Pointer cropFilter = CropFilterType::New() ;
      cropFilter->SetInput( inputImage ) ;
      cropFilter->SetExtractionRegion( extractRegion ) ;
      try 
      {
        cropFilter->Update() ;
      }
      catch (ExceptionObject & err)
      {
        cerr << "ExceptionObject caught!" << endl;
        cerr << err << endl;
        return EXIT_FAILURE;	
      }
      inputImage = cropFilter->GetOutput();
    }
  } else if (smoothOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_smooth");
    
    if (curveEvolOn) {
      if (smoothSize == -1) smoothSize = 0.125;
      cout << "smoothing  (curveEvolution): size " << smoothSize << ", iterations " << numIter << endl;
      curvFilterType::Pointer smoothFilter = curvFilterType::New();  
      
      smoothFilter->SetInput(inputImage);
      smoothFilter->SetNumberOfIterations(numIter);
      smoothFilter->SetTimeStep(smoothSize);
      try {
	smoothFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      inputImage = smoothFilter->GetOutput();
    } else if (gaussianOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (Gaussian): size " << smoothSize << ", iterations " << numIter << endl;
      gaussFilterType::Pointer smoothFilter = gaussFilterType::New();  
      for (unsigned int i=0; i < numIter; i++) {
	smoothFilter->SetInput(inputImage);
	smoothFilter->SetVariance(smoothSize);
	try {
	  smoothFilter->Update();
	}
	catch (ExceptionObject & err) {
	  cerr << "ExceptionObject caught!" << endl;
	  cerr << err << endl;
	  return EXIT_FAILURE;	
	}    
	inputImage = smoothFilter->GetOutput();
      }
	} else if (grayDilateOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayDilate): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      dilategrayFilterType::Pointer dilateFilter = dilategrayFilterType::New();  
      
      for (unsigned int i=0; i < numIter; i++) {
		dilateFilter->SetInput(inputImage);
		structuringElement.SetRadius( (int) smoothSize ); 
		structuringElement.CreateStructuringElement( );
		dilateFilter->SetKernel( structuringElement );
		try {
		  dilateFilter->Update();
		}
		catch (ExceptionObject & err) {
		  cerr << "ExceptionObject caught!" << endl;
		  cerr << err << endl;
		  return EXIT_FAILURE;	
		}    
		inputImage = dilateFilter->GetOutput();
      }
	} else if (grayErodeOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayErode): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      erodegrayFilterType::Pointer erodeFilter = erodegrayFilterType::New();  
      
      for (unsigned int i=0; i < numIter; i++) {
		erodeFilter->SetInput(inputImage);
		structuringElement.SetRadius( (int) smoothSize ); 
		structuringElement.CreateStructuringElement( );
		erodeFilter->SetKernel( structuringElement );
		try {
		  erodeFilter->Update();
		}
		catch (ExceptionObject & err) {
		  cerr << "ExceptionObject caught!" << endl;
		  cerr << err << endl;
		  return EXIT_FAILURE;	
		}    
		inputImage = erodeFilter->GetOutput();
      }
    } else if (grayCloseOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayClose): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      dilategrayFilterType::Pointer dilateFilter = dilategrayFilterType::New();  
      erodegrayFilterType::Pointer erodeFilter = erodegrayFilterType::New();  
      
      for (unsigned int i=0; i < numIter; i++) {
	dilateFilter->SetInput(inputImage);
	erodeFilter->SetInput(dilateFilter->GetOutput());
	structuringElement.SetRadius( (int) smoothSize ); 
	structuringElement.CreateStructuringElement( );
	dilateFilter->SetKernel( structuringElement );
	erodeFilter->SetKernel( structuringElement );
	try {
	  erodeFilter->Update();
	}
	catch (ExceptionObject & err) {
	  cerr << "ExceptionObject caught!" << endl;
	  cerr << err << endl;
	  return EXIT_FAILURE;	
	}    
	inputImage = erodeFilter->GetOutput();
      }
    } else if (grayOpenOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayOpen): size " << smoothSize << ", iterations " << numIter << endl;
      StructuringElementType structuringElement;
      dilategrayFilterType::Pointer dilateFilter = dilategrayFilterType::New();  
      erodegrayFilterType::Pointer erodeFilter = erodegrayFilterType::New();  
      
      for (unsigned int i=0; i < numIter; i++) {
	erodeFilter->SetInput(inputImage);
	dilateFilter->SetInput(erodeFilter->GetOutput());
	structuringElement.SetRadius( (int) smoothSize);  
	dilateFilter->SetKernel( structuringElement );
	erodeFilter->SetKernel( structuringElement );
	try {
	  dilateFilter->Update();
	}
	catch (ExceptionObject & err) {
	  cerr << "ExceptionObject caught!" << endl;
	  cerr << err << endl;
	  return EXIT_FAILURE;	
	}    
	inputImage = dilateFilter->GetOutput();
      }
	} else if (grayFillHoleOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (grayFillHole): size " << smoothSize << ", iterations " << numIter << endl;

	  fillholegrayFilterType::Pointer fillholeFilter = fillholegrayFilterType::New();
      for (unsigned int i=0; i < numIter; i++) {
		fillholeFilter->SetInput(inputImage);
		try {
		  fillholeFilter->Update();
		}
		catch (ExceptionObject & err) {
		  cerr << "ExceptionObject caught!" << endl;
		  cerr << err << endl;
		  return EXIT_FAILURE;	
		}    
		inputImage = fillholeFilter->GetOutput();
      }
    } else if (meanFilterOn) {
      if (smoothSize == -1) smoothSize = 1.0;
      cout << "smoothing  (meanFilter): size " << smoothSize << ", iterations " << numIter << endl;
      meanFilterType::Pointer smoothFilter = meanFilterType::New();  
      
      for (unsigned int i=0; i < numIter; i++) {
	smoothFilter->SetInput(inputImage);
	ImageType::SizeType size;
	size[0] = size[1] = size[2] = (int) smoothSize;
	smoothFilter->SetRadius( size );
	try {
	  smoothFilter->Update();
	}
	catch (ExceptionObject & err) {
	  cerr << "ExceptionObject caught!" << endl;
	  cerr << err << endl;
	  return EXIT_FAILURE;	
	}    
	inputImage = smoothFilter->GetOutput();
      }
    } else if (anisoDiffOn) {
      if (smoothSize == -1) smoothSize = 0.05;
      cout << "smoothing  (gradient anisotropic Diffusion ): size  " << smoothSize << ", iterations " << numIter << endl;
      anisoDiffFilterType::Pointer smoothFilter = anisoDiffFilterType::New();  
      
      smoothFilter->SetInput(inputImage);  
      smoothFilter->SetNumberOfIterations(numIter);  
      smoothFilter->SetConductanceParameter(0.75);
      smoothFilter->SetTimeStep(smoothSize);
      try {
	smoothFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      inputImage = smoothFilter->GetOutput();
    }
  } else if (normalizeEMSOn) {
    
    ImagePointer csfImage, wmImage, gmImage, restImage;
    {
      if (debug) cout << "Loading image " << csfFile  << endl;
      VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
      imageReader->SetFileName(csfFile.c_str()) ;
      try {
	imageReader->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      csfImage = imageReader->GetOutput();
    }
    {
      if (debug) cout << "Loading image " << wmFile  << endl;
      VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
      imageReader->SetFileName(wmFile.c_str()) ;
      try {
	imageReader->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      wmImage = imageReader->GetOutput();
    }
    {
      if (debug) cout << "Loading image " << gmFile  << endl;
      VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
      imageReader->SetFileName(gmFile.c_str()) ;
      try {
	imageReader->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      gmImage = imageReader->GetOutput();
    }
    if (debug) cout << "normalizing images  " << endl;

    addFilterType::Pointer addFilter = addFilterType::New();
    addFilter->SetInput1(csfImage);
    addFilter->SetInput2(wmImage);
    try {
      addFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }    
    restImage = addFilter->GetOutput(); 
    
    addFilterType::Pointer add2Filter = addFilterType::New();
    add2Filter->SetInput1(restImage);
    add2Filter->SetInput2(gmImage);
    try {
      add2Filter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }    
    inputImage = add2Filter->GetOutput();

    IteratorType normIter (inputImage, inputImage->GetBufferedRegion());
    IteratorType csfIter (csfImage, inputImage->GetBufferedRegion());
    IteratorType wmIter (wmImage, inputImage->GetBufferedRegion());
    IteratorType gmIter (gmImage, inputImage->GetBufferedRegion());
    IteratorType restIter (restImage, inputImage->GetBufferedRegion());

    while ( !normIter.IsAtEnd() )  {
      PixelType normVal =  normIter.Get();

      if (normVal < 255) {
         restIter.Set(255 - normVal);
      } else {
         double factor = 255/normVal;
         restIter.Set(1);
         gmIter.Set( gmIter.Get() * factor);
         wmIter.Set( wmIter.Get() * factor);
         csfIter.Set( csfIter.Get() * factor);
      }

      ++wmIter;
      ++gmIter;
      ++csfIter;
      ++restIter;
      ++normIter;
    }

    {
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_normWM");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(wmImage);
      try {
	castFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str()); 
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    {
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_normGM");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(gmImage);
      try {
	castFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str()); 
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    {
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_normRest");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(restImage);
      try {
	castFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str()); 
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    {
      outFileName.erase();
      outFileName.append(base_string);
      outFileName.append("_normCSF");
      outFileName.append(format);
      castShortFilterType::Pointer castFilter = castShortFilterType::New();
      castFilter->SetInput(csfImage);
      try {
	castFilter->Update();
      }
      catch (ExceptionObject & err) {
	cerr << "ExceptionObject caught!" << endl;
	cerr << err << endl;
	return EXIT_FAILURE;	
      }    
      ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
      if(!nocompOn) writer->UseCompressionOn();
      writer->SetFileName(outFileName.c_str()); 
      writer->SetInput(castFilter->GetOutput());
      writer->Write();
    }
    exit(0);
  } else if (editPixdimsOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_pixDims");

//    inputImage->SetSpacing(pixdims);
    inputBaseImage->SetSpacing(pixdims);
    inputImage = dynamic_cast< ImageType* >( inputBaseImage.GetPointer() ) ;
  } else if(imageCreationOn) {

    ImageType::Pointer image = ImageType::New();

    ImageType::SizeType ImDim;
    ImDim[0] = static_cast<int>(Dims[0]);
    ImDim[1] = static_cast<int>(Dims[1]);
    ImDim[2] = static_cast<int>(Dims[2]);
    //    inputImage->SetSize(ImDim);
    ImageType::SpacingType spacing;
    spacing[0] = Dims[3];
    spacing[1] = Dims[4];
    spacing[2] = Dims[5];
    image->SetSpacing(spacing);
    image->SetRegions(ImDim);
    image->Allocate();

    inputImage = image;

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_createIm");
  } else if (MaxOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_max");
    
    if (debug) cout << "Loading file2 " << MaxFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(MaxFile) ;
    if (debug) cout << "computing maximum  " << endl;

    MaximumImageFilterType::Pointer MaximumFilter = MaximumImageFilterType::New();
    MaximumFilter->SetInput1(inputImage);
    MaximumFilter->SetInput2(imageReader->GetOutput());
    try {
      MaximumFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }
    inputImage = MaximumFilter->GetOutput();
    
  } else if (MinOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_min");
    
    if (debug) cout << "Loading file2 " << MinFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(MinFile) ;
    if (debug) cout << "computing maximum  " << endl;

    MinimumImageFilterType::Pointer MinimumFilter = MinimumImageFilterType::New();
    MinimumFilter->SetInput1(inputImage);
    MinimumFilter->SetInput2(imageReader->GetOutput());
    try {
      MinimumFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }
    inputImage = MinimumFilter->GetOutput();

  } else if(pwrOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_pwr");

    if(debug) std::cout << "Applying power: " << pwrval << std::endl;

    //Creating the output image
    ImageType::Pointer imagepwr = ImageType::New();
    ImageType::SizeType ImDimpwr;
    ImDimpwr[0] = (inputImage->GetLargestPossibleRegion()).GetSize()[0];
    ImDimpwr[1] = (inputImage->GetLargestPossibleRegion()).GetSize()[1];
    ImDimpwr[2] = (inputImage->GetLargestPossibleRegion()).GetSize()[2];
    ImageType::SpacingType spacingpwr;
    spacingpwr[0] = inputImage->GetSpacing()[0];
    spacingpwr[1] = inputImage->GetSpacing()[0];
    spacingpwr[2] = inputImage->GetSpacing()[0];
    imagepwr->SetSpacing(spacingpwr);
    imagepwr->SetRegions(ImDimpwr);
    imagepwr->Allocate();

    IteratorType iterImage1 (inputImage, inputImage->GetLargestPossibleRegion());
    IteratorType iterImage2 (imagepwr,imagepwr->GetLargestPossibleRegion());
    iterImage2.GoToBegin();
    iterImage1.GoToBegin();
    float outval;
    
    while ( !iterImage1.IsAtEnd() )  {
      PixelType value =  iterImage1.Get();
      outval = 0;
      if (value == 0) {
	outval = 0;
      } else {
	outval = pow(static_cast<float>(value),static_cast<float>(pwrval));
      } 
      iterImage2.Set(outval);
      ++iterImage2;
      ++iterImage1;
    }
    inputImage = imagepwr;

  } else if(AvgOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_avg");

    if(debug) cout << "Computing average image " << endl;

    VolumeReaderType::Pointer ImageReader = VolumeReaderType::New();
    addFilterType::Pointer addFilter = addFilterType::New();
    for (int FileNumber = 0; FileNumber < NbFiles; FileNumber++)
      {
	// Reading image
	ImageReader->SetFileName(InputFiles[FileNumber].c_str());
	
	// Adding image
	addFilter->SetInput1(inputImage);
	addFilter->SetInput2(ImageReader->GetOutput());
	try
	  {
	    addFilter->Update();
	  }
	catch (ExceptionObject & err)
	  {
	    cerr << "ExceptionObject caught!" << endl;
	    cerr << err << endl;
	    return EXIT_FAILURE;
	  }
	inputImage = addFilter->GetOutput();
      }
    IteratorType iterImage1 (inputImage, inputImage->GetBufferedRegion());
    while ( !iterImage1.IsAtEnd() )
      {
	PixelType NewValue = iterImage1.Get() / (NbFiles+1);
	iterImage1.Set(NewValue);
	++iterImage1;
      }    
    
   } else if(MajorityVotingOn) {

    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_majVoting");

    if(debug) cout << "Majority voting " << endl;

    // Reading label images
    vector<ImagePointer> vLabelImages;
    for (int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      {
	VolumeReaderType::Pointer LabelImageReader = VolumeReaderType::New();
	if (debug) cout << "Loading file " << InputFiles[LabelFileNumber] << endl;
	LabelImageReader->SetFileName(InputFiles[LabelFileNumber].c_str());
	try 
	  {
	    LabelImageReader->Update();
	  }
	catch (ExceptionObject & err) 
	  {
	    cerr<<"ExceptionObject caught!"<<endl;
	    cerr<<err<<endl;
	    return EXIT_FAILURE;	
	  }
	vLabelImages.push_back(LabelImageReader->GetOutput());
      }
    
    // Creating output image
    ImagePointer MajVotingImage = ImageType::New();
    MajVotingImage->SetRegions(inputImage->GetRequestedRegion());
    MajVotingImage->Allocate();
    
    //  Iterators initialization
    vector<ConstIteratorType> vConstLabelIterator;
    ConstIteratorType ConstInputIterator(inputImage,inputImage->GetRequestedRegion());
    IteratorType OutputIterator(MajVotingImage, inputImage->GetRequestedRegion());
    
    ConstInputIterator.GoToBegin();
    OutputIterator.GoToBegin();      
    for (int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
      {
	ConstIteratorType ConstLabelIterator(vLabelImages[LabelFileNumber],inputImage->GetRequestedRegion());
	vConstLabelIterator.push_back(ConstLabelIterator);
	vConstLabelIterator[LabelFileNumber].GoToBegin();	  
      }
    
    //  Compute the maximum value of intensity of the first parcellation image
    ShortPixelType MaxLabel;
    MaxFilterType::Pointer maxFilter = MaxFilterType::New();
    maxFilter->SetImage(vLabelImages[0]);
    maxFilter->ComputeMaximum();
    MaxLabel = (ShortPixelType) maxFilter->GetMaximum();
    MaxLabel++;

    //  Filling output image
    while (!ConstInputIterator.IsAtEnd())
      {
	PixelType MaxVoxelValue = 0, MaxLabelValue = 0;
	PixelType *LabelArray;
	LabelArray = new PixelType[MaxLabel];
	for (int Label = 0; Label < MaxLabel; Label++)
	  LabelArray[Label] = 0;
	
	for (int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
	  LabelArray[(ShortPixelType)vConstLabelIterator[LabelFileNumber].Get()]++;
	
	for (int Label = 1; Label < MaxLabel; Label++)
	  {
	    if (LabelArray[Label] > MaxVoxelValue)
	      {
		MaxVoxelValue = LabelArray[Label];
		MaxLabelValue = Label;
	      }
	  }
	
	OutputIterator.Set(MaxLabelValue);
	delete[] LabelArray;	  
	
	++ConstInputIterator;
	++OutputIterator;
	for (int LabelFileNumber = 0; LabelFileNumber < NbFiles; LabelFileNumber++)
	  ++vConstLabelIterator[LabelFileNumber];
      }

    //  Relabeling: considering neighborhood
    NeighborhoodIteratorType::RadiusType Radius;
      Radius.Fill(1);
      NeighborhoodIteratorType NeighborhoodOutputIterator(Radius,MajVotingImage,MajVotingImage->GetRequestedRegion());
      NeighborhoodIteratorType::OffsetType offset1 = {{-1,0,0}};
      NeighborhoodIteratorType::OffsetType offset2 = {{1,0,0}};
      NeighborhoodIteratorType::OffsetType offset3 = {{0,-1,0}};
      NeighborhoodIteratorType::OffsetType offset4 = {{0,1,0 }};
      NeighborhoodIteratorType::OffsetType offset5 = {{0,0,-1}};
      NeighborhoodIteratorType::OffsetType offset6 = {{0,0,1}};

      for (OutputIterator.GoToBegin(), NeighborhoodOutputIterator.GoToBegin(); !OutputIterator.IsAtEnd(); ++OutputIterator, ++NeighborhoodOutputIterator)
	{
	  PixelType MaxVoxelValue = 0, MaxLabelValue = 0;
	  PixelType *LabelArray;
	  LabelArray = new PixelType[MaxLabel];
	  for (int Label = 0; Label < MaxLabel; Label++)
	    LabelArray[Label] = 0;
	  LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset1)]++;
	  LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset2)]++;
	  LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset3)]++;
	  LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset4)]++;
	  LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset5)]++;
	  LabelArray[(ShortPixelType)NeighborhoodOutputIterator.GetPixel(offset6)]++; 	  

	  for (int Label = 0; Label < MaxLabel; Label++)
	    {
	      if (LabelArray[Label] > MaxVoxelValue)
		{
		  MaxVoxelValue = LabelArray[Label];
		  MaxLabelValue = Label;
		}
	    }

	  if ( (MaxVoxelValue >= 4) && (MaxLabelValue != 0) && (OutputIterator.Get() != MaxLabelValue) )
	    OutputIterator.Set(MaxLabelValue);

	  delete[] LabelArray;
	}
      
    inputImage = MajVotingImage;
  } else if( flip )
  { 
    ImageType::SizeType size ;
    size = inputImage->GetLargestPossibleRegion().GetSize() ;
    ImageType::PointType origin ;
    origin = inputImage->GetOrigin() ;
    itk::Index< 3 > index ;
    for( int i = 0 ; i < 3 ; i++ )
    {
      index[ i ] = size[ i ] - 1 ;
    }
    itk::Point< double , 3 > corner ;
    inputImage->TransformIndexToPhysicalPoint( index , corner ) ;
    itk::Point< double , 3 > newOrigin ;
    for( int i = 0 ; i < 3 ; i++ )
    {
     if( (int)flipMatrix[ i ][ i ] == -1 )
     {
       newOrigin[ i ] = -corner[ i ] ;
     }
     else
     {
       newOrigin[ i ] = origin[ i ] ;
     }
    }
    typedef itk::AffineTransform< double , 3 > AffineTransformType ;
    AffineTransformType::Pointer flipTransform = AffineTransformType::New() ;
    flipTransform->SetMatrix( flipMatrix ) ;
    //Resample Image
    typedef itk::NearestNeighborInterpolateImageFunction< ImageType , double > NearestNeighborInterpolateType ;
    NearestNeighborInterpolateType::Pointer interpolator = NearestNeighborInterpolateType::New() ;
    typedef itk::ResampleImageFilter< ImageType, ImageType > ResamplerType ;
    itk::ResampleImageFilter< ImageType , ImageType >::Pointer resampler ;
    resampler = itk::ResampleImageFilter< ImageType , ImageType >::New() ;
    resampler->SetOutputParametersFromImage( inputImage ) ;
    resampler->SetOutputOrigin( newOrigin ) ;
    resampler->SetInput( inputImage ) ;
    resampler->SetInterpolator( interpolator ) ;
    resampler->SetTransform( flipTransform ) ;
    resampler->Update() ;
    inputImage = resampler->GetOutput() ;
  }
  else {
    cout << "NOTHING TO DO, no operation selected..." << endl;
    exit(1);
  }

  // add the extension to the outputfile name
  outFileName.append(format);

  // when outputFileName is set
  if (outputFileName) {
     outFileName = string(outputFileName);
  }
  
  // write image
  if (debug) cout << "writing output data " << outFileName << endl;

if( diffusionImage )
{
    typedef itk::ImageFileWriter< DiffusionImageType > DiffusionImageWriter ;
    DiffusionImageType::Pointer outputDiffusionImage ;
    outputDiffusionImage = dynamic_cast< DiffusionImageType* >( inputBaseImage.GetPointer() ) ;
    if( !outputDiffusionImage ) 
    {
       std::cerr << "Error saving output diffusion image" << std::endl ;
       return EXIT_FAILURE ;
    }
    DiffusionImageWriter::Pointer writer = DiffusionImageWriter::New() ;
    if( !nocompOn )
    {
      writer->UseCompressionOn() ;
    }
    writer->SetFileName( outFileName.c_str() ); 
    writer->SetInput( outputDiffusionImage ) ;
    try
    {
      writer->Write() ;
    }
    catch( ExceptionObject & err )
    {
      cerr << "ExceptionObject caught!" << endl ;
      cerr << err << endl ;
      return EXIT_FAILURE ;
    }

}
else
{
  if (writeByte){
    castBinaryFilterType::Pointer castFilter = castBinaryFilterType::New();
    castFilter->SetInput(inputImage);
    try {
      castFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }
    
    BinaryVolumeWriterType::Pointer writer = BinaryVolumeWriterType::New();
    if(!nocompOn) writer->UseCompressionOn();
    writer->SetFileName(outFileName.c_str()); 
    writer->SetInput(castFilter->GetOutput());
    writer->Write();
  } else if (writeFloat) {
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    if(!nocompOn) writer->UseCompressionOn();
    writer->SetFileName(outFileName.c_str()); 
    writer->SetInput(inputImage);
    writer->Write();
  } else {
    castShortFilterType::Pointer castFilter = castShortFilterType::New();
    castFilter->SetInput(inputImage);
    try {
      castFilter->Update();
    }
    catch (ExceptionObject & err) {
      cerr << "ExceptionObject caught!" << endl;
      cerr << err << endl;
      return EXIT_FAILURE;	
    }
    
    ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
    writer->SetFileName(outFileName.c_str()); 
    if(!nocompOn) writer->UseCompressionOn();
    writer->SetInput(castFilter->GetOutput());
    writer->Write();
  }
}
  return 0;
}
