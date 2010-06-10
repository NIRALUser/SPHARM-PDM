/* 
 * compute image math and combinations
 *
 * author:  Martin Styner 
 * 
 * changes:
 *
 */



#include <iostream>
#include <fstream>
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
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkCastImageFilter.h>
#include <itkMaskImageFilter.h>

#include "argio.h"
#include "ImageMath.h" 

#define DEFAULT_SAMP 2
// number of samples by default

using namespace std;
using namespace itk;

typedef unsigned char PixelType;
typedef short ShortPixelType;
enum { ImageDimension = 3 };
typedef Image<PixelType,ImageDimension>       ImageType;
typedef Image<ShortPixelType,ImageDimension>  ShortImageType;
typedef ImageType::RegionType                 ImageRegionType;
typedef ImageRegionIterator< ImageType >      Iterator;
typedef ImageType::Pointer                    ImagePointer;

typedef ImageFileReader< ImageType >          VolumeReaderType;
typedef ImageFileWriter< ImageType >          VolumeWriterType;
typedef ImageFileWriter< ShortImageType >     ShortVolumeWriterType;

typedef CastImageFilter< ImageType,  ShortImageType > castShortFilterType; 

typedef BinaryThresholdImageFilter< ImageType , ImageType > threshFilterType;
typedef MaskImageFilter< ImageType, ImageType, ImageType >  maskFilterType;


static int debug = 0;

int main(const int argc, const char **argv)
{
  if (argc <=1 || ipExistsArgument(argv, "-usage") || ipExistsArgument(argv, "-help")) {
    cout << "ImageMath 1.0 version (Oct 2004)" << endl;
    cout << " computes base statistics of an image" << endl;
    cout << "usage: ImageStat infile [-outbase outbase] [-extractLabel <label>] [-combine infile2 [-relabel]] [-type byte|short] [-threshold tmin,tmax] [-mask infile2] [-constOper opID,val] [-extension ext] [-v]" << endl;
    cout << endl;
    cout << "infile         input dataset" << endl;;
    cout << "-outbase outbase     base-outputfilename, if omitted then the same as base of input" << endl; 
    cout << "-outfile outfile     outfilename, will only be applied to main output (if multiple output are given)" << endl;
    cout << "-combine infile2     combine the inputfile by interpreting them as labelfiles. " << endl 
       << "   -relabel          labels in infile2 will be relabeled to succeed labels in infile (no label overlap)" << endl;
    cout << "        labels in infile2 overwrite only the background label in infile1" << endl; 
    cout << "-extractLabel label  extract the mentioned label from the file" << endl;
    cout << "-threshold min,max   threshold: everything I < min , I > max is set to 0, otherwise 1" << endl;
    cout << "-mask infile2        use infile2 as a binary mask (intensities > 0 ) and combine with inputfile  " << endl ;
    cout << "-constOper opID,val       apply the following operation to the image: I op val, op = +/0, -/1, */2, //3" << endl;
    cout << "-v                   verbose mode " << endl;
    cout << "-type byte|short    Type of processing image (computations are always done with unsigned char images), default is short" << endl; 
    cout << "-extension ext       Extension of the output (determines output format)" << endl; 
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
  bool writeByte = false;
  if (typeChat && !strcmp(typeChat,"byte")) writeByte = true;

  char * formatChar = ipGetStringArgument(argv, "-extension", ".gipl");
  string format;
  if (! strchr(formatChar, '.')) {
    format = string(".") + string(formatChar);
  } else {
    format = string(formatChar);
  }

  debug      = ipExistsArgument(argv, "-v");

  char *combineFile    = ipGetStringArgument(argv, "-combine", NULL);  
  bool relabelOn = ipExistsArgument(argv, "-relabel");

  bool extractLabelOn   = ipExistsArgument(argv, "-extractLabel"); 
  int extractLabel   = ipGetIntArgument(argv, "-extractLabel", 0); 

  bool maskOn   = ipExistsArgument(argv, "-mask"); 
  char *maskFile    = ipGetStringArgument(argv, "-mask", NULL);  

  bool thresholdOn    = ipExistsArgument(argv, "-threshold"); 
  char * tmp_str      = ipGetStringArgument(argv, "-threshold", NULL);
  float textend[2];
  PixelType tmin = 0;
  PixelType tmax = 0;
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

  

  ImagePointer inputImage ;

  // load image
  if (debug) cout << "Loading file " << inputFileName << endl;
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(inputFileName) ;
  imageReader->Update();
  inputImage = imageReader->GetOutput();

  // do something to InputImage
  if (combineFile) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_comb");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << combineFile << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(combineFile) ;
    imageReader->Update();
    inputImage2 = imageReader->GetOutput();
    
    Iterator iterImage1 (inputImage, inputImage->GetBufferedRegion());
    Iterator iterImage2 (inputImage2, inputImage->GetBufferedRegion());
    
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
    threshFilter->Update();
    inputImage = threshFilter->GetOutput();
    
  }  else if (thresholdOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_thresh");
    if (debug) cout << "thresholding image  " << (int) tmin << "," <<(int) tmax << endl; 
    
    threshFilterType::Pointer threshFilter = threshFilterType::New();
    threshFilter->SetInput(inputImage);
    threshFilter->SetLowerThreshold(tmin);
    threshFilter->SetUpperThreshold(tmax);
    threshFilter->SetOutsideValue (BGVAL);
    threshFilter->SetInsideValue (FGVAL);
    threshFilter->Update();
    inputImage = threshFilter->GetOutput();
  } else if (maskOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_mask");
    
    ImagePointer inputImage2 ;
    if (debug) cout << "Loading file2 " << maskFile  << endl;
    VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
    imageReader->SetFileName(maskFile) ;
    imageReader->Update();
    inputImage2 = imageReader->GetOutput();
    if (debug) cout << "masking images  " << endl;

    maskFilterType::Pointer maskFilter = maskFilterType::New();
    maskFilter->SetInput1(inputImage);
    maskFilter->SetInput2(inputImage2);
    maskFilter->Update();
    inputImage = maskFilter->GetOutput();
    
  } else if (constOperOn) {
    outFileName.erase();
    outFileName.append(base_string);
    outFileName.append("_oper");
    
    Iterator iterImage1 (inputImage, inputImage->GetBufferedRegion());
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
  } else {
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
  if (writeByte){
    VolumeWriterType::Pointer writer = VolumeWriterType::New();
    writer->SetFileName(outFileName.c_str()); 
    writer->SetInput(inputImage);
    writer->Write();
  } else {
    castShortFilterType::Pointer castFilter = castShortFilterType::New();
    castFilter->SetInput(inputImage);
    castFilter->Update();

    ShortVolumeWriterType::Pointer writer = ShortVolumeWriterType::New();
    writer->SetFileName(outFileName.c_str()); 
    writer->SetInput(castFilter->GetOutput());
    writer->Write();
  }
  return 0;
}
